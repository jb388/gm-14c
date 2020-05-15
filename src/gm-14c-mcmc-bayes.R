# MCMC script for Bayesian parameter optimization
# aut: J. Beem-Miller
# date: 15-May-2020

# script requires running all code chunks prior to "Bayesian parameter estimation" section of "gm-14c.Rmd"

# set date for saving files
date <- Sys.Date()

## Markov Chain Monte Carlo parameter optimization
# Note that the model uses a bayesian prior from the Nelder-Mead optimization

# Model iterations and how many to exclude (burn-in)
iter <- 15000 # started with 5000, not all pars stabilized; still having stabilization issues after 15000
burnin <- 5000 # set to ~1/3

# Sohi data
# same upper/lower limits as in modFit
start.mcmc.s <- Sys.time()
bayes_fit_s_3p <- modMCMC(f = s.mod.Cost, 
                          p = s.mod.fit$par, 
                          var0 = s.mod.fit$var_ms,
                          upper = c(1, .5, .1, .8, .1), 
                          lower = c(0, 0, 0, 0, 0),
                          niter = iter, 
                          burninlength = burnin)
end.mcmc.s <- Sys.time() 
end.mcmc.s - start.mcmc.s # 31 min w/ 5000 iterations

# bayes_fit_s_3p$bestpar
# s.mod.fit$par
#####

# Zimmerman
# same upper/lower limits as in modFit
start.mcmc.z <- Sys.time()
bayes_fit_z_3p <- modMCMC(f = z.mod.Cost, 
                          p = z.mod.fit$par, 
                          var0 = z.mod.fit$var_ms,
                          upper = c(5, .1, .1, .5, .1),
                          lower = c(0, 0, 0, 0, 0),
                          iter = niter, 
                          burninlength = burnin)
end.mcmc.z <- Sys.time() 
end.mcmc.z - start.mcmc.z # 24.4 mins w/ 15000 iter

# bayes_fit_z_3p$bestpar
# z.mod.fit$par

# save bayesian model output for posterity 
# *WARNING* creates new directory by date, but will overwrite if files from current date exist!
save(bayes_fit_s_3p, file = paste0("gm-14c/data/derived/bayes-par-fit-", date, "/bayes_fit_s_3p-", iter, "iter", ".RData"))
save(bayes_fit_z_3p, file = paste0("gm-14c/data/derived/bayes-par-fit-", date, "/bayes_fit_s_3p-", iter, "iter", ".RData"))


## SA and TT uncertainty

# Function to calculate system age, pool ages, and transit time for all bayesian parameter combinations
sa.tt.fx <- function(pars, iter, a31 = FALSE) {
  
  # initialize list
  ls.nms <- c("SA.ls", "TT.ls", "fast.age.ls", "intm.age.ls", "slow.age.ls", "stock.ls")
  SA.TT.ls <- lapply(ls.nms, function(ls) {
    ls <- vector(mode = "list", length = iter)
  })
  names(SA.TT.ls) <- ls.nms
  
  # set progress bar
  pb <- txtProgressBar(min = 0, max = iter, style = 3)
  
  for (i in 1:iter) {
    
    # model matrix
    A <- -1 * diag(pars$par[i, 1:3])
    A[2, 1] <- pars$par[i, 4]
    if (a31) {
      A[3, 1] <- pars$par[i, 5] 
    } else {
      A[3, 2] <- pars$par[i, 5]
    }
    
    # calculate stocks
    stock <- sum(-1 * solve(A) %*% c(ws.in, ws.in * pars$par[i, 4], ws.in * pars$par[i, 5]))
    
    # System ages and transit times
    SA <- systemAge(A = A, u = In, a = ages)
    TT <- transitTime(A = A, u = In, a = ages)
    
    # Append to list
    SA.TT.ls[["SA.ls"]][[i]] <- as.numeric(SA$meanSystemAge)
    SA.TT.ls[["TT.ls"]][[i]] <- as.numeric(TT$meanTransitTime)
    SA.TT.ls[["fast.age.ls"]][[i]] <- as.numeric(SA$meanPoolAge[1])
    SA.TT.ls[["intm.age.ls"]][[i]] <- as.numeric(SA$meanPoolAge[2])
    SA.TT.ls[["slow.age.ls"]][[i]] <- as.numeric(SA$meanPoolAge[3])
    SA.TT.ls[["stock.ls"]][[i]] <- as.numeric(stock)
    
    # tracker
    setTxtProgressBar(pb, i)
  }
  return(SA.TT.ls)
}

# note that the two following function calls are very time consuming...
s.3p.SA.TT.ls <- sa.tt.fx(pars = bayes_fit_s_3p, iter = nrow(bayes_fit_s_3p$pars), a31 = TRUE)
z.3p.SA.TT.ls <- sa.tt.fx(pars = bayes_fit_z_3p, iter = nrow(bayes_fit_z_3p$pars), a31 = FALSE)

# save results
save(s.3p.SA.TT.ls, file = paste0("gm-14c/data/derived/bayes-par-fit-", date, "/bayes.s.3p.SA.TT.", iter, "iter", ".RData"))
save(z.3p.SA.TT.ls, file = paste0("gm-14c/data/derived/bayes-par-fit-", date, "/bayes.z.3p.SA.TT.", iter, "iter", ".RData"))
