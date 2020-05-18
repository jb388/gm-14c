# MCMC script for Bayesian parameter optimization
# aut: J. Beem-Miller
# date: 15-May-2020

# script requires running all code chunks prior to "Bayesian parameter estimation" section of "gm-14c.Rmd"

# set date for saving files
date <- Sys.Date()

## Markov Chain Monte Carlo parameter optimization
# Note that the model uses a bayesian prior from the Nelder-Mead optimization

# Model iterations and how many to exclude (burn-in)
iter <- 15000 # started with 5000, not all pars stabilized;
updatecov <- 10 # good par stabilization with updatecov set to 10!

# for saving script output
save.iter <- paste0(iter, "iter", ".RData")
save.dir <- file.path(paste0("gm-14c/data/derived/bayes-par-fit-", date))
if(!dir.exists(save.dir)) {dir.create(file.path(paste0("gm-14c/data/derived/bayes-par-fit-", date)))}  

# Sohi data
# same upper/lower limits as in modFit
start.mcmc.s <- Sys.time()
bayes_fit_s_3p <- modMCMC(f = s.mod.Cost, 
                          p = s.mod.fit$par, 
                          var0 = s.mod.fit$var_ms,
                          upper = c(1, .1, .1, .5, .1), 
                          lower = c(0, 0, 0, 0, 0),
                          niter = iter, 
                          burninlength = 0,
                          updatecov = updatecov)
end.mcmc.s <- Sys.time() 
end.mcmc.s - start.mcmc.s # 31 min w/ 5000 iterations

# bayes_fit_s_3p$bestpar
# s.mod.fit$par

# Zimmerman
# same upper/lower limits as in modFit
start.mcmc.z <- Sys.time()
bayes_fit_z_3p <- modMCMC(f = z.mod.Cost, 
                          p = z.mod.fit$par, 
                          var0 = z.mod.fit$var_ms,
                          upper = c(1, .1, .1, .5, .1),
                          lower = c(0, 0, 0, 0, 0),
                          niter = iter, 
                          burninlength = 0,
                          updatecov = 10)
end.mcmc.z <- Sys.time() 
end.mcmc.z - start.mcmc.z # 24.4 mins w/ 15000 iter

# save output
# *WARNING* will overwrite if files from current date exist!
save(bayes_fit_s_3p, file = paste0(save.dir, "/bayes_fit_s_3p-", save.iter))
save(bayes_fit_z_3p, file = paste0(save.dir, "/bayes_fit_z_3p-", save.iter))


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

# average par est by every 100 to speed up following calls
brks <- seq(1, nrow(bayes_fit_s_3p$pars), by = 100)
# Sohi
bayes_fit_s_3p.100avg <- lapply(brks, function(x) matrix(nrow = 1, ncol = ncol(bayes_fit_s_3p$pars)))
for(i in seq_along(brks)) {
  if(brks[i] < brks[length(brks)]) {
  bayes_fit_s_3p.100avg[[i]] <- apply(bayes_fit_s_3p$pars[brks[i]:brks[i+1],], 2, mean) 
  } else {
    bayes_fit_s_3p.100avg[[i]] <- apply(bayes_fit_s_3p$pars[brks[i]:nrow(bayes_fit_s_3p$pars),], 2, mean)  
  }
}
bayes_fit_s_3p.100avg <- do.call(rbind, bayes_fit_s_3p.100avg)
bayes_fit_s_3p.ls <- list(pars = bayes_fit_s_3p.100avg)

# Zimmerman
bayes_fit_z_3p.100avg <- lapply(brks, function(x) matrix(nrow = 1, ncol = ncol(bayes_fit_z_3p$pars)))
for(i in seq_along(brks)) {
  if(brks[i] < brks[length(brks)]) {
    bayes_fit_z_3p.100avg[[i]] <- apply(bayes_fit_z_3p$pars[brks[i]:brks[i+1],], 2, mean) 
  } else {
    bayes_fit_z_3p.100avg[[i]] <- apply(bayes_fit_z_3p$pars[brks[i]:nrow(bayes_fit_z_3p$pars),], 2, mean)  
  }
}
bayes_fit_z_3p.100avg <- do.call(rbind, bayes_fit_z_3p.100avg)
bayes_fit_z_3p.ls <- list(pars = bayes_fit_z_3p.100avg)

# note that the two following function calls are very time consuming...
s.3p.SA.TT.ls <- sa.tt.fx(pars = bayes_fit_s_3p.ls, iter = nrow(bayes_fit_s_3p.ls$pars), a31 = TRUE)
z.3p.SA.TT.ls <- sa.tt.fx(pars = bayes_fit_z_3p.ls, iter = nrow(bayes_fit_z_3p.ls$pars), a31 = FALSE)

# save output
# *WARNING* will overwrite if files from current date exist!
save(s.3p.SA.TT.ls, file = paste0(save.dir, "/bayes_fit_s_3p.SA.TT.", save.iter))
save(z.3p.SA.TT.ls, file = paste0(save.dir, "/bayes_fit_z_3p.SA.TT.", save.iter))
