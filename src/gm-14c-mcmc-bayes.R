# MCMC script for Bayesian parameter optimization
# aut: J. Beem-Miller
# date: 14-May-2020

# script requires running all code chunks prior to "Bayesian parameter estimation" section of "gm-14c.Rmd"

# Model iterations and how many to exclude (burn-in)
niter <- 5000 # started with 5000, not all pars stabilized
burnin <- 500 # set to ~1/3

# Sohi data
#####
# run monte carlo estimatation with optimized parameters
start <- Sys.time()
bayes_fit_s_3p <- modMCMC(f = s.mod.Cost, 
                          p = s.mod.fit$par, 
                          var0 = s.mod.fit$var_ms,
                          upper = c(1, .5, .1, .8, .1), 
                          lower = c(0, 0, 0, 0, 0),
                          niter = niter, 
                          burninlength = burnin)
end <- Sys.time() 
end - start # 20.6 min w/ 5000 iterations

bayes_fit_s_3p$bestpar
s.mod.fit$par

bayesFacts <- data.frame(niter, burnin, bayes_fit_s_3p$naccepted)
colnames(bayesFacts) <- c('# Iterations', "# Burn-In", "# Accepted")

# Check that parameters have stabilized after burnin (only kfPOM and kOmin)
plot(bayes_fit_s_3p)

# Estimate parameter sensitivity and return timeseries distribution envelope
pred_uncert_s_3p <- sensRange(s.modFun_3p, parInput = bayes_fit_s_3p$pars)
sens_s_3p <- summary(pred_uncert_s_3p)

# # summary
# round(summary(bayes_fit_s_3p), 4)

# fit mod w/ bayes pars
s.F0_Delta14C.bayes <- unlist(lapply(bayes_fit_s_3p$bestpar[1:3], function(x) fm_14c(fm(x), 1883)))
s.3pc.bayes <- sohi.3p.mod(
  t = yrs,
  ks = bayes_fit_s_3p$bestpar[1:3],
  C0 = s.stock.3p.c,
  F0_Delta14C = s.F0_Delta14C.bayes,
  In = ws.in,
  a21 = bayes_fit_s_3p$bestpar[4],
  a31 = bayes_fit_s_3p$bestpar[5],
  inputFc = Datm
)

s.3p.C14m.bayes <- getF14C(s.3pc.bayes) 
s.3p.C14.bayes <- getF14(s.3pc.bayes)
s.3p.Ctot.bayes <- getC(s.3pc.bayes)

s.3ps.C14.df.bayes <- data.frame(
  years = rep(Datm$Date, 5),
  d14C = c(s.3p.C14.bayes[, 1], s.3p.C14.bayes[, 2], s.3p.C14.bayes[, 3], s.3p.C14m.bayes, Datm$NHc14),
  pool = rep(c("fPOM", "oPOM", "Omin", "total C", "atm"), each = nrow(s.3p.C14.bayes))
)
#####

# Zimmerman
#####

# run monte carlo estimatation with optimized parameters
# same upper/lower limits as in modFit
start.z <- Sys.time()
bayes_fit_z_3p <- modMCMC(f = z.mod.Cost, 
                          p = z.mod.fit$par, 
                          var0 = z.mod.fit$var_ms,
                          upper = c(5, .1, .1, .5, .1),
                          lower = c(0, 0, 0, 0, 0),
                          niter = niter, 
                          burninlength = burnin)
end.z <- Sys.time() 
end.z - start.z # 20.4 mins w/ 5000 iter

bayes_fit_z_3p$bestpar
z.mod.fit$par

bayesFacts.z <- data.frame(niter, burnin, bayes_fit_z_3p$naccepted)
colnames(bayesFacts.z) <- c('# Iterations', "# Burn-In", "# Accepted")

# Check that parameters have stabilized after burnin (only kfPOM and kOmin)
plot(bayes_fit_z_3p)

# # summary
# round(summary(bayes_fit_s_3p), 4)

# Estimate parameter sensitivity and return timeseries distribution envelope
pred_uncert_z_3p <- sensRange(z.modFun_3p, parInput = bayes_fit_z_3p$pars)
sens_z_3p <- summary(pred_uncert_z_3p)

# fit mod w/ bayes pars
z.F0_Delta14C.bayes <- unlist(lapply(bayes_fit_z_3p$bestpar[1:3], function(x) fm_14c(fm(x), 1883)))
z.3ps.bayes <- ThreepSeriesModel14(
  t = yrs,
  ks = bayes_fit_z_3p$bestpar[1:3],
  C0 = z.stock.3p.c,
  F0_Delta14C = z.F0_Delta14C.bayes,
  In = ws.in,
  a21 = bayes_fit_z_3p$bestpar[4],
  a32 = bayes_fit_z_3p$bestpar[5],
  inputFc = Datm
)

z.3p.C14m.bayes <- getF14C(z.3ps.bayes) 
z.3p.C14.bayes <- getF14(z.3ps.bayes)
z.3p.Ctot.bayes <- getC(z.3ps.bayes)

z.3ps.C14.df.bayes <- data.frame(
  years = rep(Datm$Date, 5),
  d14C = c(z.3p.C14.bayes[, 1], z.3p.C14.bayes[, 2], z.3p.C14.bayes[, 3], z.3p.C14m.bayes, Datm$NHc14),
  pool = rep(c("fast", "intm", "slow", "total C", "atm"), each = nrow(z.3p.C14.bayes))
)