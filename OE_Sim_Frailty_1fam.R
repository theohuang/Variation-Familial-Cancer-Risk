## O/E simulations for 1 family, not 4-family aggregates
## Last updated: July 6, 2020

library(dplyr)
library(abind)
library(ggplot2)
library(tidyverse)

load("penet.mmr.net.RData")
load("death.othercauses.RData")

source("OE Functions.R")
source("Estimating Functions Discrete.R")

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source(paste0(getwd(), "/Generating Families Functions/sim.simFam.R"))
source(paste0(getwd(), "/Generating Families Functions/genCancerPen.R"))
source(paste0(getwd(), "/Generating Families Functions/genOtherDeathPen.R"))
source(paste0(getwd(), "/Generating Families Functions/sim.buildGenoMat.R"))
source(paste0(getwd(), "/Generating Families Functions/sim.linkParents.R"))
source(paste0(getwd(), "/Generating Families Functions/sim.simCurAgeVar.R"))
source(paste0(getwd(), "/Generating Families Functions/sim.simCancerVars.R"))
source(paste0(getwd(), "/Generating Families Functions/sim.buildBranchOfAlleleMats.R"))
source(paste0(getwd(), "/Generating Families Functions/helpers.R"))

nfam <- 53
mutations <- c("MLH1", "MSH2", "MSH6")
cancers <- c("ColorC", "EndomC")
genos <- c("M000", "M100", "M010", "M001")

pen2hzd <- function(pen){
  apply(pen, 2, function(x) x / c(1, 1 - cumsum(x)[-length(x)]))
}

hzd2pen <- function(hzd){
  surv <- apply(1 - hzd, 2, cumprod)
  pen <- rbind(hzd[1, ], surv[-nrow(hzd), ] * hzd[-1, ])
  return(pen)
}

hzd02hzd <- function(hzd0, w){
  1 - (1 - hzd0)^(exp(w))
}


## getting baseline hazards -- assuming they are the MMRpro hazards
pen0F <- list(ColorC = penet.mmr.net$fFX[, genos],
              EndomC = penet.mmr.net$fFY[, genos])
pen0M <- list(ColorC = penet.mmr.net$fMX[, genos],
              EndomC = penet.mmr.net$fMY[, genos])
hzd0F <- lapply(pen0F, pen2hzd)
hzd0M <- lapply(pen0M, pen2hzd)


## using allele frequencies of 0.1 to increase the number of carriers
af <- setNames(rep(0.1, 3), mutations)
## baseline penetrance (MMRpro)
CP0 <- genCancerPen(mutations, cancers, pen0F, pen0M, maxK = length(mutations), age.last = 95)


ODP <- genOtherDeathPen(cancers,
                        mutate(select(death.othercauses[1:94, ],
                                      femaleCRC, maleCRC, uterus),
                               femaleColorC = femaleCRC, maleColorC = maleCRC,
                               femaleEndomC = uterus, maleEndomC = uterus),
                        age.max = 94)

cnames <- c("FamID", "W", "O.C", "E.C", "OE.C", "O.NC", "E.NC",
            "OE.NC", "OEdiff.C", "OEdiff.NC", "FamSize")



#### Frailty on all family members #####

load("~/Dropbox (Partners HealthCare)/BayesMendel Dropbox folders/Variation Familial Cancer Risk/OE_Sim_CNC_Main.RData")


res.oe.cnc.1fam <- setNames(data.frame(matrix(0, nfam, length(cnames))), cnames)
res.oe.cnc.1fam$FamID <- 1:nfam
for(k in 1:nfam){
  print(k)
  fam <- fams[[k]][[1]]
  res.exp <- expect(fam, CP0$cancerMDens[1:94, , "ColorC"], CP0$cancerFDens[1:94, , "ColorC"])
  res.obs <- observe(fam, CP0$cancerMDens[1:94, , "ColorC"], CP0$cancerFDens[1:94, , "ColorC"])
  
  ### O/E ratios
  ## using carrier probabilities
  # with censoring
  LIK <- estLik(fam, CP0, ODP)
  probs <- setNames(data.frame(matrix(0, nrow(fam), ncol(CP0$cancerFDens[, , 1]))),
                    colnames(CP0$cancerFDens[, , 1]))
  for(i in 1:nrow(fam)){
    probs[i, ] <- pp.peelingParing(fam, af, LIK, length(mutations),
                                   counselee.id = fam$ID[i])
  }
  
  res.oe.cnc.1fam$FamSize[k] <- nrow(fam)
  res.oe.cnc.1fam$O.C[k] <- sum(res.obs[, -1] * probs[, -1])
  res.oe.cnc.1fam$E.C[k] <- sum(res.exp[, -1] * probs[, -1])
  res.oe.cnc.1fam$O.NC[k] <- sum(res.obs[, 1] * probs[, 1])
  res.oe.cnc.1fam$E.NC[k] <- sum(res.exp[, 1] * probs[, 1])
  res.oe.cnc.1fam$OE.C[k] <- res.oe.cnc.1fam$O.C[k] / res.oe.cnc.1fam$E.C[k]
  res.oe.cnc.1fam$OE.NC[k] <- res.oe.cnc.1fam$O.NC[k] / res.oe.cnc.1fam$E.NC[k]
  res.oe.cnc.1fam$OEdiff.C[k] <- (res.oe.cnc.1fam$O.C[k] - res.oe.cnc.1fam$E.C[k]) / res.oe.cnc.1fam$FamSize[k]
  res.oe.cnc.1fam$OEdiff.NC[k] <- (res.oe.cnc.1fam$O.NC[k] - res.oe.cnc.1fam$E.NC[k]) / res.oe.cnc.1fam$FamSize[k]
}


## bootstrap
nboot <- 100
res.oe.cnc.1fam.boot <- vector("list", nfam)

start <- Sys.time()
set.seed(1)
boot.w <- vector("list", nfam)
for(k in 1:nfam){
  print(k)
  boot.w <- mutate(oe.boot.sim(fams[[k]][[1]], CP0, ODP, af, mutations, nboot, seed = 999),
                   FamSize = nrow(fams[[k]][[1]]))
  
  boot.w <- mutate(boot.w, FamID = k, OE.C = Obs.C / Exp.C, OE.NC = Obs.NC / Exp.NC,
                   OEdiff.C = (Obs.C - Exp.C) / FamSize, OEdiff.NC = (Obs.NC - Exp.NC) / FamSize)
  res.oe.cnc.1fam.boot[[k]] <- boot.w
}
difftime(Sys.time(), start, units = "secs")


res.oe.cnc.1fam <- mutate(res.oe.cnc.1fam, OE.C.lo = NA, OE.C.hi = NA,
                          OE.NC.lo = NA, OE.NC.hi = NA,
                          OEdiff.C.lo = NA, OEdiff.C.hi = NA,
                          OEdiff.NC.lo = NA, OEdiff.NC.hi = NA)
for(i in 1:nfam){
  res.oe.cnc.1fam$OE.C.lo[i] <- quantile(res.oe.cnc.1fam.boot[[i]]$OE.C, 0.025)
  res.oe.cnc.1fam$OE.C.hi[i] <- quantile(res.oe.cnc.1fam.boot[[i]]$OE.C, 0.975)
  res.oe.cnc.1fam$OE.NC.lo[i] <- quantile(res.oe.cnc.1fam.boot[[i]]$OE.NC, 0.025)
  res.oe.cnc.1fam$OE.NC.hi[i] <- quantile(res.oe.cnc.1fam.boot[[i]]$OE.NC, 0.975)
  res.oe.cnc.1fam$OEdiff.C.lo[i] <- quantile(res.oe.cnc.1fam.boot[[i]]$OEdiff.C, 0.025)
  res.oe.cnc.1fam$OEdiff.C.hi[i] <- quantile(res.oe.cnc.1fam.boot[[i]]$OEdiff.C, 0.975)
  res.oe.cnc.1fam$OEdiff.NC.lo[i] <- quantile(res.oe.cnc.1fam.boot[[i]]$OEdiff.NC, 0.025)
  res.oe.cnc.1fam$OEdiff.NC.hi[i] <- quantile(res.oe.cnc.1fam.boot[[i]]$OEdiff.NC, 0.975)
}
res.oe.cnc.1fam$W <- res.oe.sim$W



#### Frailty on carriers only #####

load("~/Dropbox (Partners HealthCare)/BayesMendel Dropbox folders/Variation Familial Cancer Risk/OE_Sim_C_Main.RData")


res.oe.c.1fam <- setNames(data.frame(matrix(0, nfam, length(cnames))), cnames)
res.oe.c.1fam$FamID <- 1:nfam
for(k in 1:nfam){
  print(k)
  fam <- fams[[k]][[1]]
  res.exp <- expect(fam, CP0$cancerMDens[1:94, , "ColorC"], CP0$cancerFDens[1:94, , "ColorC"])
  res.obs <- observe(fam, CP0$cancerMDens[1:94, , "ColorC"], CP0$cancerFDens[1:94, , "ColorC"])
  
  ### O/E ratios
  ## using carrier probabilities
  # with censoring
  LIK <- estLik(fam, CP0, ODP)
  probs <- setNames(data.frame(matrix(0, nrow(fam), ncol(CP0$cancerFDens[, , 1]))),
                    colnames(CP0$cancerFDens[, , 1]))
  for(i in 1:nrow(fam)){
    probs[i, ] <- pp.peelingParing(fam, af, LIK, length(mutations),
                                   counselee.id = fam$ID[i])
  }
  
  res.oe.c.1fam$FamSize[k] <- nrow(fam)
  res.oe.c.1fam$O.C[k] <- sum(res.obs[, -1] * probs[, -1])
  res.oe.c.1fam$E.C[k] <- sum(res.exp[, -1] * probs[, -1])
  res.oe.c.1fam$O.NC[k] <- sum(res.obs[, 1] * probs[, 1])
  res.oe.c.1fam$E.NC[k] <- sum(res.exp[, 1] * probs[, 1])
  res.oe.c.1fam$OE.C[k] <- res.oe.c.1fam$O.C[k] / res.oe.c.1fam$E.C[k]
  res.oe.c.1fam$OE.NC[k] <- res.oe.c.1fam$O.NC[k] / res.oe.c.1fam$E.NC[k]
  res.oe.c.1fam$OEdiff.C[k] <- (res.oe.c.1fam$O.C[k] - res.oe.c.1fam$E.C[k]) / res.oe.c.1fam$FamSize[k]
  res.oe.c.1fam$OEdiff.NC[k] <- (res.oe.c.1fam$O.NC[k] - res.oe.c.1fam$E.NC[k]) / res.oe.c.1fam$FamSize[k]
}


## bootstrap
nboot <- 100
res.oe.c.1fam.boot <- vector("list", nfam)

start <- Sys.time()
set.seed(1)
boot.w <- vector("list", nfam)
for(k in 1:nfam){
  print(k)
  boot.w <- mutate(oe.boot.sim(fams[[k]][[1]], CP0, ODP, af, mutations, nboot, seed = 999),
                   FamSize = nrow(fams[[k]][[1]]))
  
  boot.w <- mutate(boot.w, FamID = k, OE.C = Obs.C / Exp.C, OE.NC = Obs.NC / Exp.NC,
                   OEdiff.C = (Obs.C - Exp.C) / FamSize, OEdiff.NC = (Obs.NC - Exp.NC) / FamSize)
  res.oe.c.1fam.boot[[k]] <- boot.w
}
difftime(Sys.time(), start, units = "secs")


res.oe.c.1fam <- mutate(res.oe.c.1fam, OE.C.lo = NA, OE.C.hi = NA,
                        OE.NC.lo = NA, OE.NC.hi = NA,
                        OEdiff.C.lo = NA, OEdiff.C.hi = NA,
                        OEdiff.NC.lo = NA, OEdiff.NC.hi = NA)
for(i in 1:nfam){
  res.oe.c.1fam$OE.C.lo[i] <- quantile(res.oe.c.1fam.boot[[i]]$OE.C, 0.025)
  res.oe.c.1fam$OE.C.hi[i] <- quantile(res.oe.c.1fam.boot[[i]]$OE.C, 0.975)
  res.oe.c.1fam$OE.NC.lo[i] <- quantile(res.oe.c.1fam.boot[[i]]$OE.NC, 0.025)
  res.oe.c.1fam$OE.NC.hi[i] <- quantile(res.oe.c.1fam.boot[[i]]$OE.NC, 0.975)
  res.oe.c.1fam$OEdiff.C.lo[i] <- quantile(res.oe.c.1fam.boot[[i]]$OEdiff.C, 0.025)
  res.oe.c.1fam$OEdiff.C.hi[i] <- quantile(res.oe.c.1fam.boot[[i]]$OEdiff.C, 0.975)
  res.oe.c.1fam$OEdiff.NC.lo[i] <- quantile(res.oe.c.1fam.boot[[i]]$OEdiff.NC, 0.025)
  res.oe.c.1fam$OEdiff.NC.hi[i] <- quantile(res.oe.c.1fam.boot[[i]]$OEdiff.NC, 0.975)
}
res.oe.c.1fam$W <- res.oe.sim$W


#### Frailty on noncarriers only #####

load("~/Dropbox (Partners HealthCare)/BayesMendel Dropbox folders/Variation Familial Cancer Risk/OE_Sim_NC_Main.RData")


res.oe.nc.1fam <- setNames(data.frame(matrix(0, nfam, length(cnames))), cnames)
res.oe.nc.1fam$FamID <- 1:nfam
for(k in 1:nfam){
  print(k)
  fam <- fams[[k]][[1]]
  res.exp <- expect(fam, CP0$cancerMDens[1:94, , "ColorC"], CP0$cancerFDens[1:94, , "ColorC"])
  res.obs <- observe(fam, CP0$cancerMDens[1:94, , "ColorC"], CP0$cancerFDens[1:94, , "ColorC"])
  
  ### O/E ratios
  ## using carrier probabilities
  # with censoring
  LIK <- estLik(fam, CP0, ODP)
  probs <- setNames(data.frame(matrix(0, nrow(fam), ncol(CP0$cancerFDens[, , 1]))),
                    colnames(CP0$cancerFDens[, , 1]))
  for(i in 1:nrow(fam)){
    probs[i, ] <- pp.peelingParing(fam, af, LIK, length(mutations),
                                   counselee.id = fam$ID[i])
  }
  
  res.oe.nc.1fam$FamSize[k] <- nrow(fam)
  res.oe.nc.1fam$O.C[k] <- sum(res.obs[, -1] * probs[, -1])
  res.oe.nc.1fam$E.C[k] <- sum(res.exp[, -1] * probs[, -1])
  res.oe.nc.1fam$O.NC[k] <- sum(res.obs[, 1] * probs[, 1])
  res.oe.nc.1fam$E.NC[k] <- sum(res.exp[, 1] * probs[, 1])
  res.oe.nc.1fam$OE.C[k] <- res.oe.nc.1fam$O.C[k] / res.oe.nc.1fam$E.C[k]
  res.oe.nc.1fam$OE.NC[k] <- res.oe.nc.1fam$O.NC[k] / res.oe.nc.1fam$E.NC[k]
  res.oe.nc.1fam$OEdiff.C[k] <- (res.oe.nc.1fam$O.C[k] - res.oe.nc.1fam$E.C[k]) / res.oe.nc.1fam$FamSize[k]
  res.oe.nc.1fam$OEdiff.NC[k] <- (res.oe.nc.1fam$O.NC[k] - res.oe.nc.1fam$E.NC[k]) / res.oe.nc.1fam$FamSize[k]
}


## bootstrap
nboot <- 100
res.oe.nc.1fam.boot <- vector("list", nfam)

start <- Sys.time()
set.seed(1)
boot.w <- vector("list", nfam)
for(k in 1:nfam){
  print(k)
  boot.w <- mutate(oe.boot.sim(fams[[k]][[1]], CP0, ODP, af, mutations, nboot, seed = 999),
                   FamSize = nrow(fams[[k]][[1]]))
  
  boot.w <- mutate(boot.w, FamID = k, OE.C = Obs.C / Exp.C, OE.NC = Obs.NC / Exp.NC,
                   OEdiff.C = (Obs.C - Exp.C) / FamSize, OEdiff.NC = (Obs.NC - Exp.NC) / FamSize)
  res.oe.nc.1fam.boot[[k]] <- boot.w
}
difftime(Sys.time(), start, units = "secs")


res.oe.nc.1fam <- mutate(res.oe.nc.1fam, OE.C.lo = NA, OE.C.hi = NA,
                         OE.NC.lo = NA, OE.NC.hi = NA,
                         OEdiff.C.lo = NA, OEdiff.C.hi = NA,
                         OEdiff.NC.lo = NA, OEdiff.NC.hi = NA)
for(i in 1:nfam){
  res.oe.nc.1fam$OE.C.lo[i] <- quantile(res.oe.nc.1fam.boot[[i]]$OE.C, 0.025)
  res.oe.nc.1fam$OE.C.hi[i] <- quantile(res.oe.nc.1fam.boot[[i]]$OE.C, 0.975)
  res.oe.nc.1fam$OE.NC.lo[i] <- quantile(res.oe.nc.1fam.boot[[i]]$OE.NC, 0.025)
  res.oe.nc.1fam$OE.NC.hi[i] <- quantile(res.oe.nc.1fam.boot[[i]]$OE.NC, 0.975)
  res.oe.nc.1fam$OEdiff.C.lo[i] <- quantile(res.oe.nc.1fam.boot[[i]]$OEdiff.C, 0.025)
  res.oe.nc.1fam$OEdiff.C.hi[i] <- quantile(res.oe.nc.1fam.boot[[i]]$OEdiff.C, 0.975)
  res.oe.nc.1fam$OEdiff.NC.lo[i] <- quantile(res.oe.nc.1fam.boot[[i]]$OEdiff.NC, 0.025)
  res.oe.nc.1fam$OEdiff.NC.hi[i] <- quantile(res.oe.nc.1fam.boot[[i]]$OEdiff.NC, 0.975)
}
res.oe.nc.1fam$W <- res.oe.sim$W


#### Frailty on none ####

load("~/Dropbox (Partners HealthCare)/BayesMendel Dropbox folders/Variation Familial Cancer Risk/OE_Sim_None_Main.RData")


res.oe.none.1fam <- setNames(data.frame(matrix(0, nfam, length(cnames))), cnames)
res.oe.none.1fam$FamID <- 1:nfam
for(k in 1:nfam){
  print(k)
  fam <- fams[[k]][[1]]
  res.exp <- expect(fam, CP0$cancerMDens[1:94, , "ColorC"], CP0$cancerFDens[1:94, , "ColorC"])
  res.obs <- observe(fam, CP0$cancerMDens[1:94, , "ColorC"], CP0$cancerFDens[1:94, , "ColorC"])
  
  ### O/E ratios
  ## using carrier probabilities
  # with censoring
  LIK <- estLik(fam, CP0, ODP)
  probs <- setNames(data.frame(matrix(0, nrow(fam), ncol(CP0$cancerFDens[, , 1]))),
                    colnames(CP0$cancerFDens[, , 1]))
  for(i in 1:nrow(fam)){
    probs[i, ] <- pp.peelingParing(fam, af, LIK, length(mutations),
                                   counselee.id = fam$ID[i])
  }
  
  res.oe.none.1fam$FamSize[k] <- nrow(fam)
  res.oe.none.1fam$O.C[k] <- sum(res.obs[, -1] * probs[, -1])
  res.oe.none.1fam$E.C[k] <- sum(res.exp[, -1] * probs[, -1])
  res.oe.none.1fam$O.NC[k] <- sum(res.obs[, 1] * probs[, 1])
  res.oe.none.1fam$E.NC[k] <- sum(res.exp[, 1] * probs[, 1])
  res.oe.none.1fam$OE.C[k] <- res.oe.none.1fam$O.C[k] / res.oe.none.1fam$E.C[k]
  res.oe.none.1fam$OE.NC[k] <- res.oe.none.1fam$O.NC[k] / res.oe.none.1fam$E.NC[k]
  res.oe.none.1fam$OEdiff.C[k] <- (res.oe.none.1fam$O.C[k] - res.oe.none.1fam$E.C[k]) / res.oe.none.1fam$FamSize[k]
  res.oe.none.1fam$OEdiff.NC[k] <- (res.oe.none.1fam$O.NC[k] - res.oe.none.1fam$E.NC[k]) / res.oe.none.1fam$FamSize[k]
}


## bootstrap
nboot <- 100
res.oe.none.1fam.boot <- vector("list", nfam)

start <- Sys.time()
set.seed(1)
boot.w <- vector("list", nfam)
for(k in 1:nfam){
  print(k)
  boot.w <- mutate(oe.boot.sim(fams[[k]][[1]], CP0, ODP, af, mutations, nboot, seed = 999),
                   FamSize = nrow(fams[[k]][[1]]))
  
  boot.w <- mutate(boot.w, FamID = k, OE.C = Obs.C / Exp.C, OE.NC = Obs.NC / Exp.NC,
                   OEdiff.C = (Obs.C - Exp.C) / FamSize, OEdiff.NC = (Obs.NC - Exp.NC) / FamSize)
  res.oe.none.1fam.boot[[k]] <- boot.w
}
difftime(Sys.time(), start, units = "secs")


res.oe.none.1fam <- mutate(res.oe.none.1fam, W = NA, OE.C.lo = NA, OE.C.hi = NA,
                           OE.NC.lo = NA, OE.NC.hi = NA,
                           OEdiff.C.lo = NA, OEdiff.C.hi = NA,
                           OEdiff.NC.lo = NA, OEdiff.NC.hi = NA)
for(i in 1:nfam){
  res.oe.none.1fam$W
  res.oe.none.1fam$OE.C.lo[i] <- quantile(res.oe.none.1fam.boot[[i]]$OE.C, 0.025)
  res.oe.none.1fam$OE.C.hi[i] <- quantile(res.oe.none.1fam.boot[[i]]$OE.C, 0.975)
  res.oe.none.1fam$OE.NC.lo[i] <- quantile(res.oe.none.1fam.boot[[i]]$OE.NC, 0.025)
  res.oe.none.1fam$OE.NC.hi[i] <- quantile(res.oe.none.1fam.boot[[i]]$OE.NC, 0.975)
  res.oe.none.1fam$OEdiff.C.lo[i] <- quantile(res.oe.none.1fam.boot[[i]]$OEdiff.C, 0.025)
  res.oe.none.1fam$OEdiff.C.hi[i] <- quantile(res.oe.none.1fam.boot[[i]]$OEdiff.C, 0.975)
  res.oe.none.1fam$OEdiff.NC.lo[i] <- quantile(res.oe.none.1fam.boot[[i]]$OEdiff.NC, 0.025)
  res.oe.none.1fam$OEdiff.NC.hi[i] <- quantile(res.oe.none.1fam.boot[[i]]$OEdiff.NC, 0.975)
}


####### ashr #######

library(ashr)

#### Frailty model on all ####

oe.cnc.c <- oe.cnc.nc <- s.oe.cnc.c <- s.oe.cnc.nc <- rep(NA, nfam)
for(i in 1:nfam){
  oe.cnc.c[i] <- res.oe.cnc.1fam$OE.C[i] - 1
  oe.cnc.nc[i] <- res.oe.cnc.1fam$OE.NC[i] - 1
  s.oe.cnc.c[i] <- sd(res.oe.cnc.1fam.boot[[i]]$OE.C)
  s.oe.cnc.nc[i] <- sd(res.oe.cnc.1fam.boot[[i]]$OE.NC)
}


##### half-uniform #####

res.sim.hu.c <- ash(oe.cnc.c, s.oe.cnc.c, mixcompdist = "halfuniform")
res.sim.hu.nc <- ash(oe.cnc.nc, s.oe.cnc.nc, mixcompdist = "halfuniform")

pi.sim.hu.c <- res.sim.hu.c$fitted_g$pi
mean.sim.hu.c <- (res.sim.hu.c$fitted_g$a + res.sim.hu.c$fitted_g$b) / 2
var.sim.hu.c <- (res.sim.hu.c$fitted_g$b - res.sim.hu.c$fitted_g$a)^2 / 12

pi.sim.hu.nc <- res.sim.hu.nc$fitted_g$pi
mean.sim.hu.nc <- (res.sim.hu.nc$fitted_g$a + res.sim.hu.nc$fitted_g$b) / 2
var.sim.hu.nc <- (res.sim.hu.nc$fitted_g$b - res.sim.hu.nc$fitted_g$a)^2 / 12

## variances of g
var.sim.g.hu.c <- sum(pi.sim.hu.c[-1] * var.sim.hu.c[-1]) +
  sum(pi.sim.hu.c[-1] * mean.sim.hu.c[-1]^2) -
  sum(pi.sim.hu.c[-1] * mean.sim.hu.c[-1])^2
var.sim.g.hu.nc <- sum(pi.sim.hu.nc[-1] * var.sim.hu.nc[-1]) +
  sum(pi.sim.hu.nc[-1] * mean.sim.hu.nc[-1]^2) -
  sum(pi.sim.hu.nc[-1] * mean.sim.hu.nc[-1])^2

round(c(1-pi.sim.hu.c[1],
        1-pi.sim.hu.nc[1],
        sqrt(var.sim.g.hu.c),
        sqrt(var.sim.g.hu.nc)), 3)


##### O-E difference #####

oediff.cnc.c <- oediff.cnc.nc <- s.oediff.cnc.c <- s.oediff.cnc.nc <- rep(NA, nfam)
for(i in 1:nfam){
  oediff.cnc.c[i] <- res.oe.cnc.1fam$OEdiff.C[i]
  oediff.cnc.nc[i] <- res.oe.cnc.1fam$OEdiff.NC[i]
  s.oediff.cnc.c[i] <- sd(res.oe.cnc.1fam.boot[[i]]$OEdiff.C)
  s.oediff.cnc.nc[i] <- sd(res.oe.cnc.1fam.boot[[i]]$OEdiff.NC)
}

##### half-uniform #####

res.sim.diff.hu.c <- ash(oediff.cnc.c, s.oediff.cnc.c, mixcompdist = "halfuniform")
res.sim.diff.hu.nc <- ash(oediff.cnc.nc, s.oediff.cnc.nc, mixcompdist = "halfuniform")

pi.sim.diff.hu.c <- res.sim.diff.hu.c$fitted_g$pi
mean.sim.diff.hu.c <- (res.sim.diff.hu.c$fitted_g$a + res.sim.diff.hu.c$fitted_g$b) / 2
var.sim.diff.hu.c <- (res.sim.diff.hu.c$fitted_g$b - res.sim.diff.hu.c$fitted_g$a)^2 / 12

pi.sim.diff.hu.nc <- res.sim.diff.hu.nc$fitted_g$pi
mean.sim.diff.hu.nc <- (res.sim.diff.hu.nc$fitted_g$a + res.sim.diff.hu.nc$fitted_g$b) / 2
var.sim.diff.hu.nc <- (res.sim.diff.hu.nc$fitted_g$b - res.sim.diff.hu.nc$fitted_g$a)^2 / 12

## variances of g
var.sim.diff.g.hu.c <- sum(pi.sim.diff.hu.c[-1] * var.sim.diff.hu.c[-1]) +
  sum(pi.sim.diff.hu.c[-1] * mean.sim.diff.hu.c[-1]^2) -
  sum(pi.sim.diff.hu.c[-1] * mean.sim.diff.hu.c[-1])^2
var.sim.diff.g.hu.nc <- sum(pi.sim.diff.hu.nc[-1] * var.sim.diff.hu.nc[-1]) +
  sum(pi.sim.diff.hu.nc[-1] * mean.sim.diff.hu.nc[-1]^2) -
  sum(pi.sim.diff.hu.nc[-1] * mean.sim.diff.hu.nc[-1])^2

round(c(1-pi.sim.diff.hu.c[1],
        1-pi.sim.diff.hu.nc[1],
        sqrt(var.sim.diff.g.hu.c),
        sqrt(var.sim.diff.g.hu.nc)), 3)


## Frailty model on noncarriers ##

oe.nc.c <- oe.nc.nc <- s.oe.nc.c <- s.oe.nc.nc <- rep(NA, nfam)
for(i in 1:nfam){
  oe.nc.c[i] <- res.oe.nc.1fam$OE.C[i] - 1
  oe.nc.nc[i] <- res.oe.nc.1fam$OE.NC[i] - 1
  s.oe.nc.c[i] <- sd(res.oe.nc.1fam.boot[[i]]$OE.C)
  s.oe.nc.nc[i] <- sd(res.oe.nc.1fam.boot[[i]]$OE.NC)
}


##### half-uniform #####

res.sim.hu.c <- ash(oe.nc.c, s.oe.nc.c, mixcompdist = "halfuniform")
res.sim.hu.nc <- ash(oe.nc.nc, s.oe.nc.nc, mixcompdist = "halfuniform")

pi.sim.hu.c <- res.sim.hu.c$fitted_g$pi
mean.sim.hu.c <- (res.sim.hu.c$fitted_g$a + res.sim.hu.c$fitted_g$b) / 2
var.sim.hu.c <- (res.sim.hu.c$fitted_g$b - res.sim.hu.c$fitted_g$a)^2 / 12

pi.sim.hu.nc <- res.sim.hu.nc$fitted_g$pi
mean.sim.hu.nc <- (res.sim.hu.nc$fitted_g$a + res.sim.hu.nc$fitted_g$b) / 2
var.sim.hu.nc <- (res.sim.hu.nc$fitted_g$b - res.sim.hu.nc$fitted_g$a)^2 / 12

## variances of g
var.sim.g.hu.c <- sum(pi.sim.hu.c[-1] * var.sim.hu.c[-1]) +
  sum(pi.sim.hu.c[-1] * mean.sim.hu.c[-1]^2) -
  sum(pi.sim.hu.c[-1] * mean.sim.hu.c[-1])^2
var.sim.g.hu.nc <- sum(pi.sim.hu.nc[-1] * var.sim.hu.nc[-1]) +
  sum(pi.sim.hu.nc[-1] * mean.sim.hu.nc[-1]^2) -
  sum(pi.sim.hu.nc[-1] * mean.sim.hu.nc[-1])^2

round(c(1-pi.sim.hu.c[1],
        1-pi.sim.hu.nc[1],
        sqrt(var.sim.g.hu.c),
        sqrt(var.sim.g.hu.nc)), 3)


##### O-E difference #####

oediff.nc.c <- oediff.nc.nc <- s.oediff.nc.c <- s.oediff.nc.nc <- rep(NA, nfam)
for(i in 1:nfam){
  oediff.nc.c[i] <- res.oe.nc.1fam$OEdiff.C[i]
  oediff.nc.nc[i] <- res.oe.nc.1fam$OEdiff.NC[i]
  s.oediff.nc.c[i] <- sd(res.oe.nc.1fam.boot[[i]]$OEdiff.C)
  s.oediff.nc.nc[i] <- sd(res.oe.nc.1fam.boot[[i]]$OEdiff.NC)
}

##### half-uniform #####

res.sim.diff.hu.c <- ash(oediff.nc.c, s.oediff.nc.c, mixcompdist = "halfuniform")
res.sim.diff.hu.nc <- ash(oediff.nc.nc, s.oediff.nc.nc, mixcompdist = "halfuniform")

pi.sim.diff.hu.c <- res.sim.diff.hu.c$fitted_g$pi
mean.sim.diff.hu.c <- (res.sim.diff.hu.c$fitted_g$a + res.sim.diff.hu.c$fitted_g$b) / 2
var.sim.diff.hu.c <- (res.sim.diff.hu.c$fitted_g$b - res.sim.diff.hu.c$fitted_g$a)^2 / 12

pi.sim.diff.hu.nc <- res.sim.diff.hu.nc$fitted_g$pi
mean.sim.diff.hu.nc <- (res.sim.diff.hu.nc$fitted_g$a + res.sim.diff.hu.nc$fitted_g$b) / 2
var.sim.diff.hu.nc <- (res.sim.diff.hu.nc$fitted_g$b - res.sim.diff.hu.nc$fitted_g$a)^2 / 12

## variances of g
var.sim.diff.g.hu.c <- sum(pi.sim.diff.hu.c[-1] * var.sim.diff.hu.c[-1]) +
  sum(pi.sim.diff.hu.c[-1] * mean.sim.diff.hu.c[-1]^2) -
  sum(pi.sim.diff.hu.c[-1] * mean.sim.diff.hu.c[-1])^2
var.sim.diff.g.hu.nc <- sum(pi.sim.diff.hu.nc[-1] * var.sim.diff.hu.nc[-1]) +
  sum(pi.sim.diff.hu.nc[-1] * mean.sim.diff.hu.nc[-1]^2) -
  sum(pi.sim.diff.hu.nc[-1] * mean.sim.diff.hu.nc[-1])^2

round(c(1-pi.sim.diff.hu.c[1],
        1-pi.sim.diff.hu.nc[1],
        sqrt(var.sim.diff.g.hu.c),
        sqrt(var.sim.diff.g.hu.nc)), 3)


#### No frailty ####

oe.none.c <- oe.none.nc <- s.oe.none.c <- s.oe.none.nc <- rep(NA, nfam)
for(i in 1:nfam){
  oe.none.c[i] <- res.oe.none.1fam$OE.C[i] - 1
  oe.none.nc[i] <- res.oe.none.1fam$OE.NC[i] - 1
  s.oe.none.c[i] <- sd(res.oe.none.1fam.boot[[i]]$OE.C)
  s.oe.none.nc[i] <- sd(res.oe.none.1fam.boot[[i]]$OE.NC)
}


##### half-uniform #####

res.sim.hu.c <- ash(oe.none.c, s.oe.none.c, mixcompdist = "halfuniform")
res.sim.hu.nc <- ash(oe.none.nc, s.oe.none.nc, mixcompdist = "halfuniform")

pi.sim.hu.c <- res.sim.hu.c$fitted_g$pi
mean.sim.hu.c <- (res.sim.hu.c$fitted_g$a + res.sim.hu.c$fitted_g$b) / 2
var.sim.hu.c <- (res.sim.hu.c$fitted_g$b - res.sim.hu.c$fitted_g$a)^2 / 12

pi.sim.hu.nc <- res.sim.hu.nc$fitted_g$pi
mean.sim.hu.nc <- (res.sim.hu.nc$fitted_g$a + res.sim.hu.nc$fitted_g$b) / 2
var.sim.hu.nc <- (res.sim.hu.nc$fitted_g$b - res.sim.hu.nc$fitted_g$a)^2 / 12

## variances of g
var.sim.g.hu.c <- sum(pi.sim.hu.c[-1] * var.sim.hu.c[-1]) +
  sum(pi.sim.hu.c[-1] * mean.sim.hu.c[-1]^2) -
  sum(pi.sim.hu.c[-1] * mean.sim.hu.c[-1])^2
var.sim.g.hu.nc <- sum(pi.sim.hu.nc[-1] * var.sim.hu.nc[-1]) +
  sum(pi.sim.hu.nc[-1] * mean.sim.hu.nc[-1]^2) -
  sum(pi.sim.hu.nc[-1] * mean.sim.hu.nc[-1])^2

round(c(1-pi.sim.hu.c[1],
        1-pi.sim.hu.nc[1],
        sqrt(var.sim.g.hu.c),
        sqrt(var.sim.g.hu.nc)), 3)


##### O-E difference #####

oediff.none.c <- oediff.none.nc <- s.oediff.none.c <- s.oediff.none.nc <- rep(NA, nfam)
for(i in 1:nfam){
  oediff.none.c[i] <- res.oe.none.1fam$OEdiff.C[i]
  oediff.none.nc[i] <- res.oe.none.1fam$OEdiff.NC[i]
  s.oediff.none.c[i] <- sd(res.oe.none.1fam.boot[[i]]$OEdiff.C)
  s.oediff.none.nc[i] <- sd(res.oe.none.1fam.boot[[i]]$OEdiff.NC)
}

##### half-uniform #####

res.sim.diff.hu.c <- ash(oediff.none.c, s.oediff.none.c, mixcompdist = "halfuniform")
res.sim.diff.hu.nc <- ash(oediff.none.nc, s.oediff.none.nc, mixcompdist = "halfuniform")

pi.sim.diff.hu.c <- res.sim.diff.hu.c$fitted_g$pi
mean.sim.diff.hu.c <- (res.sim.diff.hu.c$fitted_g$a + res.sim.diff.hu.c$fitted_g$b) / 2
var.sim.diff.hu.c <- (res.sim.diff.hu.c$fitted_g$b - res.sim.diff.hu.c$fitted_g$a)^2 / 12

pi.sim.diff.hu.nc <- res.sim.diff.hu.nc$fitted_g$pi
mean.sim.diff.hu.nc <- (res.sim.diff.hu.nc$fitted_g$a + res.sim.diff.hu.nc$fitted_g$b) / 2
var.sim.diff.hu.nc <- (res.sim.diff.hu.nc$fitted_g$b - res.sim.diff.hu.nc$fitted_g$a)^2 / 12

## variances of g
var.sim.diff.g.hu.c <- sum(pi.sim.diff.hu.c[-1] * var.sim.diff.hu.c[-1]) +
  sum(pi.sim.diff.hu.c[-1] * mean.sim.diff.hu.c[-1]^2) -
  sum(pi.sim.diff.hu.c[-1] * mean.sim.diff.hu.c[-1])^2
var.sim.diff.g.hu.nc <- sum(pi.sim.diff.hu.nc[-1] * var.sim.diff.hu.nc[-1]) +
  sum(pi.sim.diff.hu.nc[-1] * mean.sim.diff.hu.nc[-1]^2) -
  sum(pi.sim.diff.hu.nc[-1] * mean.sim.diff.hu.nc[-1])^2

round(c(1-pi.sim.diff.hu.c[1],
        1-pi.sim.diff.hu.nc[1],
        sqrt(var.sim.diff.g.hu.c),
        sqrt(var.sim.diff.g.hu.nc)), 3)



save(res.oe.cnc.1fam, res.oe.cnc.1fam.boot,
     res.oe.c.1fam, res.oe.c.1fam.boot,
     res.oe.nc.1fam, res.oe.nc.1fam.boot,
     res.oe.none.1fam, res.oe.none.1fam.boot,
     file = "OE_Sim_1Fam.RData")


