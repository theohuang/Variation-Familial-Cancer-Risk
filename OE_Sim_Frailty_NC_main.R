## Simulation for O/E ratios for a frailty model
## Frailties on non-carriers only
## Last updated: July 1, 2020

rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(dplyr)
library(abind)
library(ggplot2)
library(tidyverse)

## loading BayesMendel colorectal and endometrial cancer penetrances and
## hazard for death from other causes
load("penet.mmr.net.RData")
load("death.othercauses.RData")

source("OE Functions.R")
source("Estimating Functions Discrete.R")

source(paste0(getwd(), "/Generating Families Functions/sim.simFam.R"))
source(paste0(getwd(), "/Generating Families Functions/genCancerPen.R"))
source(paste0(getwd(), "/Generating Families Functions/genOtherDeathPen.R"))
source(paste0(getwd(), "/Generating Families Functions/sim.buildGenoMat.R"))
source(paste0(getwd(), "/Generating Families Functions/sim.linkParents.R"))
source(paste0(getwd(), "/Generating Families Functions/sim.simCurAgeVar.R"))
source(paste0(getwd(), "/Generating Families Functions/sim.simCancerVars.R"))
source(paste0(getwd(), "/Generating Families Functions/sim.buildBranchOfAlleleMats.R"))
source(paste0(getwd(), "/Generating Families Functions/helpers.R"))

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


nfam <- 53; nfam.w <- 4
fams <- vector("list", nfam)
for(i in 1:nfam){fams[[i]] <- vector("list", nfam.w)}
cnames <- c("FamID", "W", "O.C", "E.C", "OE.C", "O.NC", "E.NC",
            "OE.NC", "OEdiff.C", "OEdiff.NC", "FamSize")
res.oe.sim <- setNames(data.frame(matrix(0, nfam, length(cnames))), cnames)
res.oe.sim$FamID <- 1:nfam

start <- Sys.time()
set.seed(999)
for(k in 1:nfam){
  print(k)
  
  res.oe.sim$W[k] <- sample(seq(-2, 2, 0.1), 1)
  ## frailty only for colorectal cancer (and only carriers)
  hzdF <- hzd0F; hzdF$ColorC[, 1] <- hzd02hzd(hzd0F$ColorC, w = res.oe.sim$W[k])[, 1]
  hzdM <- hzd0M; hzdM$ColorC[, 1] <- hzd02hzd(hzd0M$ColorC, w = res.oe.sim$W[k])[, 1]
  penF <- lapply(hzdF, hzd2pen)
  penM <- lapply(hzdM, hzd2pen)
  CP <- genCancerPen(mutations, cancers, penF, penM, maxK = length(mutations), age.last = 95)
  
  ## getting "base" families (10 per frailty) to obtain genotypes and current ages
  # keep generating until we have a family with at least one carrier
  fam0 <- genoMat <- genderPro <- geno <- fam <- res.exp <-res.obs <- 
    LIK <- probs <- vector("list", nfam.w)
  for(j in 1:nfam.w){
    nSibsPatern <- sample(1:4, 2, replace = TRUE)
    nSibsMatern <- sample(1:4, 2, replace = TRUE)
    nSibs <- sample(1:4, 2, replace = TRUE)
    nGrandchild <- matrix(sample(1:4, 2 * (sum(nSibs) + 1), replace = TRUE),
                          sum(nSibs) + 1, 2)
    repeat{
      fam0[[j]] <- sim.simFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild,
                              af, CP0, includeGeno = TRUE, age.max = 94, age.min = 2)
      genoMat[[j]] <- select(fam0[[j]], MLH1, MSH2, MSH6)
      if(sum(genoMat[[j]]) > 0) break
    }
    
    genderPro[[j]] <- ifelse(filter(fam0[[j]], isProband == 1)$Gender == 1, "Male", "Female")
    geno[[j]] <- rep(0, nrow(genoMat[[j]]))
    for(i in 1:nrow(genoMat[[j]])){
      geno[[j]][i] <- ifelse(all(genoMat[[j]][i, ] == c(0, 0, 0)), 1,
                             ifelse(all(genoMat[[j]][i, ] == c(1, 0, 0)), 2,
                                    ifelse(all(genoMat[[j]][i, ] == c(0, 1, 0)), 3,
                                           ifelse(all(genoMat[[j]][i, ] == c(0, 0, 1)), 4,
                                                  ifelse(all(genoMat[[j]][i, ] == c(1, 1, 0)), 5,
                                                         ifelse(all(genoMat[[j]][i, ] == c(1, 0, 1)), 6,
                                                                ifelse(all(genoMat[[j]][i, ] == c(0, 1, 1)), 7, 8)))))))
    }
    
    fam[[j]] <- sim.simFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild,
                           af, CP, includeGeno = TRUE, age.max = 94, age.min = 2,
                           genderPro = genderPro[[j]], genoMat = genoMat[[j]],
                           CurAge = fam0[[j]]$CurAge, affTime = TRUE)
    fam[[j]] <- mutate(fam[[j]], AffectedColon = isAffColorC, AffectedEndometrium = isAffEndomC,
                       AgeColon = AgeColorC, AgeEndometrium = AgeEndomC, FamID = j)
    fams[[k]][[j]] <- fam[[j]]
    
    ## getting expected and observed for each individual
    # with censoring
    res.exp[[j]] <- expect(fam[[j]], CP0$cancerMDens[1:94, , "ColorC"], CP0$cancerFDens[1:94, , "ColorC"])
    res.obs[[j]] <- observe(fam[[j]], CP0$cancerMDens[1:94, , "ColorC"], CP0$cancerFDens[1:94, , "ColorC"])
    
    ### O/E ratios
    ## using carrier probabilities
    # with censoring
    LIK[[j]] <- estLik(fam[[j]], CP0, ODP)
    probs[[j]] <- setNames(data.frame(matrix(0, nrow(fam[[j]]), ncol(CP0$cancerFDens[, , 1]))),
                           colnames(CP0$cancerFDens[, , 1]))
    for(i in 1:nrow(fam[[j]])){
      probs[[j]][i, ] <- pp.peelingParing(fam[[j]], af, LIK[[j]], length(mutations),
                                          counselee.id = fam[[j]]$ID[i])
    }
  }
  
  ## combining the nfam.w families together
  res.oe.sim$FamSize[k] <- mean(unlist(lapply(fam, nrow)))
  res.oe.sim$O.C[k] <- mean(unlist(lapply(Map("*", res.obs, probs), function(x) sum(x[, -1]))))
  res.oe.sim$E.C[k] <- mean(unlist(lapply(Map("*", res.exp, probs), function(x) sum(x[, -1]))))
  res.oe.sim$O.NC[k] <- mean(unlist(lapply(Map("*", res.obs, probs), function(x) sum(x[, 1]))))
  res.oe.sim$E.NC[k] <- mean(unlist(lapply(Map("*", res.exp, probs), function(x) sum(x[, 1]))))
  res.oe.sim$OE.C[k] <- res.oe.sim$O.C[k] / res.oe.sim$E.C[k]
  res.oe.sim$OE.NC[k] <- res.oe.sim$O.NC[k] / res.oe.sim$E.NC[k]
  res.oe.sim$OEdiff.C[k] <- (res.oe.sim$O.C[k] - res.oe.sim$E.C[k]) / res.oe.sim$FamSize[k]
  res.oe.sim$OEdiff.NC[k] <- (res.oe.sim$O.NC[k] - res.oe.sim$E.NC[k]) / res.oe.sim$FamSize[k]
}
difftime(Sys.time(), start, units = "secs")


#### Bootstraping ####

nboot <- 100
res.oe.sim.boot <- vector("list", nfam)

start <- Sys.time()
set.seed(1)
boot.w <- vector("list", nfam)
for(i in 1:nfam) boot.w[[i]] <- vector("list", nfam.w)
for(i in 1:nfam){
  print(i)
  for(j in 1:nfam.w){
    boot.w[[i]][[j]] <- mutate(oe.boot.sim(fams[[i]][[j]], CP0, ODP, af, mutations, nboot, seed = 999),
                               FamSize = nrow(fams[[i]][[j]]))
  }
  boot.sum <- boot.w[[i]][[1]]
  for(l in 1:nboot){
    boot.sum[l, ] <- (boot.w[[i]][[1]][l, ] + boot.w[[i]][[2]][l, ] + boot.w[[i]][[3]][l, ] +
                        boot.w[[i]][[4]][l, ]) / 4
  }
  boot.sum <- mutate(boot.sum, FamID = i, OE.C = Obs.C / Exp.C, OE.NC = Obs.NC / Exp.NC,
                     OEdiff.C = (Obs.C - Exp.C) / FamSize, OEdiff.NC = (Obs.NC - Exp.NC) / FamSize)
  res.oe.sim.boot[[i]] <- boot.sum
}
difftime(Sys.time(), start, units = "secs")


res.oe.sim <- mutate(res.oe.sim, OE.C.lo = NA, OE.C.hi = NA,
                     OE.NC.lo = NA, OE.NC.hi = NA,
                     OEdiff.C.lo = NA, OEdiff.C.hi = NA,
                     OEdiff.NC.lo = NA, OEdiff.NC.hi = NA)
for(i in 1:nfam){
  res.oe.sim$OE.C.lo[i] <- quantile(res.oe.sim.boot[[i]]$OE.C, 0.025)
  res.oe.sim$OE.C.hi[i] <- quantile(res.oe.sim.boot[[i]]$OE.C, 0.975)
  res.oe.sim$OE.NC.lo[i] <- quantile(res.oe.sim.boot[[i]]$OE.NC, 0.025)
  res.oe.sim$OE.NC.hi[i] <- quantile(res.oe.sim.boot[[i]]$OE.NC, 0.975)
  res.oe.sim$OEdiff.C.lo[i] <- quantile(res.oe.sim.boot[[i]]$OEdiff.C, 0.025)
  res.oe.sim$OEdiff.C.hi[i] <- quantile(res.oe.sim.boot[[i]]$OEdiff.C, 0.975)
  res.oe.sim$OEdiff.NC.lo[i] <- quantile(res.oe.sim.boot[[i]]$OEdiff.NC, 0.025)
  res.oe.sim$OEdiff.NC.hi[i] <- quantile(res.oe.sim.boot[[i]]$OEdiff.NC, 0.975)
}

# res.oe.sim$FamSize <- unlist(lapply(fams, nrow))

# res.oe.sim <- mutate(res.oe.sim,
#                      OEdiff.C.cp = (O.C.cp - E.C.cp) / FamSize,
#                      OEdiff.NC.cp = (O.NC.cp - E.NC.cp) / FamSize)
# ## O-E difference
# for(i in 1:nfam){
#   res.oe.sim$OEdiff.C.cp.lo[i] <- quantile((res.oe.sim.boot[[i]]$Obs.C - res.oe.sim.boot[[i]]$Exp.C) / res.oe.sim$FamSize[i], 0.025)
#   res.oe.sim$OEdiff.C.cp.hi[i] <- quantile((res.oe.sim.boot[[i]]$Obs.C - res.oe.sim.boot[[i]]$Exp.C) / res.oe.sim$FamSize[i], 0.975)
#   res.oe.sim$OEdiff.NC.cp.lo[i] <- quantile((res.oe.sim.boot[[i]]$Obs.NC - res.oe.sim.boot[[i]]$Exp.NC) / res.oe.sim$FamSize[i], 0.025)
#   res.oe.sim$OEdiff.NC.cp.hi[i] <- quantile((res.oe.sim.boot[[i]]$Obs.NC - res.oe.sim.boot[[i]]$Exp.NC) / res.oe.sim$FamSize[i], 0.975)
# }


#### Plots #####

## plot of O/E ratios with confidence intervals

ggplot(res.oe.sim, aes(OE.C, OE.NC)) +
  geom_point(aes(color = W), size = 3) +
  geom_errorbar(aes(ymin = OE.NC.lo, ymax = OE.NC.hi), width = 0.005, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OE.C.lo, xmax = OE.C.hi), height = 0.005, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "O/E, carriers",
       y = "O/E, non-carriers")

# ggplot(mutate(res.oe.sim, diff.c = O.C - E.C,
#               diff.nc = O.NC - E.NC), aes(diff.c, diff.nc)) +
#   geom_point(aes(color = W), size = 3) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
#   labs(x = "O-E, carriers",
#        y = "O-E, non-carriers")


## O-E difference
ggplot(res.oe.sim, aes(OEdiff.C, OEdiff.NC)) +
  geom_point(aes(color = W), size = 3) +
  geom_errorbar(aes(ymin = OEdiff.NC.lo, ymax = OEdiff.NC.hi),
                width = 0.001, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OEdiff.C.lo, xmax = OEdiff.C.hi),
                 height = 0.0003, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "(O-E)/n, carriers",
       y = "(O-E)/n, non-carriers")


## plot of confidence interval
ggplot(gather(mutate(res.oe.sim, Width.C = OE.C.hi - OE.C.lo,
                     Width.NC = OE.NC.hi - OE.NC.lo),
              Type, Width, c("Width.C", "Width.NC")),
       aes(FamSize, Width)) +
  geom_point(aes(color = Type)) +
  labs(x = "Family Size", y = "Confidence Interval Width") +
  scale_color_discrete(name = "", labels = c("Carriers", "Non-carriers"))


#### ashr ####

library(ashr)

oe.sim.c <- oe.sim.nc <- s.oe.sim.c <- s.oe.sim.nc <- rep(NA, nfam)
for(i in 1:nfam){
  oe.sim.c[i] <- res.oe.sim$OE.C[i] - 1
  oe.sim.nc[i] <- res.oe.sim$OE.NC[i] - 1
  s.oe.sim.c[i] <- sd(res.oe.sim.boot[[i]]$OE.C)
  s.oe.sim.nc[i] <- sd(res.oe.sim.boot[[i]]$OE.NC)
}

# ##### normal #####
# 
# res.sim.norm.c <- ash(oe.sim.c, s.oe.sim.c, mixcompdist = "normal")
# res.sim.norm.nc <- ash(oe.sim.nc, s.oe.sim.nc, mixcompdist = "normal")
# 
# pi.sim.norm.c <- res.sim.norm.c$fitted_g$pi
# var.sim.norm.c <- res.sim.norm.c$fitted_g$sd^2
# 
# pi.sim.norm.nc <- res.sim.norm.nc$fitted_g$pi
# var.sim.norm.nc <- res.sim.norm.nc$fitted_g$sd^2
# 
# ## var.simiances of g
# var.sim.g.norm.c <- sum(pi.sim.norm.c[-1] * var.sim.norm.c[-1])
# var.sim.g.norm.nc <- sum(pi.sim.norm.nc[-1] * var.sim.norm.nc[-1])
# 
# sqrt(var.sim.g.norm.c)
# sqrt(var.sim.g.norm.nc)
# 
# 


##### half-uniform #####

res.sim.hu.c <- ash(oe.sim.c, s.oe.sim.c, mixcompdist = "halfuniform")
res.sim.hu.nc <- ash(oe.sim.nc, s.oe.sim.nc, mixcompdist = "halfuniform")

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

round(c(1 - pi.sim.hu.c[1],
        1 - pi.sim.hu.nc[1],
        sqrt(var.sim.g.hu.c),
        sqrt(var.sim.g.hu.nc)), 3)


##### O-E difference #####

oediff.sim.c <- oediff.sim.nc <- s.oediff.sim.c <- s.oediff.sim.nc <- rep(NA, nfam)
for(i in 1:nfam){
  oediff.sim.c[i] <- res.oe.sim$OEdiff.C[i]
  oediff.sim.nc[i] <- res.oe.sim$OEdiff.NC[i]
  s.oediff.sim.c[i] <- sd(res.oe.sim.boot[[i]]$OEdiff.C)
  s.oediff.sim.nc[i] <- sd(res.oe.sim.boot[[i]]$OEdiff.NC)
}


# #### normal ####
# 
# res.sim.diff.norm.c <- ash(oediff.sim.c, s.oediff.sim.c, mixcompdist = "normal")
# res.sim.diff.norm.nc <- ash(oediff.sim.nc, s.oediff.sim.nc, mixcompdist = "normal")
# 
# pi.sim.diff.norm.c <- res.sim.diff.norm.c$fitted_g$pi
# var.sim.diff.norm.c <- res.sim.diff.norm.c$fitted_g$sd^2
# 
# pi.sim.diff.norm.nc <- res.sim.diff.norm.nc$fitted_g$pi
# var.sim.diff.norm.nc <- res.sim.diff.norm.nc$fitted_g$sd^2
# 
# ## var.simiances of g
# var.sim.diff.g.norm.c <- sum(pi.sim.diff.norm.c[-1] * var.sim.diff.norm.c[-1])
# var.sim.diff.g.norm.nc <- sum(pi.sim.diff.norm.nc[-1] * var.sim.diff.norm.nc[-1])
# 
# pi.sim.diff.norm.c[1]
# pi.sim.diff.norm.nc[1]
# sqrt(var.sim.diff.g.norm.c)
# sqrt(var.sim.diff.g.norm.nc)


##### half-uniform #####

res.sim.diff.hu.c <- ash(oediff.sim.c, s.oediff.sim.c, mixcompdist = "halfuniform")
res.sim.diff.hu.nc <- ash(oediff.sim.nc, s.oediff.sim.nc, mixcompdist = "halfuniform")

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

round(c(1 - pi.sim.diff.hu.c[1],
        1 - pi.sim.diff.hu.nc[1],
        sqrt(var.sim.diff.g.hu.c),
        sqrt(var.sim.diff.g.hu.nc)), 3)


save(res.oe.sim, res.oe.sim.boot, fams, boot.w, res.sim.hu.c, res.sim.hu.nc,
     res.sim.diff.hu.c, res.sim.diff.hu.nc, file = "OE_Sim_NC_Main.RData")
