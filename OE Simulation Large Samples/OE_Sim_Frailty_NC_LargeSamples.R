## Simulation for O/E ratios for a frailty model
## Large samples
## Frailties on only non-carriers
## Last updated: July 16, 2020

rm(list = ls())
a1 <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(999 + a1)

library(dplyr)
library(abind)

## loading BayesMendel colorectal and endometrial cancer penetrances and
## hazard for death from other causes
load("penet.mmr.net.RData")
load("death.othercauses.RData")

source("OE Functions.R")
source("Estimating Functions Discrete.R")

source(paste(getwd(), "/Generating Families Functions/sim.simFam.R", sep = ""))
source(paste(getwd(), "/Generating Families Functions/genCancerPen.R", sep = ""))
source(paste(getwd(), "/Generating Families Functions/genOtherDeathPen.R", sep = ""))
source(paste(getwd(), "/Generating Families Functions/sim.buildGenoMat.R", sep = ""))
source(paste(getwd(), "/Generating Families Functions/sim.linkParents.R", sep = ""))
source(paste(getwd(), "/Generating Families Functions/sim.simCurAgeVar.R", sep = ""))
source(paste(getwd(), "/Generating Families Functions/sim.simCancerVars.R", sep = ""))
source(paste(getwd(), "/Generating Families Functions/sim.buildBranchOfAlleleMats.R", sep = ""))
source(paste(getwd(), "/Generating Families Functions/helpers.R", sep = ""))

mutations <- c("MLH1", "MSH2", "MSH6")
cancers <- c("ColorC", "EndomC")
genos <- c("M000", "M100", "M010", "M001")

pen2hzd <- function(pen){
  apply(pen, 2, function(x) x / c(1, 1 - cumsum(x)[-length(x)]))
  # hzd <- pen
  # for(i in 1:ncol(pen)){
  #   hzd[, i] <- pen[, i] / c(1, 1 - cumsum(pen[, i])[-1])
  # }
  # return(hzd)
}

hzd2pen <- function(hzd){
  surv <- apply(1 - hzd, 2, cumprod)
  pen <- rbind(hzd[1, ], surv[-nrow(hzd), ] * hzd[-1, ])
  # for(i in 1:ncol(hzd)){
  #   surv[, i] <- cumprod(1 - hzd[, i])
  #   hzd0[, i] <- pen[, i] / c(1, 1 - cumsum(pen[, i])[-1])
  # }
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

# af <- setNames(unlist(lapply(MMRparams()$allef, "[[", 2)),
#                mutations)
## using allele frequencies of 0.1 to increase the number of carriers
af <- setNames(rep(0.1, 3), mutations)
## baseline penetrance (MMRpro)
CP0 <- genCancerPen(mutations, cancers, pen0F, pen0M, maxK = length(mutations), age.last = 95)

w.list <- seq(-2, 2, 0.1)

nSibsPatern <- c(2, 2)
nSibsMatern <- c(2, 2)
nSibs <- c(2, 2)
nGrandchild <- matrix(rep(2, 10), 5, 2)

ODP <- genOtherDeathPen(cancers,
                        mutate(select(death.othercauses[1:94, ],
                                      femaleCRC, maleCRC, uterus),
                               femaleColorC = femaleCRC, maleColorC = maleCRC,
                               femaleEndomC = uterus, maleEndomC = uterus),
                        age.max = 94)


nboot <- 20
res.oe <- fams <- vector("list", nboot)

start <- Sys.time()
for(k in 1:nboot){
  print(k)
  ## getting "base" family to obtain genotypes and current ages
  # keep generating until we have a family with at least one carrier
  repeat{
    fam0 <- sim.simFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild,
                       af, CP0, includeGeno = TRUE, age.max = 94, age.min = 2)
    genoMat <- select(fam0, MLH1, MSH2, MSH6)
    if(sum(genoMat) > 0) break
  }
  genderPro <- ifelse(filter(fam0, isProband == 1)$Gender == 1, "Male", "Female")
  geno <- rep(0, nrow(genoMat))
  for(i in 1:nrow(genoMat)){
    geno[i] <- ifelse(all(genoMat[i, ] == c(0, 0, 0)), 1,
                      ifelse(all(genoMat[i, ] == c(1, 0, 0)), 2,
                             ifelse(all(genoMat[i, ] == c(0, 1, 0)), 3,
                                    ifelse(all(genoMat[i, ] == c(0, 0, 1)), 4,
                                           ifelse(all(genoMat[i, ] == c(1, 1, 0)), 5,
                                                  ifelse(all(genoMat[i, ] == c(1, 0, 1)), 6,
                                                         ifelse(all(genoMat[i, ] == c(0, 1, 1)), 7, 8)))))))
  }
  cnames <- c("W", "O.C", "E.C", "OE.C", "O.NC", "E.NC", "OE.NC",
              paste(c("O.C", "E.C", "OE.C", "O.NC", "E.NC", "OE.NC"), ".cp", sep = ""),
              paste(c("O.C", "E.C", "OE.C", "O.NC", "E.NC", "OE.NC"), ".nocens", sep = ""),
              paste(c("O.C", "E.C", "OE.C", "O.NC", "E.NC", "OE.NC"), ".cp.nocens", sep = ""))
  res.oe[[k]] <- setNames(data.frame(matrix(0, length(w.list), length(cnames))),
                          cnames)
  res.oe[[k]]$W <- w.list
  for(i in 1:length(w.list)){
    ## frailty only for colorectal cancer (and only non-carriers)
    hzdF <- hzd0F; hzdF$ColorC[, 1] <- hzd02hzd(hzd0F$ColorC, w = w.list[i])[, 1]
    hzdM <- hzd0M; hzdM$ColorC[, 1] <- hzd02hzd(hzd0M$ColorC, w = w.list[i])[, 1]
    penF <- lapply(hzdF, hzd2pen)
    penM <- lapply(hzdM, hzd2pen)
    CP <- genCancerPen(mutations, cancers, penF, penM, maxK = length(mutations), age.last = 95)
    fam <- sim.simFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild,
                      af, CP, includeGeno = TRUE, age.max = 94, age.min = 2,
                      genderPro = genderPro, genoMat = genoMat,
                      CurAge = fam0$CurAge, affTime = TRUE)
    fam <- mutate(fam, AffectedColon = isAffColorC, AffectedEndometrium = isAffEndomC,
                  AgeColon = AgeColorC, AgeEndometrium = AgeEndomC)
    fams[[k]] <- fam
    
    # family if there was no censoring (we observe everyone until age 94)
    fam.nocens <- mutate(fam, AgeColon = replace(AffAgeColorC, AffAgeColorC == 95, 94),
                         AgeEndometrium = replace(AffAgeEndomC, AffAgeEndomC == 95, 94),
                         AffectedColon = ifelse(AffAgeColorC < 95, 1, 0),
                         AffectedEndometrium = ifelse(AffAgeEndomC < 95, 1, 0))
    
    ## getting expected and observed for each individual
    # with censoring
    res.exp <- expect(fam, CP0$cancerMDens[1:94, , "ColorC"], CP0$cancerFDens[1:94, , "ColorC"])
    res.obs <- observe(fam, CP0$cancerMDens[1:94, , "ColorC"], CP0$cancerFDens[1:94, , "ColorC"])
    # without censoring
    res.exp.nocens <- expect(fam.nocens, CP0$cancerMDens[1:94, , "ColorC"], CP0$cancerFDens[1:94, , "ColorC"])
    res.obs.nocens <- observe(fam.nocens, CP0$cancerMDens[1:94, , "ColorC"], CP0$cancerFDens[1:94, , "ColorC"])
    
    ### O/E ratios
    ## using the true genotypes
    fam.o <- fam.e <- fam.o.nocens <- fam.e.nocens <- rep(0, nrow(fam))
    for(j in 1:nrow(fam)){
      fam.o[j] <- res.obs[j, geno[j]]
      fam.e[j] <- res.exp[j, geno[j]]
      fam.o.nocens[j] <- res.obs.nocens[j, geno[j]]
      fam.e.nocens[j] <- res.exp.nocens[j, geno[j]]
    }
    
    # with censoring
    res.oe[[k]]$O.C[i] <- sum(fam.o[geno != 1])
    res.oe[[k]]$E.C[i] <- sum(fam.e[geno != 1])
    res.oe[[k]]$O.NC[i] <- sum(fam.o[geno == 1])
    res.oe[[k]]$E.NC[i] <- sum(fam.e[geno == 1])
    res.oe[[k]]$OE.C[i] <- res.oe[[k]]$O.C[i] / res.oe[[k]]$E.C[i]
    res.oe[[k]]$OE.NC[i] <- res.oe[[k]]$O.NC[i] / res.oe[[k]]$E.NC[i]
    
    # no censoring
    res.oe[[k]]$O.C.nocens[i] <- sum(fam.o.nocens[geno != 1])
    res.oe[[k]]$E.C.nocens[i] <- sum(fam.e.nocens[geno != 1])
    res.oe[[k]]$O.NC.nocens[i] <- sum(fam.o.nocens[geno == 1])
    res.oe[[k]]$E.NC.nocens[i] <- sum(fam.e.nocens[geno == 1])
    res.oe[[k]]$OE.C.nocens[i] <- res.oe[[k]]$O.C.nocens[i] / res.oe[[k]]$E.C.nocens[i]
    res.oe[[k]]$OE.NC.nocens[i] <- res.oe[[k]]$O.NC.nocens[i] / res.oe[[k]]$E.NC.nocens[i]
    
    ## using carrier probabilities
    # with censoring
    LIK <- estLik(fam, CP0, ODP)
    probs <- setNames(data.frame(matrix(0, nrow(fam), ncol(CP0$cancerFDens[, , 1]))),
                      colnames(CP0$cancerFDens[, , 1]))
    for(j in 1:nrow(fam)){
      probs[j, ] <- pp.peelingParing(fam, af, LIK, length(mutations),
                                     counselee.id = fam$ID[j])
    }
    
    res.oe[[k]]$O.C.cp[i] <- sum(res.obs[, -1] * probs[, -1])
    res.oe[[k]]$E.C.cp[i] <- sum(res.exp[, -1] * probs[, -1])
    res.oe[[k]]$O.NC.cp[i] <- sum(res.obs[, 1] * probs[, 1])
    res.oe[[k]]$E.NC.cp[i] <- sum(res.exp[, 1] * probs[, 1])
    res.oe[[k]]$OE.C.cp[i] <- res.oe[[k]]$O.C.cp[i] / res.oe[[k]]$E.C.cp[i]
    res.oe[[k]]$OE.NC.cp[i] <- res.oe[[k]]$O.NC.cp[i] / res.oe[[k]]$E.NC.cp[i]
    
    # without censoring
    LIK.nocens <- estLik(fam.nocens, CP0, ODP)
    probs.nocens <- setNames(data.frame(matrix(0, nrow(fam.nocens), ncol(CP0$cancerFDens[, , 1]))),
                             colnames(CP0$cancerFDens[, , 1]))
    for(j in 1:nrow(fam.nocens)){
      probs.nocens[j, ] <- pp.peelingParing(fam.nocens, af, LIK.nocens, length(mutations),
                                            counselee.id = fam.nocens$ID[j])
    }
    
    res.oe[[k]]$O.C.cp.nocens[i] <- sum(res.obs.nocens[, -1] * probs.nocens[, -1])
    res.oe[[k]]$E.C.cp.nocens[i] <- sum(res.exp.nocens[, -1] * probs.nocens[, -1])
    res.oe[[k]]$O.NC.cp.nocens[i] <- sum(res.obs.nocens[, 1] * probs.nocens[, 1])
    res.oe[[k]]$E.NC.cp.nocens[i] <- sum(res.exp.nocens[, 1] * probs.nocens[, 1])
    res.oe[[k]]$OE.C.cp.nocens[i] <- res.oe[[k]]$O.C.cp.nocens[i] / res.oe[[k]]$E.C.cp.nocens[i]
    res.oe[[k]]$OE.NC.cp.nocens[i] <- res.oe[[k]]$O.NC.cp.nocens[i] / res.oe[[k]]$E.NC.cp.nocens[i]
  }
}
difftime(Sys.time(), start, units = "secs")


save(res.oe, fams, file = paste(getwd(), "/Extended Frailty/OE Simulation/NonCarriers/oe_sim_nc_", a1, ".RData", sep = ""))

