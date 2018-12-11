## Looking at observed/expected ratios in Creighton data
## Last updated: December 6, 2018

rm(list = ls())
a1 <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))


library(BayesMendel)
library(dplyr)

load("creighton_clean.RData")
source("MMRpro.cp.R")


expect <- function(fam, penet.m, penet.f){
  ex <- matrix(NA, nrow(fam), ncol(penet.m))
  for(i in 1:nrow(fam)){
    for(j in 1:ncol(penet.m)){
      ex[i, j] <- ifelse(fam$Gender[i] == 1, sum(penet.m[, j]), sum(penet.f[, j]))
    }
  }
  return(ex)
}


observe <- function(fam, penet.m, penet.f){
  obs <- matrix(NA, nrow(fam), ncol(penet.m))
  for(i in 1:nrow(fam)){
    for(j in 1:ncol(penet.m)){
      obs[i, j] <- ifelse(fam$AffectedColon[i] == 1, 1,
                          ifelse(fam$Gender[i] == 1, sum(penet.m[(fam$AgeColon[i] + 1):94, j]) / (1 - sum(penet.m[1:fam$AgeColon[i], j])),
                                 sum(penet.f[(fam$AgeColon[i] + 1):94, j]) / (1 - sum(penet.f[1:fam$AgeColon[i], j]))))
    }
  }
  return(obs)
}



fam.new <- function(fam){
  dat.new <- fam
  rel <- dat.new$Relation
  id.rel <- dat.new$ID
  id.pro <- dat.new$ID[dat.new$Relation == 1]
  dat.new$ethnic <- rep("Lynch", nrow(fam))
  ## imputing missing ages
  dat.new <- ImputeAge(
    fff = CheckFamStructure(model = "MMRpro", fff = as.data.frame(dat.new), counselee.id = id.pro,
                            germline.testing = NULL, marker.testing = NULL,
                            oophorectomy = NULL, mastectomy = NULL,
                            imputeAges = TRUE, imputeRelatives = TRUE,
                            params = MMRparams()),
    params = MMRparams(), model = "MMRpro"
  )$fff
  ## CheckFamStructure will remove the some columns if it needs to add relatives
  ## Here we add them back in because they are used in the estimating equation function
  ## and in brcapro
  dat.new$Relation <- 999; dat.new$Relation[dat.new$ID %in% id.rel] <- rel
  return(dat.new)
}






penet.m <- penet.mmr.net$fMX
penet.m <- penet.m + 1e-5
if(any(colSums(penet.m) > 1)){
  penet.m[, colSums(penet.m) > 1] <- penet.m[, colSums(penet.m) > 1] / colSums(penet.m)[colSums(penet.m) > 1]
}
penet.f <- penet.mmr.net$fFX
penet.f <- penet.f + 1e-5
if(any(colSums(penet.f) > 1)){
  penet.f[, colSums(penet.f) > 1] <- penet.f[, colSums(penet.f) > 1] / colSums(penet.f)[colSums(penet.f) > 1]
}


fam <- fam.new(mutate(filter(creighton, Famid == unique(creighton$Famid)[a1], !is.na(Gender)),
                      Twins = 0))
cp <- MMRpro.cp(fam, filter(fam, Relation == 1)$ID)
res.exp <- expect(fam, penet.m, penet.f)
res.obs <- observe(fam, penet.m, penet.f)
res.eo <- c(fam$Famid[1],
            sum(res.exp[, -1] * cp[, -1]),
            sum(res.exp[, 1] * cp[, 1]),
            sum(res.obs[, -1] * cp[, -1]),
            sum(res.obs[, 1] * cp[, 1]))

save(res.eo, file = paste(getwd(), "/Extended Frailty/OE Results/OE_Creighton_", a1, ".RData", sep = ""))
                        
                        
                        
                        
                        
                        