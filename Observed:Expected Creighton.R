## Looking at observed/expected ratios in Creighton data
## Last updated: December 5, 2018

library(BayesMendel)
library(dplyr)
library(data.table)
library(ggplot2)

load("creighton_clean.RData")
# creighton$Gender[creighton$Famid == 19 & creighton$ID == 19000377] <- 0

source("MMRpro.cp.R")


expect <- function(fam, penet.m, penet.f){
  # geno <- select(fam, MLH1, MSH2, MSH6)
  ex <- matrix(NA, nrow(fam), ncol(penet.m))
  for(i in 1:nrow(fam)){
    for(j in 1:ncol(penet.m)){
      # geno.i <- geno[i, ]
      # geno.j <- as.numeric(substring(colnames(penet.m)[j], 2:4, 2:4))
      # pos <- which(geno.i == 1)
      # neg <- which(geno.i == 2)
      # if(any(geno.j[pos] == 0) | any(geno.j[neg] == 1)){
      #   ex[i, j] <- 0
      # } else{
      #   ex[i,]
      # }
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


# exp.c <- function(fam, penet.m, penet.f, cp){
#   ex <- expect(fam, penet.m, penet.f)
#   sum(ex[, -1] * cp[, -1])
# }
# 
# exp.nc <- function(fam, penet.m, penet.f, cp){
#   ex <- expect(fam, penet.m, penet.f)
#   sum(ex[, 1] * cp[, 1])
# }
# 
# obs.c <- function(fam, penet.m, penet.f, cp){
#   obs <- observe(fam, penet.m, penet.f)
#   sum(obs[, -1] * cp[, -1])
# }
# 
# obs.nc <- function(fam, penet.m, penet.f, cp){
#   obs <- observe(fam, penet.m, penet.f)
#   sum(obs[, 1] * cp[, 1])
# }

nfam <- length(unique(creighton$Famid))

res.eo <- setNames(data.frame(matrix(NA, nfam, 5)),
                   c("Famid", "Exp.C", "Obs.C", "Exp.NC", "Obs.NC"))
res.eo$Famid <- unique(creighton$Famid)
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
cp <- res.exp <- res.obs <- vector("list", nfam)
for(i in 1:nfam){
  print(i)
  fam <- fam.new(mutate(filter(creighton, Famid == unique(creighton$Famid)[i], !is.na(Gender)),
                        Twins = 0))
  cp[[i]] <- MMRpro.cp(fam, filter(fam, Relation == 1)$ID)
  res.exp[[i]] <- expect(fam, penet.m, penet.f)
  res.obs[[i]] <- observe(fam, penet.m, penet.f)
  res.eo$Exp.C[i] <- sum(res.exp[[i]][, -1] * cp[[i]][, -1])
  res.eo$Exp.NC[i] <- sum(res.exp[[i]][, 1] * cp[[i]][, 1])
  res.eo$Obs.C[i] <- sum(res.obs[[i]][, -1] * cp[[i]][, -1])
  res.eo$Obs.NC[i] <- sum(res.obs[[i]][, 1] * cp[[i]][, 1])
}


res.eo.tot <- setNames(data.frame(matrix(NA, nfam, 9)),
                   c("Famid", "Obs.C", "Exp.C", "Obs.NC", "Exp.NC",
                     "OE.C", "OE.NC", "FamSize", "PropMMR"))
cp.tot <- res.exp.tot <- res.obs.tot <- vector("list", nfam)
for(i in setdiff(1:nfam, c(12, 55))){
  print(i)
  load(paste(getwd(), "/OE Results/OE_Creighton_", i, ".RData", sep = ""))
  # cp.tot[[i]] <- cp
  # res.exp.tot[[i]] <- res.exp
  # res.obs.tot[[i]] <- res.obs
  res.eo.tot[i, 1:5] <- res.eo[c(1, 4, 2, 5, 3)]
}
res.eo.tot$OE.C <- res.eo.tot$Obs.C / res.eo.tot$Exp.C
res.eo.tot$OE.NC <- res.eo.tot$Obs.NC / res.eo.tot$Exp.NC
dt <- data.table(creighton); dt <- dt[, nrow(.SD), by = Famid]; res.eo.tot$FamSize <- dt$V1
dt <- data.table(creighton); dt <- dt[, sum(.SD$MMR) / nrow(.SD), by = Famid]; res.eo.tot$PropMMR <- dt$V1



plot(res.eo.tot$OE.C, res.eo.tot$OE.NC,
     xlab = "O/E Carriers", ylab = "O/E Non-carriers")
abline(0,1); abline(h = 1); abline(v = 1)



ggplot(dplyr::arrange(res.eo.tot, PropMMR), aes(OE.C, OE.NC)) +
  geom_point(aes(color = PropMMR), size = 3) +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = 1) + 
  geom_vline(xintercept = 1) +
  labs(x = "O/E, Carriers", y = "O/E, Non-carriers")
cor(res.eo.tot$OE.C, res.eo.tot$OE.NC)

ggplot(res.eo.tot, aes(PropMMR, OE.C)) +
  geom_point(size = 3) +
  labs(x = "Proportion of family members who are MMR carriers",
       y = "O/E, Carriers")
cor(res.eo.tot$PropMMR, res.eo.tot$OE.C)

ggplot(res.eo.tot, aes(PropMMR, OE.NC)) +
  geom_point(size = 3) +
  labs(x = "Proportion of family members who are MMR carriers",
       y = "O/E, Non-arriers")
cor(res.eo.tot$PropMMR, res.eo.tot$OE.NC)

ggplot(res.eo.tot, aes(FamSize, OE.C)) +
  geom_point(size = 3) +
  labs(x = "Family size",
       y = "O/E, Carriers")
cor(res.eo.tot$FamSize, res.eo.tot$OE.C)

ggplot(res.eo.tot, aes(FamSize, OE.NC)) +
  geom_point(size = 3) +
  labs(x = "Family size",
       y = "O/E, Non-carriers")
cor(res.eo.tot$FamSize, res.eo.tot$OE.NC)
