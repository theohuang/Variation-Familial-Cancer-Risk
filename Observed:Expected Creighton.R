## Looking at observed/expected ratios in Creighton data
## Last updated: December 18, 2018

library(BayesMendel)
library(dplyr)
library(data.table)
library(ggplot2)
library(gplots)

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

nfam.creighton <- length(unique(creighton$Famid))

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
cp <- res.exp <- res.obs <- vector("list", nfam.creighton)
for(i in 1:nfam.creighton){
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


oe.dat <- setNames(data.frame(matrix(NA, nfam.creighton, 9)),
                   c("Famid", "Obs.C", "Exp.C", "Obs.NC", "Exp.NC",
                     "OE.C", "OE.NC", "FamSize", "PropMMR"))
cp.tot <- res.exp.tot <- res.obs.tot <- vector("list", nfam)
for(i in setdiff(1:nfam.creighton, c(12, 55))){
  if(file.exists(paste(getwd(), "/OE Results/OE_Creighton_", i, ".RData", sep = ""))){
    load(paste(getwd(), "/OE Results/OE_Creighton_", i, ".RData", sep = ""))
    # cp.tot[[i]] <- cp
    # res.exp.tot[[i]] <- res.exp
    # res.obs.tot[[i]] <- res.obs
    oe.dat[i, 1:5] <- res.eo[c(1, 4, 2, 5, 3)]
  }
}
oe.dat$OE.C <- oe.dat$Obs.C / oe.dat$Exp.C
oe.dat$OE.NC <- oe.dat$Obs.NC / oe.dat$Exp.NC
dt <- data.table(creighton); dt <- dt[, nrow(.SD), by = Famid]; oe.dat$FamSize <- dt$V1
dt <- data.table(creighton); dt <- dt[, sum(.SD$MMR) / nrow(.SD), by = Famid]; oe.dat$PropMMR <- dt$V1


oe.dat <- oe.dat[-c(12, 55), ]

cor(oe.dat$OE.C, oe.dat$OE.NC)


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




dt <- data.table(filter(creighton, Famid %in% oe.dat$Famid))
dt <- dt[, list(sum(.SD$AffectedEndometrium) / nrow(.SD),
                sum(.SD$AffectedBreast) / nrow(.SD),
                sum(.SD$AffectedOvary) / nrow(.SD),
                sum(.SD$AffectedGastric) / nrow(.SD),
                sum(.SD$AffectedSB) / nrow(.SD),
                sum(.SD$AffectedPancreas) / nrow(.SD),
                sum(.SD$AffectedProstate) / nrow(.SD),
                sum(.SD$AffectedUT) / nrow(.SD),
                sum(.SD$AffectedLiver) / nrow(.SD),
                sum(.SD$AffectedKidney) / nrow(.SD),
                sum(.SD$AffectedBD) / nrow(.SD),
                nrow(filter(.SD, AffectedEndometrium == 1 |
                              AffectedBreast == 1 |
                              AffectedOvary == 1 |
                              AffectedGastric == 1 |
                              AffectedSB == 1 |
                              AffectedPancreas == 1 |
                              AffectedProstate == 1 |
                              AffectedUT == 1 |
                              AffectedLiver == 1 |
                              AffectedKidney == 1 |
                              AffectedBD == 1)) / nrow(.SD)), by = Famid]
dt$V1

oe.dat$PropNonCRC <- dt$V12
oe.dat$LowessC <- lowess(dt$V12, oe.dat$OE.C)$y
oe.dat$LowessNC <- lowess(dt$V12, oe.dat$OE.NC)$y
oe.dat$LowessX <- lowess(dt$V12, oe.dat$OE.C)$x

save(oe.dat, file = "creighton_oe.RData")

plot(dt$V2, res.eo.tot$OE.C, xlab = "PropEC", ylab = "OE", col = "blue")
points(dt$V2, res.eo.tot$OE.NC, col = "red")
lines(lowess(dt$V2, res.eo.tot$OE.C)$x,
      lowess(dt$V2, res.eo.tot$OE.C)$y, col = "blue")
lines(lowess(dt$V2, res.eo.tot$OE.NC)$x,
      lowess(dt$V2, res.eo.tot$OE.NC)$y, col = "red")
abline(h = 1)

plot(dt$V12, res.eo.tot$OE.C, xlab = "PropNonCRC", ylab = "OE", col = "blue")
points(dt$V12, res.eo.tot$OE.NC, col = "red")
lines(lowess(dt$V12, res.eo.tot$OE.C)$x,
      lowess(dt$V12, res.eo.tot$OE.C)$y, col = "blue")
lines(lowess(dt$V12, res.eo.tot$OE.NC)$x,
      lowess(dt$V12, res.eo.tot$OE.NC)$y, col = "red")
abline(h = 1)


ggplot(data.frame(PropNonCRC = rep(res.eo.tot$PropNonCRC, 2),
                  OE = c(res.eo.tot$OE.C, res.eo.tot$OE.NC),
                  Carrier = rep(c("Yes", "No"), each = nrow(res.eo.tot)),
                  Lowess = c(res.eo.tot$LowessC, res.eo.tot$LowessNC),
                  LowessX = rep(res.eo.tot$LowessX, 2)),
       aes(PropNonCRC, OE)) +
  geom_point(aes(color = Carrier)) +
  geom_hline(yintercept = 1) +
  geom_line(aes(LowessX, Lowess, color = Carrier)) +
  labs(x = "Proportion of family with at least 1 non-CRC Lynch syndrome cancer",
       y = "O/E")




## bootstrap confidence interval results
res.ci2 <- setNames(data.frame(matrix(NA, 53, 7)),
                    c("Famid", "OE.C.mean", "OE.C.lo", "OE.C.hi",
                      "OE.NC.mean", "OE.NC.lo", "OE.NC.hi"))
res.ci2$Famid <- res.eo.tot$Famid
for(i in 1:53){
  if(file.exists(paste(getwd(), "/OE Results/Bootstrap/OE_Creighton_Boot_", i, ".RData", sep = ""))){
    load(paste(getwd(), "/OE Results/Bootstrap/OE_Creighton_Boot_", i, ".RData", sep = ""))
    res.eo <- mutate(res.eo, OE.C = Obs.C / Exp.C, OE.NC = Obs.NC / Exp.NC)
    res.ci2[res.ci2$Famid == res.eo$FamID[1], 2:7] <-
      c(mean(res.eo$OE.C), quantile(res.eo$OE.C, 0.025), quantile(res.eo$OE.C, 0.975),
        mean(res.eo$OE.NC), quantile(res.eo$OE.NC, 0.025), quantile(res.eo$OE.NC, 0.975))
  }
}
res.ci2 <- mutate(res.ci2, Ind = 1:53)
res.ci2 <- merge(res.ci2, select(res.eo.tot, Famid, PropMMR, PropNonCRC),
                 by = "Famid")

# plot(res.ci2$OE.C.mean, type = "l")
# lines(res.ci2$OE.C.lo, col = "blue")
# lines(res.ci2$OE.C.hi, col = "red")
# 
# 
# res.ci <- NULL
# for(i in 1:53){
#   if(file.exists(paste(getwd(), "/OE Results/Bootstrap/OE_Creighton_Boot_", i, ".RData", sep = ""))){
#     load(paste(getwd(), "/OE Results/Bootstrap/OE_Creighton_Boot_", i, ".RData", sep = ""))
#     res.eo <- mutate(res.eo, OE.C = Obs.C / Exp.C, OE.NC = Obs.NC / Exp.NC)
#     if(is.null(res.ci)){
#       res.ci <- res.eo
#     } else{
#       res.ci <- rbind(res.ci, res.eo)
#     }
#   }
# }
# res.ci <- mutate(res.ci, Ind = rep(1:length(unique(res.ci$FamID)), each = 50))


plotmeans(OE.C ~ Ind, data = res.ci, connect = FALSE, n.label = FALSE)

plotCI(res.ci2$OE.C.mean, uiw = res.ci2$OE.C.hi, liw = res.ci2$OE.C.lo)

ggplot(res.ci2, aes(PropNonCRC, OE.C.mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = OE.C.lo, ymax = OE.C.hi), width = 0.002) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Proportion of family members with at least 1 non-CRC Lynch syndrome cancer",
       y = "O/E, carriers")

ggplot(res.ci2, aes(PropNonCRC, OE.NC.mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = OE.NC.lo, ymax = OE.NC.hi), width = 0.002) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Proportion of family members with at least 1 non-CRC Lynch syndrome cancer",
       y = "O/E, non-carriers")

ggplot(res.ci2, aes(PropMMR, OE.C.mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = OE.C.lo, ymax = OE.C.hi), width = 0.002) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Proportion of family members with at least 1 MMR mutation",
       y = "O/E, carriers")

ggplot(res.ci2, aes(PropMMR, OE.NC.mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = OE.NC.lo, ymax = OE.NC.hi), width = 0.002) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Proportion of family members with at least 1 MMR mutation",
       y = "O/E, non-carriers")



