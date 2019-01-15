## Looking at observed/expected ratios in Creighton data
## Last updated: December 13, 2018

rm(list = ls())
a1 <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))


library(BayesMendel)
library(dplyr)

load("creighton_clean.RData")
source("MMRpro.cp.R")
source("OE Functions.R")




penet.m <- penet.mmr.net$fMX
# penet.m <- penet.m + 1e-5
# if(any(colSums(penet.m) > 1)){
#   penet.m[, colSums(penet.m) > 1] <- penet.m[, colSums(penet.m) > 1] / colSums(penet.m)[colSums(penet.m) > 1]
# }
penet.f <- penet.mmr.net$fFX
# penet.f <- penet.f + 1e-5
# if(any(colSums(penet.f) > 1)){
#   penet.f[, colSums(penet.f) > 1] <- penet.f[, colSums(penet.f) > 1] / colSums(penet.f)[colSums(penet.f) > 1]
# }


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





