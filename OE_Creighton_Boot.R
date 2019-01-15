## Looking at observed/expected ratios in Creighton data
## Bootstrap confidence intervals
## Last updated: December 17, 2018

rm(list = ls())
a1 <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

library(BayesMendel)
library(dplyr)

load("creighton_clean.RData")
source("MMRpro.cp.R")
source("ImputeAge.cp.R")
source("OE Functions.R")

famid.list <- setdiff(unique(creighton$Famid), c(19, 704))

penet.m <- penet.mmr.net$fMX
penet.f <- penet.mmr.net$fFX



fam <- fam.new(mutate(filter(creighton, Famid == famid.list[a1], !is.na(Gender)),
                      Twins = 0))
nboot <- 1500
# cp <- MMRpro.cp(fam, filter(fam, Relation == 1)$ID)
# res.exp <- expect(fam, penet.m, penet.f)
# res.obs <- observe(fam, penet.m, penet.f)
# res.eo <- c(fam$Famid[1],
#             sum(res.exp[, -1] * cp[, -1]),
#             sum(res.exp[, 1] * cp[, 1]),
#             sum(res.obs[, -1] * cp[, -1]),
#             sum(res.obs[, 1] * cp[, 1]))

start <- Sys.time()
res.oe <- oe.boot(fam, penet.m, penet.f, nboot, seed = 999)
difftime(Sys.time(), start, units = "secs")

save(res.oe, file = paste(getwd(), "/Extended Frailty/OE Results/Bootstrap/OE_Creighton_Boot_", a1, ".RData", sep = ""))





