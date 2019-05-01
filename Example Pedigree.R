## Generating example pedigree
## Lsat updated: April 5, 2019

library(BayesMendel)
library(dplyr)

load("creighton_clean.RData")
source("MMRpro.cp.R")
source("ImputeAge.cp.R")
source("OE Functions.R")
source("plot.BayesMendel.id.R")

set.seed(12345)

fam <- filter(creighton, Famid == 24)
fam <- fam %>%
  mutate(Twins = 0, ID = ID - 2.4e7,
         MotherID = ifelse(MotherID != 0, MotherID - 2.4e7, 0),
         FatherID = ifelse(FatherID != 0, FatherID - 2.4e7, 0)) %>%
  mutate(FatherID = replace(FatherID, ID == 9, 0),
         MotherID = replace(MotherID, ID == 9, 0)) %>%
  mutate(Relation = ifelse(ID == 3, 1, ifelse(ID == 7, 0, Relation)))

fam <- add_row(fam, Famid = 24, ID = 100:136,
               FatherID = c(69, 69, 69, 3, 75, 72, 13, 32, 85, 100, 101, 102,
                            110, 110, 111, 31, 0, 0, 88, 88, 84, 84, 70, 56, 56, 45, 45, 0,
                            130, 130, 0, 67, 0, 111, 111, 100, 100),
               MotherID = c(68, 68, 68, 0, 6, 12, 0, 33, 86, 0, 0, 0,
                            116, 116, 117, 127, 0, 0, 89, 89, 83, 83, 71, 57, 57, 46, 46, 0,
                            65, 65, 0, 132, 0, 117, 117, 138, 138),
               Twins = 0, AffectedColon = 0, AffectedEndometrium = 0,
               AgeColon = NA, AgeEndometrium = NA,
               Gender = ifelse(ID %in% c(116:117, 127, 132), 0,
                               ifelse(ID %in% 110:111, 1,
                                      sample(0:1, length(100:136), replace = TRUE))))

fam <- fam %>% filter(!(ID %in% c(78, 81, 105, 90, 48, 49, 50, 51, 90, 26, 62, 87, 91,
                                  10, 11, 12, 13, 72, 113, 73, 105, 106, 85, 108))) %>%
  mutate(AffectedColon = ifelse(ID %in% c(69, 54, 101, 102, 63, 64, 52, 110,
                                          8, 19, 40, 1, 18, 37, 6, 60), 1, 0),
         AffectedEndometrium = ifelse(ID %in% c(40, 66, 4, 83), 1, 0)) %>% 
  mutate(AgeColon = ifelse(is.na(AgeColon) | AgeColon %in% 0:1,
                           sample(50:80, nrow(fam), replace = TRUE), 
                           ifelse(AffectedColon == 1,
                                  sample(40:70, nrow(fam), replace = TRUE), AgeColon)),
         AgeEndometrium = ifelse(is.na(AgeEndometrium) | AgeEndometrium %in% 0:1,
                                 sample(50:80, nrow(fam), replace = TRUE),
                                 ifelse(AffectedEndometrium == 1,
                                        sample(40:70, nrow(fam), replace = TRUE), AgeEndometrium))) %>%
  mutate(AgeEndometrium = ifelse(AffectedColon == 1 & AffectedEndometrium == 0 & AgeColon > AgeEndometrium,
                                 sample(50:70, nrow(fam), replace = TRUE), AgeColon),
         AgeColon = ifelse(AffectedColon == 1 & AffectedEndometrium == 0 & AgeColon > AgeEndometrium,
                           AgeEndometrium - 5, AgeColon),
         AgeColon = ifelse(AffectedColon == 0 & AffectedEndometrium == 1 & AgeColon < AgeEndometrium,
                           sample(50:70, nrow(fam), replace = TRUE), AgeEndometrium),
         AgeEndometrium = ifelse(AffectedColon == 0 & AffectedEndometrium == 1 & AgeColon < AgeEndometrium,
                                 AgeColon - 5, AgeEndometrium),
         AgeColon = ifelse(AffectedColon == 0 & AffectedEndometrium == 0 & AgeColon != AgeEndometrium,
                           AgeEndometrium, AgeColon))

plot.BayesMendel.id(fam)
plot.BayesMendel.nid(fam)


cp <- MMRpro.cp(fam, filter(fam, Relation == 1)$ID)
penet.m <- penet.mmr.net$fMX
penet.f <- penet.mmr.net$fFX



weights <- rep(1, nrow(fam))
res.exp <- expect(fam, penet.m, penet.f)
res.obs <- observe(fam, penet.m, penet.f)
res.oe <- c(sum(apply(res.obs[, -1] * cp[, -1], 2, function(x) x * weights)),
            sum(apply(res.exp[, -1] * cp[, -1], 2, function(x) x * weights)),
            sum(res.obs[, 1] * cp[, 1] * weights),
            sum(res.exp[, 1] * cp[, 1] * weights))
res.oe[1] / res.oe[2]; res.oe[3] / res.oe[4]

start <- Sys.time()
res.oe.boot <- oe.boot(fam, penet.m, penet.f, nboot = 1000, seed = 12345)
difftime(Sys.time(), start, units = "secs")

res.oe.boot <- mutate(res.oe.boot, OE.C = Obs.C / Exp.C, OE.NC = Obs.NC / Exp.NC)

c(quantile(res.oe.boot$OE.C, 0.025), quantile(res.oe.boot$OE.C, 0.975))
c(quantile(res.oe.boot$OE.NC, 0.025), quantile(res.oe.boot$OE.NC, 0.975))

save(fam, res.oe, res.oe.boot, file = "example_ped.RData")
