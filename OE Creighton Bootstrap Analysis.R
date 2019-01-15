## O/E Creighton Bootstrap Analysis
## Last updated: January 3, 2019

library(BayesMendel)
library(dplyr)
library(data.table)
library(ggplot2)
library(gplots)
library(tidyr)

load("creighton_clean.RData")
load("creighton_oe.RData")

source("MMRpro.cp.R")

famid.list <- setdiff(unique(creighton$Famid), c(19, 704))
nfam <- length(famid.list)

## getting family size, proportion of members who are MMR carriers,
## and proportion of members who have a non-CRC Lynch syndrome cancer
dt.fs <- data.table(filter(creighton, Famid %in% famid.list))
dt.fs <- dt.fs[, nrow(.SD), by = Famid]
dt.pm <- data.table(filter(creighton, Famid %in% famid.list))
dt.pm <- dt.pm[, sum(.SD$MMR) / nrow(.SD), by = Famid]
dt.c <- data.table(filter(creighton, Famid %in% res.eo.tot$Famid))
dt.c <- dt.c[, list(sum(.SD$AffectedEndometrium) / nrow(.SD),
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



## getting results
res.oe.tot <- setNames(vector("list", nfam), famid.list)
for(i in 1:nfam){
  load(paste(getwd(), "/OE Results/Bootstrap/OE_Creighton_Boot_", i, ".RData", sep = ""))
  res.oe <- mutate(res.oe, OE.C = Obs.C / Exp.C, OE.NC = Obs.NC / Exp.NC,
                   FamSize = filter(dt.fs, Famid == famid.list[i])$V1,
                   PropMMR = filter(dt.pm, Famid == famid.list[i])$V1,
                   PropNonCRC = filter(dt.c, Famid == famid.list[i])$V12)
  res.oe.tot[[i]] <- res.oe
}


## bootstrap confidence intervals for O/E ratios
ci.oe <- setNames(data.frame(matrix(NA, nfam, 7)),
                  c("Famid", "OE.C.mean", "OE.C.lo", "OE.C.hi",
                    "OE.NC.mean", "OE.NC.lo", "OE.NC.hi"))
ci.oe$Famid <- famid.list
for(i in 1:53){
  load(paste(getwd(), "/OE Results/Bootstrap/OE_Creighton_Boot_", i, ".RData", sep = ""))
  ci.oe[i, 2:7] <-
    c(mean(res.oe.tot[[i]]$OE.C), quantile(res.oe.tot[[i]]$OE.C, 0.025), quantile(res.oe.tot[[i]]$OE.C, 0.975),
      mean(res.oe.tot[[i]]$OE.NC), quantile(res.oe.tot[[i]]$OE.NC, 0.025), quantile(res.oe.tot[[i]]$OE.NC, 0.975))
}

## adding in information from the whole data
ci.oe <- merge(ci.oe, select(oe.dat, Famid, OE.C, OE.NC, FamSize, PropMMR, PropNonCRC), by = "Famid")
ci.oe <- mutate(ci.oe, Width.C = OE.C.hi - OE.C.lo,
                Width.NC = OE.NC.hi - OE.NC.lo)
dt <- data.table(filter(creighton, Famid %in% famid.list))
dt <- dt[, sum(.SD$AffectedColon) / nrow(.SD), by = Famid]
ci.oe$PropCRC <- dt$V1

gather(ci.oe, Data, Mean, c("OE.C", "OE.C.mean"))

## plots of confidence intervals for the O/E ratios
ggplot(ci.oe, aes(OE.NC, OE.C)) +
  geom_point() +
  geom_point(aes(OE.C.mean, OE.NC.mean), color = "blue", alpha = 0.5) +
  geom_errorbar(aes(ymin = OE.C.lo, ymax = OE.C.hi), width = 0.02, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "O/E, non-carriers",
       y = "O/E, carriers")

ggplot(ci.oe, aes(OE.C, OE.NC)) +
  geom_point() +
  geom_point(aes(OE.C.mean, OE.NC.mean), color = "blue", alpha = 0.5) +
  geom_errorbar(aes(ymin = OE.NC.lo, ymax = OE.NC.hi), width = 0.02, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "O/E, carriers",
       y = "O/E, non-carriers")

ggplot(ci.oe, aes(OE.C, OE.NC)) +
  geom_point() +
  geom_point(aes(OE.C.mean, OE.NC.mean), color = "blue", alpha = 0.5) +
  geom_errorbarh(aes(xmin = OE.C.lo, xmax = OE.C.hi), height = 0.02, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "O/E, carriers",
       y = "O/E, non-carriers")

ggplot(ci.oe, aes(OE.C, OE.NC)) +
  geom_point() +
  geom_point(aes(OE.C.mean, OE.NC.mean), color = "blue", alpha = 0.5) +
  geom_errorbar(aes(ymin = OE.NC.lo, ymax = OE.NC.hi), width = 0.02, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OE.C.lo, xmax = OE.C.hi), height = 0.03, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "O/E, carriers",
       y = "O/E, non-carriers")


# ggplot(gather(ci.oe, Data, Mean, c("OE.C", "OE.C.mean")), aes(Mean, OE.NC)) +
#   geom_point(aes(color = Data)) +
#   geom_errorbar(aes(ymin = OE.NC.lo, ymax = OE.NC.hi), width = 0.02, alpha = 0.5) +
#   geom_errorbarh(aes(xmin = OE.C.lo, xmax = OE.C.hi), height = 0.02, alpha = 0.5) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
#   geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
#   labs(x = "O/E, carriers",
#        y = "O/E, non-carriers")

plot(ci.oe$OE.C, ci.oe$OE.C.mean)
abline(0,1)

plot(ci.oe$OE.NC, ci.oe$OE.NC.mean)
abline(0,1)

ggplot(ci.oe, aes(Width.C, Width.NC)) +
  geom_point(aes(color = FamSize), size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Width of CI, carriers", y = "Width of CI, non-carriers")

ggplot(ci.oe, aes(Width.C, FamSize)) +
  geom_point(aes(color = Width.NC)) +
  labs(x = "Width of CI, carriers", y = "FamSize") +
  scale_color_continuous(name = "Width, NC")

ggplot(ci.oe, aes(Width.NC, FamSize)) +
  geom_point(aes(color = Width.C)) +
  labs(x = "Width of CI, non-carriers", y = "FamSize") +
  scale_color_continuous(name = "Width, C")

ggplot(ci.oe, aes(Width.C, Width.NC)) +
  geom_point(aes(color = PropNonCRC), size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Width of CI, carriers", y = "Width of CI, non-carriers")

ggplot(ci.oe, aes(Width.NC, PropNonCRC)) +
  geom_point(aes(color = Width.C)) +
  labs(x = "Width of CI, non-carriers", y = "Prop of family with non-CRC LS cancers") +
  scale_color_continuous(name = "Width, C")

ggplot(ci.oe, aes(Width.C, PropNonCRC)) +
  geom_point(aes(color = Width.NC)) +
  labs(x = "Width of CI, carriers", y = "Prop of family with non-CRC LS cancers") +
  scale_color_continuous(name = "Width, NC")

ggplot(ci.oe, aes(Width.NC, PropNonCRC)) +
  geom_point(aes(color = Width.C)) +
  labs(x = "Width of CI, non-carriers", y = "Prop of family with non-CRC LS cancers") +
  scale_color_continuous(name = "Width, C")

ggplot(ci.oe, aes(Width.C, Width.NC)) +
  geom_point(aes(color = PropCRC), size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Width of CI, carriers", y = "Width of CI, non-carriers")

ggplot(ci.oe, aes(Width.C, PropCRC)) +
  geom_point(aes(color = Width.NC)) +
  labs(x = "Width of CI, carriers", y = "Prop of family with CRC") +
  scale_color_continuous(name = "Width, NC")

ggplot(ci.oe, aes(Width.NC, PropCRC)) +
  geom_point(aes(color = Width.C)) +
  labs(x = "Width of CI, non-carriers", y = "Prop of family with CRC") +
  scale_color_continuous(name = "Width, C")

ggplot(ci.oe, aes(Width.C, Width.NC)) +
  geom_point(aes(color = PropMMR), size = 3) +
  labs(x = "Width of CI, carriers", y = "Width of CI, non-carriers") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

ggplot(ci.oe, aes(Width.C, PropMMR)) +
  geom_point(aes(color = Width.NC)) +
  labs(x = "Width of CI, carriers", y = "Prop of family who are carriers") +
  scale_color_continuous(name = "Width, NC")

ggplot(ci.oe, aes(Width.NC, PropMMR)) +
  geom_point(aes(color = Width.C)) +
  labs(x = "Width of CI, non-carriers", y = "Prop of family who are carriers") +
  scale_color_continuous(name = "Width, C")


cor(select(ci.oe, Width.C, FamSize))[1, 2]
cor(select(ci.oe, Width.NC, FamSize))[1, 2]
cor(select(ci.oe, Width.C, PropCRC))[1, 2]
cor(select(ci.oe, Width.NC, PropCRC))[1, 2]
cor(select(ci.oe, Width.C, PropNonCRC))[1, 2]
cor(select(ci.oe, Width.NC, PropNonCRC))[1, 2]
cor(select(ci.oe, Width.C, PropMMR))[1, 2]
cor(select(ci.oe, Width.NC, PropMMR))[1, 2]

plot(ci.oe$PropNonCRC, ci.oe$Width.)

### confidence interval for correlation between O/E for carriers and non-carriers
set.seed(999)
nboot.cor <- 1000
cor.boot <- rep(0, nboot.cor)
# famid.ct <- rep(0, nfam)
for(i in 1:nboot.cor){
  ind.boot <- sample(1:nfam, nfam, replace = TRUE)
  # famid.boot <- famid.list[sample(1:nfam, nfam, replace = TRUE)]
  # famid.boot.sort <- sort(unique(famid.boot))
  # oe.boot.c <- oe.boot.nc <- vector()
  # for(j in 1:length(famid.boot.sort)){
  #   nj <- length(which(famid.boot == famid.boot.sort[j]))
  #   for(k in 1:nj){
  #     oe.boot.c <- c(oe.boot.c, res.oe.tot[[as.character(famid.boot.sort[j])]]$OE.C[famid.ct[famid.list == famid.boot.sort[j]] + k])
  #     oe.boot.nc <- c(oe.boot.nc, res.oe.tot[[as.character(famid.boot.sort[j])]]$OE.NC[famid.ct[famid.list == famid.boot.sort[j]] + k])
  #   }
  #   famid.ct[famid.list == famid.boot.sort[j]] <- famid.ct[famid.list == famid.boot.sort[j]] + nj
  # }
  # cor.boot[i] <- cor(oe.boot.c, oe.boot.nc)
  cor.boot[i] <- cor(oe.dat$OE.C[ind.boot], oe.dat$OE.NC[ind.boot])
}


summary(cor.boot)
c(cor(oe.dat$OE.C, oe.dat$OE.NC), quantile(cor.boot, 0.025), quantile(cor.boot, 0.975))

hist(cor.boot, xlab = "Correlation between O/E for carriers and O/E for non-carriers",
     main = "")
abline(v = cor(oe.dat$OE.C, oe.dat$OE.NC), col = "red", lwd = 2)
abline(v = mean(cor.boot), col = "blue", lwd = 2)
legend("topright", c("Cor for data", "Mean cor of bootstrap sample"),
       col = c("red", "blue"), lwd = 2)


save(res.oe.tot, ci.oe, cor.boot, file = "creighton_boot.RData")

