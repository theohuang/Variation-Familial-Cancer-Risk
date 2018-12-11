## Exploring Creighton Data
## Last updated: December 3, 2018

library(dplyr)
library(data.table)

load("creighton_clean.RData")
famid.list <- unique(creighton$Famid)

## loading data
# dir.ct <- "/Users/Theo/Dropbox (Partners HealthCare)/BayesMendel/Data/Colon-SPORE/Creighton"
# load(paste(dir.ct, "/Creighton.RData", sep = ""))
# 
# load(paste(dir.ct, "/CreightonFull.RData", sep = ""))
# load("creighton_clean.RData")

# dat <- mmrpro.mat


# famid.list <- vector()
# for(i in 1:length(unique(dat$FamilyName))){
#   famid.list.i <- unique(dat$Famid[dat$FamilyName == as.character(unique(dat$FamilyName)[i])])
#   famid.list[i] <- famid.list.i[1]
# }
# dat <- filter(dat, Famid %in% famid.list)

creighton.pro <- filter(creighton, Relation == 1)


## number of families
nfam <- length(unique(creighton$Famid)) #55

## average family size
nrow(creighton) / nfam # 236.13
dt <- data.table(creighton)
dt <- dt[, nrow(.SD), by = Famid]
hist(dt$V1, xlab = "Family size", main = "", breaks = 15)

## proband gender
table(creighton.pro$Gender)

## number of probands with CRC
table(creighton.pro$AffectedColon) # 79

## number of probands with EC
table(creighton.pro$AffectedEndometrium) # 20

## proband genetic testing
table(creighton.pro$MLH1)
table(creighton.pro$MSH2)
table(creighton.pro$MSH6)
creighton$MMR <- ifelse(creighton$MLH1 == 1 | creighton$MSH2 == 1 | creighton$MSH6 == 1, 1,
                  ifelse(creighton$MLH1 == 2 & creighton$MSH2 == 2 & creighton$MSH6 == 2, 2, 0))
creighton.pro <- filter(creighton, Relation == 1)
table(creighton.pro$MMR) # 162 carriers, 296 non-carriers (all probands were tested)

## family member genetic testing
sum(filter(creighton, Relation != 1)$MMR) / nrow(filter(creighton, Relation != 1))

dt <- data.table(filter(creighton, Famid %in% famid.list[-c(12, 55)]))
dt <- data.table(creighton)
dt <- dt[, nrow(filter(.SD, MLH1 %in% 1:2 | MSH2 %in% 1:2 | MSH6 %in% 1:2)), by = Famid]
dt$V1
hist(dt$V1, xlab = "Number of family members with genetic testing", main = "")

dt <- data.table(filter(creighton, Famid %in% famid.list[-c(12, 55)]))
dt <- data.table(creighton)
dt <- dt[, nrow(filter(.SD, MLH1 %in% 1:2 | MSH2 %in% 1:2 | MSH6 %in% 1:2)) / nrow(.SD), by = Famid]
dt$V1
hist(dt$V1, xlab = "Proportion of family members with genetic testing", main = "")

dt <- data.table(filter(creighton, Famid %in% famid.list[-c(12, 55)]))
dt <- data.table(creighton)
dt <- dt[, nrow(filter(.SD, MLH1 == 1 | MSH2 == 1 | MSH6 == 1)), by = Famid]
dt$V1
hist(dt$V1, xlab = "Number of family members who are MMR carriers", main = "")

dt <- data.table(filter(creighton, Famid %in% famid.list[-c(12, 55)]))
dt <- data.table(creighton)
dt <- dt[, nrow(filter(.SD, MLH1 == 2 | MSH2 == 2 | MSH6 == 2)), by = Famid]
dt$V1
hist(dt$V1, xlab = "Number of family members who are non-carriers of at least 1 gene", main = "")


## average number of affected relatives
nrow(filter(creighton, Relation != 1, AffectedColon == 1)) / nfam
nrow(filter(creighton, Relation != 1, AffectedEndometrium == 1)) / nfam


## average proband age







