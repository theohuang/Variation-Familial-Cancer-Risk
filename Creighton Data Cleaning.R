## Creighton data cleaning
## Last updated: December 13, 2018

library(dplyr)
library(data.table)
library(doParallel)
registerDoParallel(cores = 4)

## loading data
dir.ct <- "/Users/Theo/Dropbox (Partners HealthCare)/BayesMendel/Data/Colon-SPORE/Creighton"
load(paste(dir.ct, "/Creighton.RData", sep = ""))
dat <- mmrpro.mat

family.multiplier <- 1000000
# Famid <- round(family$id/family.multiplier)
# Source <- rep("Creighton", nrow(family))
# MSH6 <- MSI <- rep(0, nrow(family))
# FamilyName <- paste("Creighton",Famid, sep="")
# MLH1 <- family$MLH1
# MSH2 <- family$MSH2

individuals <- read.table("/Users/Theo/Dropbox (Partners HealthCare)/BayesMendel/Projects/MMRPRO_Validation/Sources/Creighton/data/individuals1214.csv",
                          header = TRUE, sep = ",", na.strings = "")
individuals$Result <- as.character(individuals$Result)
temp <- !is.na(individuals$Result) & individuals$Result=='n'
individuals$Result[temp] <- 'N'
individuals$gender <- as.character(individuals$gender)
individuals$ID <- individuals$family*family.multiplier + individuals$individ
individuals$fath.id <- individuals$family*family.multiplier + individuals$FATH.NUM
individuals$moth.id <- individuals$family*family.multiplier + individuals$MOTH.NUM


diagnoses <- read.table("/Users/Theo/Dropbox (Partners HealthCare)/BayesMendel/Projects/MMRPRO_Validation/Sources/Creighton/data/diagnoses1214.csv",
                        header = TRUE, sep = ",", na.strings = "")
family70 <- diagnoses[diagnoses$family==70,c("family","individ","end")]
family70$end <- NA
for(i in 1:length(family70$family))
{
  # end stores the new family id
  if(family70$individ[i]!=244) # no corresponding record for 244
  {
    family70$end[i] <- individuals$family[individuals$family>700 & individuals$individ==family70$individ[i]]
  }else{
    family70$end[i] <- 70
  }
}
diagnoses$family.old <- diagnoses$family
diagnoses$family[diagnoses$family.old==70] <- family70$end
diagnoses$ID <- diagnoses$family*family.multiplier + diagnoses$individ

# diagnoses$Dx1st indicates 1st Dx, but be aware there're people with multiple dx at the same time.
# don't waste time to find the smart way, just write a stupid loop! loooooooooooooop!

diagnoses$Dx1st <- 0
min.age <- rep(NA,length(unique(diagnoses$ID)))

for(i.ID in 1:length(unique(diagnoses$ID)))
{
  min.age[i.ID] <- min(diagnoses$DxAge[diagnoses$ID==unique(diagnoses$ID)[i.ID]],na.rm=T)
  if(min.age[i.ID]!=Inf)
  {
    diagnoses$Dx1st[diagnoses$ID==unique(diagnoses$ID)[i.ID] & diagnoses$DxAge==min.age[i.ID] & !is.na(diagnoses$DxAge)] <- 1
  }
}


## fixing ages of diagnosis and current ages
start <- Sys.time()
for(i in 1:nrow(dat)){
  if(dat$AffectedColon[i] == 1){
    dat$AgeColon[i] <- min(filter(diagnoses, family == dat$Famid[i], ID == dat$ID[i],
                                  SITECODE == 1, INVASIVE == 4)$DxAge)
  } else{
    dat$AgeColon[i] <- min(select(filter(individuals, family == dat$Famid[i], ID == dat$ID[i]),
                                  agedeath, agelfu), na.rm = TRUE)
  }
  if(dat$AffectedEndometrium[i] == 1){
    dat$AgeEndometrium[i] <- min(filter(diagnoses, family == dat$Famid[i], ID == dat$ID[i],
                                        SITECODE == 3, INVASIVE == 4)$DxAge)
  } else{
    dat$AgeEndometrium[i] <- min(select(filter(individuals, family == dat$Famid[i], ID == dat$ID[i]),
                                        agedeath, agelfu), na.rm = TRUE)
  }
}
difftime(Sys.time(), start, units = "secs")

dat <- mutate(dat, AgeColon = replace(AgeColon, AgeColon == Inf, NA),
              AgeEndometrium = replace(AgeEndometrium, AgeEndometrium == Inf, NA))
dat <- mutate(dat, AgeColon = floor(AgeColon), AgeEndometrium = floor(AgeEndometrium))

dat$MMR <- ifelse(dat$MLH1 == 1 | dat$MSH2 == 1 | dat$MSH6 == 1, 1,
                  ifelse(dat$MLH1 == 2 & dat$MSH2 == 2 & dat$MSH6 == 2, 2, 0))


# for(i in 1:length(which(dat$AffectedEndometrium == 1))){
#   ind.i <- which(dat$AffectedEndometrium == 1)[i]
#   dat$AgeEndometrium[ind.i] <-
#     min(filter(diagnoses, family == dat$Famid[ind.i], ID == dat$ID[ind.i],
#                SITECODE == 3, INVASIVE == 4)$DxAge)
# }

### Adding other cancers
dat <- mutate(dat,
              AffectedBreast = 0, AgeBreast = NA,
              AffectedGastric = 0, AgeGastric = NA,
              AffectedOvary = 0, AgeOvary = NA,
              AffectedSB = 0, AgeSB = NA,
              AffectedPancreas = 0, AgePancreas = NA,
              AffectedProstate = 0, AgeProstate = NA,
              AffectedUT = 0, AgeUT = NA,
              AffectedLiver = 0, AgeLiver = NA,
              AffectedKidney = 0, AgeKidney = NA,
              AffectedBD = 0, AgeBD = NA)
start <- Sys.time()
for(i in 1:nrow(dat)){
  if(i %% 1000 == 0) print(i)
  # breast
  br <- filter(diagnoses, family == dat$Famid[i], ID == dat$ID[i],
         SITECODE == 50, INVASIVE == 4)
  if(nrow(br) > 0){
    dat$AffectedBreast[i] <- 1
    dat$AgeBreast[i] <- min(br$DxAge)
  } else{
    dat$AgeBreast[i] <- ifelse(dat$AffectedEndometrium[i] == 0, dat$AgeEndometrium[i],
                               min(select(filter(individuals, family == dat$Famid[i], ID == dat$ID[i]),
                                   agedeath, agelfu), na.rm = TRUE))
  }
  
  # ovarian
  ov <- filter(diagnoses, family == dat$Famid[i], ID == dat$ID[i],
               SITECODE == 35, INVASIVE == 4)
  if(nrow(ov) > 0){
    dat$AffectedOvary[i] <- 1
    dat$AgeOvary[i] <- min(ov$DxAge)
  } else{
    dat$AgeOvary[i] <- ifelse(dat$AffectedEndometrium[i] == 0, dat$AgeEndometrium[i],
                               min(select(filter(individuals, family == dat$Famid[i], ID == dat$ID[i]),
                                          agedeath, agelfu), na.rm = TRUE))
  }
  
  # liver
  lv <- filter(diagnoses, family == dat$Famid[i], ID == dat$ID[i],
               SITECODE == 9, INVASIVE == 4)
  if(nrow(lv) > 0){
    dat$AffectedLiver[i] <- 1
    dat$AgeLiver[i] <- min(lv$DxAge)
  } else{
    dat$AgeLiver[i] <- ifelse(dat$AffectedEndometrium[i] == 0, dat$AgeEndometrium[i],
                              min(select(filter(individuals, family == dat$Famid[i], ID == dat$ID[i]),
                                         agedeath, agelfu), na.rm = TRUE))
  }
  
  # gastric
  gs <- filter(diagnoses, family == dat$Famid[i], ID == dat$ID[i],
               SITECODE == 11, INVASIVE == 4)
  if(nrow(gs) > 0){
    dat$AffectedGastric[i] <- 1
    dat$AgeGastric[i] <- min(gs$DxAge)
  } else{
    dat$AgeGastric[i] <- ifelse(dat$AffectedEndometrium[i] == 0, dat$AgeEndometrium[i],
                              min(select(filter(individuals, family == dat$Famid[i], ID == dat$ID[i]),
                                         agedeath, agelfu), na.rm = TRUE))
  }
  
  # kidney
  kd <- filter(diagnoses, family == dat$Famid[i], ID == dat$ID[i],
               SITECODE %in% 24:25, INVASIVE == 4)
  if(nrow(kd) > 0){
    dat$AffectedKidney[i] <- 1
    dat$AgeKidney[i] <- min(kd$DxAge)
  } else{
    dat$AgeKidney[i] <- ifelse(dat$AffectedEndometrium[i] == 0, dat$AgeEndometrium[i],
                              min(select(filter(individuals, family == dat$Famid[i], ID == dat$ID[i]),
                                         agedeath, agelfu), na.rm = TRUE))
  }
  
  # prostate
  pr <- filter(diagnoses, family == dat$Famid[i], ID == dat$ID[i],
               SITECODE == 6, INVASIVE == 4)
  if(nrow(pr) > 0){
    dat$AffectedProstate[i] <- 1
    dat$AgeProstate[i] <- min(pr$DxAge)
  } else{
    dat$AgeProstate[i] <- ifelse(dat$AffectedEndometrium[i] == 0, dat$AgeEndometrium[i],
                              min(select(filter(individuals, family == dat$Famid[i], ID == dat$ID[i]),
                                         agedeath, agelfu), na.rm = TRUE))
  }
  
  # pancreas
  pc <- filter(diagnoses, family == dat$Famid[i], ID == dat$ID[i],
               SITECODE == 18, INVASIVE == 4)
  if(nrow(pc) > 0){
    dat$AffectedPancreas[i] <- 1
    dat$AgePancreas[i] <- min(pc$DxAge)
  } else{
    dat$AgePancreas[i] <- ifelse(dat$AffectedEndometrium[i] == 0, dat$AgeEndometrium[i],
                              min(select(filter(individuals, family == dat$Famid[i], ID == dat$ID[i]),
                                         agedeath, agelfu), na.rm = TRUE))
  }
  
  # small bowel
  sb <- filter(diagnoses, family == dat$Famid[i], ID == dat$ID[i],
               SITECODE == 88, INVASIVE == 4)
  if(nrow(sb) > 0){
    dat$AffectedSB[i] <- 1
    dat$AgeSB[i] <- min(sb$DxAge)
  } else{
    dat$AgeSB[i] <- ifelse(dat$AffectedEndometrium[i] == 0, dat$AgeEndometrium[i],
                              min(select(filter(individuals, family == dat$Famid[i], ID == dat$ID[i]),
                                         agedeath, agelfu), na.rm = TRUE))
  }
  
  # bile duct
  bd <- filter(diagnoses, family == dat$Famid[i], ID == dat$ID[i],
               SITECODE == 43, INVASIVE == 4)
  if(nrow(bd) > 0){
    dat$AffectedBD[i] <- 1
    dat$AgeBD[i] <- min(bd$DxAge)
  } else{
    dat$AgeBD[i] <- ifelse(dat$AffectedEndometrium[i] == 0, dat$AgeEndometrium[i],
                              min(select(filter(individuals, family == dat$Famid[i], ID == dat$ID[i]),
                                         agedeath, agelfu), na.rm = TRUE))
  }
  
  # urinary tract
  ut <- filter(diagnoses, family == dat$Famid[i], ID == dat$ID[i],
               SITECODE == 30, INVASIVE == 4)
  if(nrow(ut) > 0){
    dat$AffectedUT[i] <- 1
    dat$AgeUT[i] <- min(ut$DxAge)
  } else{
    dat$AgeUT[i] <- ifelse(dat$AffectedEndometrium[i] == 0, dat$AgeEndometrium[i],
                              min(select(filter(individuals, family == dat$Famid[i], ID == dat$ID[i]),
                                         agedeath, agelfu), na.rm = TRUE))
  }
}
difftime(Sys.time(), start, units = "secs")






creighton <- dat
## fixing one mother whose gender is male
creighton$Gender[creighton$Famid == 19 & creighton$ID == 19000377] <- 0




save(creighton, file = "creighton_clean.RData")

