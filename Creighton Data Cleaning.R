## Creighton data cleaning
## Last updated: December 6, 2018

library(dplyr)
library(data.table)
library(doParallel)
registerDoParallel(cores = 4)

## loading data
dir.ct <- "/Users/Theo/Dropbox (Partners HealthCare)/BayesMendel/Data/Colon-SPORE/Creighton"
load(paste(dir.ct, "/Creighton.RData", sep = ""))
dat <- mmrpro.mat

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
                                        SITECODE %in% 3:4, INVASIVE == 4)$DxAge)
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



creighton <- dat
## fixing one mother whose gender is male
creighton$Gender[creighton$Famid == 19 & creighton$ID == 19000377] <- 0




save(creighton, file = "creighton_clean.RData")

