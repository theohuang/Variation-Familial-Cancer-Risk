# Output from 'brcapro' and 'MMRpro' - family matrix and probabilities - used to plot pedigree structure
# using generic function 'plot'


plot.BayesMendel.id <- function(family) {
  
  if (is.null(family$Death)) {status=rep(0, nrow(family))}
  if (!is.null(family$Death)) {status=family$Death}	
  
  canc <- "colon"
  
  id <- family$ID
  dadid <- family$FatherID
  momid <- family$MotherID
  sex <- family$Gender
  sex[sex==0] <- 2
  
  refs <- 1:length(id)
  get1 <- refs[dadid==0 & momid > 0]
  get2 <- refs[dadid > 0 & momid==0]
  rem <- sort(c(get1, get2))
  
  num <- length(unique(momid[get1])) + length(unique(dadid[get2]))
  
  if (num > 0) {
    
    newid <- max(id)+(1:num)
    newid1 <- newid[1:length(unique(momid[get1]))]
    newid2 <- newid[(length(unique(momid[get1]))+1):(length(newid))]
    
    id <- c(id, newid)
    sex <- c(sex, rep(1, length(unique(momid[get1]))), rep(2, length(unique(dadid[get2]))))
    
    dadid <- c(dadid, rep(0, num))
    momid <- c(momid, rep(0, num))	
    
    if (length(get1) > 0) {
      for (i in 1:length(get1)) {
        dadid[id==id[get1[i]]] <- newid1[unique(momid[get1])==momid[get1[i]]]
      }
    }
    
    if (length(get2) > 0) {
      for (i in 1:length(get2)) {
        momid[id==id[get2[i]]] <- newid2[unique(dadid[get2])==dadid[get2[i]]]
      }
    }
    
  }
  
  affected <- cbind(family$AffectedColon, family$AffectedEndometrium)
  affected <- affected + 1
  affected[affected == 3] <- 2
  
  ages.cc <- c(family$AgeColon, rep(1, num))
  ages.ec <- c(family$AgeEndometrium, rep(1, num))
  cc <- family$AffectedColon
  ec <- family$AffectedEndometrium
  
  affected <- rbind(affected, matrix(1, nrow=num, ncol=2))
  status <- c(status, rep(0, num))
  relations <- NULL
  
  if (!is.null(family$Twins) & sum(family$Twins)>0) {
    twin.ids <- unique(family$Twins[family$Twins > 0])
    id1 <- id2 <- code <- NULL
    for (ttt in 1:length(twin.ids)) {
      id1 <- c(id1,family$ID[family$Twins==unique(twin.ids)[ttt]][1])
      id2 <- c(id2,family$ID[family$Twins==unique(twin.ids)[ttt]][2])
      code <- c(code, 1)
    }
    relations <- matrix(c(id1,id2,code),byrow=F,ncol=3)
    ped.y <- kinship2::pedigree(id, dadid, momid, sex, affected, status, relation=relations)
  }
  
  if ((!is.null(family$Twins) & sum(family$Twins)==0) | is.null(family$Twins)) {
    ped.y <- kinship2::pedigree(id, dadid, momid, sex, affected, status)
  }
  
  par(xpd=TRUE)
  plt <- kinship2::plot.pedigree(ped.y, id=id, align=FALSE)
  # xp <- plt$x
  # yp <- plt$y
  # 
  # inc <- (abs(min(yp)) - abs(max(yp)))/12
  # xp.inc <- (abs(min(xp)) - abs(max(xp)))/50
  # 
  # couns <- filter(fam, Relation == 1)$ID
  # xp.couns <- xp[id==couns]
  # yp.couns <- yp[id==couns]
  # arrows(xp.couns-.35, yp.couns-.15, xp.couns-0.1, yp.couns-0.1, length=0.1)
  # 
}

plot.BayesMendel.nid <- function(family) {
  
  if (is.null(family$Death)) {status=rep(0, nrow(family))}
  if (!is.null(family$Death)) {status=family$Death}	
  
  canc <- "colon"
  
  id <- family$ID
  dadid <- family$FatherID
  momid <- family$MotherID
  sex <- family$Gender
  sex[sex==0] <- 2
  
  refs <- 1:length(id)
  get1 <- refs[dadid==0 & momid > 0]
  get2 <- refs[dadid > 0 & momid==0]
  rem <- sort(c(get1, get2))
  
  num <- length(unique(momid[get1])) + length(unique(dadid[get2]))
  
  if (num > 0) {
    
    newid <- max(id)+(1:num)
    newid1 <- newid[1:length(unique(momid[get1]))]
    newid2 <- newid[(length(unique(momid[get1]))+1):(length(newid))]
    
    id <- c(id, newid)
    sex <- c(sex, rep(1, length(unique(momid[get1]))), rep(2, length(unique(dadid[get2]))))
    
    dadid <- c(dadid, rep(0, num))
    momid <- c(momid, rep(0, num))	
    
    if (length(get1) > 0) {
      for (i in 1:length(get1)) {
        dadid[id==id[get1[i]]] <- newid1[unique(momid[get1])==momid[get1[i]]]
      }
    }
    
    if (length(get2) > 0) {
      for (i in 1:length(get2)) {
        momid[id==id[get2[i]]] <- newid2[unique(dadid[get2])==dadid[get2[i]]]
      }
    }
    
  }
  
  affected <- cbind(family$AffectedColon, family$AffectedEndometrium)
  affected <- affected + 1
  affected[affected == 3] <- 2
  
  ages.cc <- c(family$AgeColon, rep(1, num))
  ages.ec <- c(family$AgeEndometrium, rep(1, num))
  cc <- family$AffectedColon
  ec <- family$AffectedEndometrium
  
  affected <- rbind(affected, matrix(1, nrow=num, ncol=2))
  status <- c(status, rep(0, num))
  relations <- NULL
  
  if (!is.null(family$Twins) & sum(family$Twins)>0) {
    twin.ids <- unique(family$Twins[family$Twins > 0])
    id1 <- id2 <- code <- NULL
    for (ttt in 1:length(twin.ids)) {
      id1 <- c(id1,family$ID[family$Twins==unique(twin.ids)[ttt]][1])
      id2 <- c(id2,family$ID[family$Twins==unique(twin.ids)[ttt]][2])
      code <- c(code, 1)
    }
    relations <- matrix(c(id1,id2,code),byrow=F,ncol=3)
    ped.y <- kinship2::pedigree(id, dadid, momid, sex, affected, status, relation=relations)
  }
  
  if ((!is.null(family$Twins) & sum(family$Twins)==0) | is.null(family$Twins)) {
    ped.y <- kinship2::pedigree(id, dadid, momid, sex, affected, status)
  }
  
  par(xpd=TRUE)
  
  df.col <- data.frame(ID = ped.y$id)
  df.col <- mutate(merge(df.col, select(family, ID, MLH1), by = "ID", all.x = TRUE),
                   MLH1 = ifelse(is.na(MLH1), 0, MLH1),
                   col = ifelse(MLH1 == 1, "red", ifelse(MLH1 == 2, "blue", "black")))
  
  ## plotting colors: red is carriers and blue is non-carriers. black is unknown
  plt <- kinship2::plot.pedigree(ped.y, id=rep("", nrow(ped.y)), align = FALSE,
                                 col = df.col$col)
  legend("topright", legend = expression(italic("MLH1")*" Carrier", italic("MLH1")*" Non-carrier", italic("MLH1")*" Carrier Status Unknown"),
         col = c("red", "blue", "black"), pch = 19)
  # xp <- plt$x
  # yp <- plt$y
  # 
  # inc <- (abs(min(yp)) - abs(max(yp)))/12
  # xp.inc <- (abs(min(xp)) - abs(max(xp)))/50
  # 
  # couns <- filter(fam, Relation == 1)$ID
  # xp.couns <- xp[id==couns]
  # yp.couns <- yp[id==couns]
  # arrows(xp.couns-.35, yp.couns-.15, xp.couns-0.1, yp.couns-0.1, length=0.1)
  
}


