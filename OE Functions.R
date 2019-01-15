## O/E Functions
## Last updated: January 3, 2019

expect <- function(fam, penet.m, penet.f){
  ex <- matrix(NA, nrow(fam), ncol(penet.m))
  for(i in 1:nrow(fam)){
    for(j in 1:ncol(penet.m)){
      ex[i, j] <- ifelse(fam$Gender[i] == 1, sum(penet.m[, j]), sum(penet.f[, j]))
    }
  }
  return(ex)
}


observe <- function(fam, penet.m, penet.f){
  obs <- matrix(NA, nrow(fam), ncol(penet.m))
  for(i in 1:nrow(fam)){
    for(j in 1:ncol(penet.m)){
      if(is.na(fam$AffectedColon)[i]){
        obs[i, j] <- 0
      } else{
        if(fam$AffectedColon[i] == 1){
          obs[i, j] <- 1
        } else{
          if(fam$Gender[i] == 1){
            if(sum(penet.m[1:fam$AgeColon[i], j]) == 1){
              obs[i, j] <- 0
            } else{
              obs[i, j] <- ifelse(fam$AgeColon[i] == 94, 0, sum(penet.m[(fam$AgeColon[i] + 1):94, j]) / (1 - sum(penet.m[1:fam$AgeColon[i], j])))
            }
          } else{
            if(sum(penet.f[1:fam$AgeColon[i], j]) == 1){
              obs[i, j] <- 0
            } else{
              obs[i, j] <- ifelse(fam$AgeColon[i] == 94, 0, sum(penet.f[(fam$AgeColon[i] + 1):94, j]) / (1 - sum(penet.f[1:fam$AgeColon[i], j])))
            }
          }
        }
      }
      # obs[i, j] <- ifelse(fam$AffectedColon[i] == 1, 1,
      #                     ifelse(fam$Gender[i] == 1,
      #                            ifelse(sum(penet.m[1:fam$AgeColon[i], j]) == 1, 0,
      #                                   ifelse(sum(penet.m[(fam$AgeColon[i] + 1):94, j]) / (1 - sum(penet.m[1:fam$AgeColon[i], j])),
      #                                   ifelse(sum(penet.f[1:fam$AgeColon[i], j]) == 1, 0,
      #                                          sum(penet.f[(fam$AgeColon[i] + 1):94, j]) / (1 - sum(penet.f[1:fam$AgeColon[i], j]))))))
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


oe.boot.cp <- function(fam, penet.m, penet.f, nboot, seed = NULL){
  ## bootstrapping where we obtain the carrier probabilities for each
  ## bootstrap replicate of the family
  if(!is.null(seed)){
    set.seed(seed)
  }
  res.eo <- setNames(data.frame(matrix(NA, nboot, 5)),
                     c("FamID", "Obs.C", "Exp.C", "Obs.NC", "Exp.NC"))
  for(i in 1:nboot){
    print(i)
    ind.smp <- sample(1:nrow(fam), nrow(fam), replace = TRUE)
    id.nsmp <- setdiff(fam$ID, unique(fam$ID[ind.smp]))
    which.nsmp <- fam$ID %in% id.nsmp
    weights <- rep(0, nrow(fam))
    weights[sort(unique(ind.smp))] <- table(ind.smp)
    fam.boot <- fam
    fam.boot <- mutate(fam.boot,
                       AffectedColon = ifelse(which.nsmp, 0, AffectedColon),
                       AffectedEndometrium = ifelse(which.nsmp, 0, AffectedEndometrium),
                       AgeColon = ifelse(which.nsmp, 1, AgeColon),
                       AgeEndometrium = ifelse(which.nsmp, 1, AgeEndometrium),
                       MLH1 = ifelse(which.nsmp, 0, MLH1),
                       MSH2 = ifelse(which.nsmp, 0, MSH2),
                       MSH6 = ifelse(which.nsmp, 0, MSH6))
    cp.boot <- MMRpro.cp(fam.boot, fam.boot$ID[1])
    res.exp.boot <- expect(fam.boot, penet.m, penet.f)
    res.obs.boot <- observe(fam.boot, penet.m, penet.f)
    res.eo[i, ] <- c(fam$Famid[1],
                     sum(apply(res.obs.boot[, -1] * cp.boot[, -1], 2, function(x) x * weights)),
                     sum(apply(res.exp.boot[, -1] * cp.boot[, -1], 2, function(x) x * weights)),
                     sum(res.obs.boot[, 1] * cp.boot[, 1] * weights),
                     sum(res.exp.boot[, 1] * cp.boot[, 1] * weights))
  }
  return(res.eo)
}


oe.boot <- function(fam, penet.m, penet.f, nboot, seed = NULL){
  ## bootstrapping where we only obtain the carrier probabilities once
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  # getting the carrier probabilities
  cp <- MMRpro.cp(fam, filter(fam, Relation == 1)$ID)
  
  res.oe <- setNames(data.frame(matrix(NA, nboot, 5)),
                     c("FamID", "Obs.C", "Exp.C", "Obs.NC", "Exp.NC"))
  for(i in 1:nboot){
    if(i %% 100 == 0) print(i)
    ind.smp <- sample(1:nrow(fam), nrow(fam), replace = TRUE)
    weights <- rep(0, nrow(fam))
    weights[sort(unique(ind.smp))] <- table(ind.smp)
    res.exp.boot <- expect(fam, penet.m, penet.f)
    res.obs.boot <- observe(fam, penet.m, penet.f)
    res.oe[i, ] <- c(fam$Famid[1],
                     sum(apply(res.obs.boot[, -1] * cp[, -1], 2, function(x) x * weights)),
                     sum(apply(res.exp.boot[, -1] * cp[, -1], 2, function(x) x * weights)),
                     sum(res.obs.boot[, 1] * cp[, 1] * weights),
                     sum(res.exp.boot[, 1] * cp[, 1] * weights))
  }
  return(res.oe)
}



