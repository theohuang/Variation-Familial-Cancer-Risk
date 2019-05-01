#' Subset Family Matrix by Pattern Matching and by Cancer
#' 
#' @param pattern character string containing a regular expression to be matched 
#' in the column names of \code{fam}. 
#' @param fam a family matrix
#' @param cancers names of cancers to use
#' @details This is a helper function for the \code{estLik} function that was written
#' to pull out the \code{isAff} and \code{Age} variables for the cancers of interest. 
#' @family multigene
subsetCancers = function(pattern, fam, cancers) {
  # Pull out columns with names that match `pattern`
  matches = data.frame(fam[, grep(pattern, names(fam)), drop=FALSE])
  names(matches) = base::sub(pattern, "", names(matches))
  
  # Subset the `pattern`-matched columns for cancers of interest
  if (length(intersect(names(matches), cancers)) < length(cancers)) {
    stop("Family matrix is missing information for one or more 
         cancers in penetrance matrices.")
  } else {
    matches = matches[,cancers, drop=FALSE]
  }
}



#' Estimate Likelihood Matrix
#' 
#' The likelihood matrix is estimated using prevalence, cancer penetrance, other death 
#' penetrance, and pedigree matrix. It can then be passed into the peeling-paring 
#' functiton. It can be used for any arbitrary pedigree with the N family members 
#' and does not depend on the number of cancers Y 
#' \itemize{
#'   \item Dim 1 = N family members (counting proband)
#'   \item Dim 2 = 1 + X + (X choose 2) genotypes
#' }
#' Assumes naming conventions for cancers are consistent between arguments
#' BRCAPRO modifies baseline (noncarrier) penetrance by race
#' and penetrance based on whether or not you had a bilateral mastectomy
#' @param fam a family matrix with the following columns: 
#' \itemize{
#'   \item \code{ID} = Member identifier
#'   \item \code{MotherID} = Mother's identifier number
#'   \item \code{FatherID} = Father's identifier number
#'   \item \code{Gender} = 1 for males, 0 for females
#'   \item \code{isProband} = Indicator for whether or not person is proband
#'   \item \code{CurAge} = Family member's current age, or age at death if dead
#'   \item \code{isAffCancer} variables = Indicators for whether or not person 
#'   had cancer type \code{Cancer}, e.g. \code{isAffBC}, \code{isAffO}, etc.
#'   \item \code{AgeCancer} variables = Age of cancer diagnosis for cancer type 
#'   \code{Cancer}, e.g. \code{AgeBC}, \code{AgeO}, etc.
#'   \item \code{isDead} = Indicator for whether or not person is dead
#' }
#' @param CP list of cancer penetrance matrices, separated by each gender
#' @param ODP list of other causes penetrance matrices, separated by each gender
#' @param comprisk list of competing risk matrices, separated by each gender 
#' (the column names must correspond to the possible genotypes, `PG`, in `CP`). 
#' @param max.cancers numeric value indicating the maximum number of cancers an 
#' individual can have. If an individual's number of cancers exceeds max.cancers, 
#' a warning message will be printed. Defaults to 3
#' @family multigene exported
#' @export 
estLik = function(fam, CP, ODP, comprisk=NULL, max.cancers=3){
  N = nrow(fam) # number of people in family
  PG = CP$PG # possible genotypes
  lik = matrix(NA, nrow = N, ncol = length(PG), dimnames=list(1:N, PG))
  
  # Structural checks
  if (length(CP$cancers) != length(ODP$cancers)) {
    stop("Cancer and other causes penetrance matrices don't have same number of cancers")
  }
  
  # Build the matrices of isAff and Age variables for the cancers of interest 
  isAff = subsetCancers("^isAff", fam, ODP$cancers)
  Age = subsetCancers("^Age", fam, ODP$cancers)
  
  # Break up the survival and density penetrance matrices into lists
  cancerDens = list(CP$cancerFDens, CP$cancerMDens)
  cancerSur = list(CP$cancerFSur, CP$cancerMSur)
  deathDens = list(ODP$deathFDens, ODP$deathMDens)
  deathSur = list(ODP$deathFSur, ODP$deathMSur)
  
  for(j in 1:length(PG)){ # Possible genotypes
    for(i in 1:N){ # People
      # Find out how many cancers person j has using isAffCancer indicator 
      idx = which(isAff[i,] != 0) 
      
      # Print a warning message of person j has more than the maximum allowed number of cancers
      if (length(idx) > max.cancers) {
        warning(paste("Person", fam$ID[i], "has more than", max.cancers, "cancers."))
      }
      
      if (length(idx)==0) { # No cancers 
        # Dead
        # For each cancer: multiply 
        # - probability of NOT getting the cancer by person's current age, given genotype;  
        # - with probability of dying from other causes for person's current age
        # Then take the product of all the cancers's probablities
        if (fam$isDead[i]==1) {
          lik[i,j]=prod(cancerSur[[fam$Gender[i]+1]][fam$CurAge[i],j,]
                        *deathDens[[fam$Gender[i]+1]][fam$CurAge[i],])
          # Not dead
          # For each cancer: multiply 
          # - probability of NOT getting the cancer by person's current age, given genotype;  
          # - with probability of surviving other causes for person's current age
          # Then take the product of all the cancers's probablities
        } else if(fam$isDead[i]==0){
          lik[i,j]=prod(cancerSur[[fam$Gender[i]+1]][fam$CurAge[i],j,]
                        *deathSur[[fam$Gender[i]+1]][fam$CurAge[i],])
        } 
      } else if (length(idx)==ncol(Age)) { # Person has all the cancers
        # Multiply 
        # - probability of getting the cancer by the age the person got it, given genotype;  
        # - with probability of surviving other causes for that age
        lik[i,j]=prod(cancerDens[[fam$Gender[i]+1]][cbind(t(Age[i,idx]),j,idx)]
                      *deathSur[[fam$Gender[i]+1]][cbind(t(Age[i,idx]),idx)])
      } else { # At least one cancer, assuming more than one possible cancers
        # Indices for cancers the person doesn't have
        notIdx = which(!(1:ncol(Age) %in% idx)) 
        # Dead
        # For each cancer: multiply 
        # - probability of getting the cancer by the age the person got it, given genotype;  
        # - with probability of surviving other causes for that age
        # For each cancer the person DOESN'T have, multiply 
        # - probability of NOT getting the cancer by person's current age, given genotype;  
        # - with probability of dying from other causes for person's current age
        # Then take the product of all the cancers's probablities
        if (fam$isDead[i]==1){ 
          lik[i,j]=prod(cancerDens[[fam$Gender[i]+1]][cbind(t(Age[i,idx]),j,idx)]
                        *deathSur[[fam$Gender[i]+1]][cbind(t(Age[i,idx]),idx)],
                        cancerSur[[fam$Gender[i]+1]][cbind(fam$CurAge[i],j,notIdx)]
                        *deathDens[[fam$Gender[i]+1]][cbind(fam$CurAge[i],notIdx)])
          # Not dead
          # For each cancer: multiply 
          # - probability of getting the cancer by the age the person got it, given genotype;  
          # - with probability of surviving other causes for that age
          # For each cancer the person DOESN'T have, 
          # - probability of NOT getting the cancer by person's current age, given genotype;  
          # - with probability of surviving other causes for person's current age
          # Then take the product of all the cancers's probablities
        } else if (fam$isDead[i]==0){ 
          lik[i,j]=prod(cancerDens[[fam$Gender[i]+1]][cbind(t(Age[i,idx]),j,idx)]
                        *deathSur[[fam$Gender[i]+1]][cbind(t(Age[i,idx]),idx)],
                        cancerSur[[fam$Gender[i]+1]][cbind(fam$CurAge[i],j,notIdx)]
                        *deathSur[[fam$Gender[i]+1]][cbind(fam$CurAge[i],notIdx)])
        }
      }
    }
  }
  
  # Competing risk
  if (!is.null(comprisk)) {
    # females 
    lik[fam$Gender==0,] = lik[fam$Gender==0,] * comprisk[[1]][fam$CurAge[fam$Gender==0],PG] 
    # males
    lik[fam$Gender==1,] = lik[fam$Gender==1,] * comprisk[[2]][fam$CurAge[fam$Gender==1],PG] 
  }
  
  return(lik)
}
