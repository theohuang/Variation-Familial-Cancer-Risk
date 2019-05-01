#' Generate a List of Other Death Penetrance Matrices
#' 
#' Returns a list of penetrance matrices (survival functions and density 
#' functions) for non-cancer related deaths, separated by sex. The output 
#' list also includes the input vector of cancer names.
#' @param cancers names of Y cancers to use
#' @param deathOtherCauses data frame of death by other causes
#' @param age.max maximum age to consider
#' @details Assumes that the naming conventions for the different cancers in 
#' \code{cancers} and \code{deathOtherCauses} are consistent, and that 
#' \code{deathOtherCauses} distinguishes between sexes by prefixing the column 
#' names by "female" and "male. 
#' @family multigene exported
#' @export 
genOtherDeathPen = function(cancers, deathOtherCauses, age.max = 110){
  # Survival matrices (1 column corresponding to each cancer)
  deathFSur = sapply(cancers, function(canc){
    exp(-cumsum(deathOtherCauses[,paste0("female",canc)])) 
  })
  deathMSur = sapply(cancers, function(canc){
    exp(-cumsum(deathOtherCauses[,paste0("male",canc)])) 
  })
  # Density matrices (1 column corresponding to each cancer)
  deathFDens = sapply(cancers, function(canc){
    deathFSur[,canc]*deathOtherCauses[,paste0("female",canc)]
  })
  deathFDens[age.max,] = 1-safeApply(deathFDens[1:(age.max - 1),], 2, sum)
  deathMDens = sapply(cancers, function(canc){
    deathMSur[,canc]*deathOtherCauses[,paste0("male",canc)]
  })
  deathMDens[age.max,] = 1-safeApply(deathMDens[1:(age.max - 1),], 2, sum)
  
  # Output as a list
  ODP = list() 
  ODP$cancers = cancers # Names of cancers, supplied as input
  ODP$deathFSur = deathFSur # Female survival matrix
  ODP$deathMSur = deathMSur # Male survival matrix
  ODP$deathFDens = deathFDens # Female density matrix
  ODP$deathMDens = deathMDens # Male density matrix
  
  return(ODP)
}