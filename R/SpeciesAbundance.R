SpeciesAbundance <-
function(data, method = c("all", "Homogeneous", "Chao", "CE", "Jackknife"), k = 10, conf = 0.95){
  #method <- match.arg(method)
  return(SPECIES.TABLE = round(SpecAbunOut(data, method, k, conf), 3))
}
