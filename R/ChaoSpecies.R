ChaoSpeciesOnline <-
function(data, datatype = c("abundance", "incidence"), method = c("all", "Homogeneous", "Chao", "CE", "Jackknife"),
         k = 10, conf = 0.95, detail = TRUE){
  if (sum(is.null(method)) != 0)
    stop("Error: You haven't chosen the estimated method!")
  #method <- match.arg(method, several.ok=T)
  if (k != round(k) || k < 0) 
    stop("Error: The cutoff t to define less abundant species must be non-negative integer!")
  if (is.numeric(conf) == FALSE || conf > 1 || conf < 0) 
    stop("Error: confidence level must be a numerical value between 0 and 1, e.g. 0.95")
  
  if (is.matrix(data) == T || is.data.frame(data) == T){
    if (ncol(data) != 1 & nrow(data) != 1)
      stop("Error: The data format is wrong.")
    if (ncol(data) == 1){
      data <- data[, 1]
    } else {
      data <- data[1, ]
    }
  }
  
  data <- as.numeric(round(data))
  
  if (datatype == "abundance"){
      f <- function(i, data){length(data[which(data == i)])}
      if (f(1, data) == sum(data)){
        stop("Error: The information of data is not enough.")}
      if (detail == T) basicAbuncat(data, k)
      SpeciesAbundance(data, method = method, k = k, conf = conf)
    } else {
      dat <- data[-1]; Q <- function(i, data){length(data[which(data == i)])}
      if (Q(1, dat) == sum(dat)){
        stop("Error: The information of data is not enough.")}
      if (detail == T) basicIncicat(data, k)
      SpeciesIncidence(data, method = method, k = k, conf = conf)  
    }
}
