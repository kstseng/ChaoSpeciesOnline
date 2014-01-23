InfreqSpeciesGroupImprove <- function(data, k){
  if (is.matrix(data) == T || is.data.frame(data) == T){
    if (ncol(data) != 1 & nrow(data) != 1)
      stop("Error: The data format is wrong.")
    if (ncol(data) == 1){
      data <- data[, 1]
    } else {
      data <- data[1, ]
    }
  }
  data <- as.numeric(data)
  
  data <- data[-1]
  Q <- function(i, data){length(data[which(data == i)])}
  
  x <- data[which(data != 0)]
  r <- c(1:k)
  Rare.Species.Group <- matrix(sapply(r, function(r)Q(r, x)), 1, k)
  rownames(Rare.Species.Group) <- ""
  colnames(Rare.Species.Group) <- paste("Q", r, sep="")
  return(Rare.Species.Group)
}
