basicAbuntype <- function(data, k){
  data <- as.numeric(data)
  
  x <- data[which(data != 0)]
  n <- sum(x)
  D <- length(x)
  n_rare <- sum(x[which(x <= k)])
  D_rare <- length(x[which(x <= k)])
  f <- function(i, data){length(data[which(data == i)])}
  if (n_rare != 0){
    C_rare <- 1 - f(1, x)/n_rare
  } else {
    C_rare = 1
  } 
  n_abun <- n - n_rare
  D_abun <- length(x[which(x > k)])
  
  j <- c(1:k)
  a1 <- sum(sapply(j, function(j)j*(j - 1)*f(j, x)))
  a2 <- sum(sapply(j, function(j)j*f(j, x)))
  if (C_rare != 0){
    gamma_rare_hat_square <- max(D_rare/C_rare*a1/a2/(a2 - 1) - 1, 0)
    gamma_rare_1_square <- max(gamma_rare_hat_square*(1 + (1 - C_rare)/C_rare*a1/(a2 - 1)), 0)
  }else{
    gamma_rare_hat_square <- 0
    gamma_rare_1_square <- 0
  }
  CV_rare <- sqrt(gamma_rare_hat_square)
  CV1_rare <- sqrt(gamma_rare_1_square)

  basicInfo <- matrix(c(n, D, k), nrow = 1); colnames(basicInfo) <- c("n", "S.obs", "cutpt")
    rownames(basicInfo) <- "BasicInfo"
  rareInfo <- matrix(c(n_rare, D_rare), nrow = 1); colnames(rareInfo) <- c("n_rare", "S.obs_rare")
    rownames(rareInfo) <- "RareInfo"
  rareInfo2 <- matrix(c(C_rare, CV_rare, CV1_rare), nrow = 1); colnames(rareInfo2) <- c("C_rare", "CV_rare", "CV1_rare")
    rownames(rareInfo2) <- "OtherInfo"
  abunInfo <- matrix(c(n_abun, D_abun), nrow = 1); colnames(abunInfo) <- c("n_abun", "S.obs_abun")
    rownames(abunInfo) <- "AbunInfo"
  
  return(list(cbind(basicInfo, rareInfo, abunInfo), rareInfo2))
  
}
