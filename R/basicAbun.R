basicAbun <-
function(data, k){
  data <- as.numeric(data)
  
  x <- data[which(data != 0)]
  n <- sum(x)
  D <- length(x)
  n_rare <- sum(x[which(x <= k)])
  D_rare <- length(x[which(x <= k)])
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
  
  BASIC.DATA <- matrix(paste(c("n", "D", "k", "n_rare", "D_rare", "C_rare", "CV_rare", "CV1_rare", "n_abun", "D_abun"),
                             round(c(n, D, k, n_rare, D_rare, C_rare, CV_rare, CV1_rare, n_abun, D_abun),3),
                             sep = "="), ncol = 1)
  colnames(BASIC.DATA) <- c("Value")
  rownames(BASIC.DATA) <- c("Number of observed individuals", "Number of observed species","Cut-off point",
                            "Number of observed in dividuals for rare species", "Number of observed species for rare species",
                            "Estimation of the sample converage for rare species",
                            "Estimation of CV for rare species in ACE", "Estimation of CV1 for rare species in ACE-1",
                            "Number of observed individuals for abundant species", "Number of observed species for abundant species")
  BASIC.DATA <- data.frame(BASIC.DATA)
  return(list(BASIC.DATA, n, D, n_rare, D_rare, C_rare, CV_rare, CV1_rare, n_abun, D_abun))
}
