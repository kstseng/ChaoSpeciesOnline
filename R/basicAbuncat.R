basicAbuncat <- function(data, k){
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
    
    r <- c(1:k)
    rsg <- matrix(sapply(r, function(r)f(r, x)), 1, k)
    
    cat("(1) BASIC DATA INFORMATION  \n ", fill = TRUE)
    cat("             (Number of observed individuals)  n ", sep = " = ", n, fill = TRUE)  
    cat("             (Number of observed species)      D ", sep = " = ", D, fill = TRUE)    
    cat("             (Cut-off point)                   k ", sep = "  = ", k, fill = TRUE)    
    cat("\n", fill = TRUE)    
    cat("   \"Rare\" Shared Species Group: (Frequencies counts up to cut-off point k)", fill = TRUE)
    cat("          Some Statistics:", fill = TRUE)
    cat("          ---------------------------------------------------------------------", fill = TRUE)
    cat("         ", paste("f1 =", rsg[1], "; ","f2 =", rsg[2], "; ",
                           "f3 =", rsg[3], "; ","f4 =", rsg[4], "; ",
                           "f5 =", rsg[5], "; ","f6 =", rsg[6], "; ",
                           "f7 =", rsg[7], "; ","f8 =", rsg[8], "; ",
                           "f9 =", rsg[9], "; ","f10 =", rsg[10]), fill = TRUE)
    cat("          ---------------------------------------------------------------------", fill = TRUE)
    cat("               (Number of observed individuals for rare species)      n_rare",  sep = " = ", n_rare, fill = TRUE)
    cat("               (Number of observed species for rare species)          D_rare",  sep = " = ", D_rare, fill = TRUE)
    cat("               (Estimation of the sample converage for rare species)  C_rare" , sep = " = ", C_rare, fill = TRUE)
    cat("               (Estimation of CV for rare species in ACE)             CV_rare", sep = " = ", CV_rare, fill = TRUE)
    cat("               (Estimation of CV1 for rare species in ACE-1)          CV1_rare", sep = " = ", CV1_rare, fill = TRUE)  
    cat("\n", fill = TRUE)    
    cat("   \"Abundant\" Species Group: (Frequencies beyond the cut-off point) ", fill = TRUE)
    cat("               (Number of observed individuals for abundant species)      n_abun",  sep = " = ", n_abun, fill = TRUE)
    cat("               (Number of observed species for abundant species)          D_abun",  sep = " = ", D_abun, fill = TRUE)
    cat("\n", fill = TRUE)    
  }
