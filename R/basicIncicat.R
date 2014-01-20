basicIncicat <- function(data, k){
  data <- as.numeric(data)
  t <- data[1]
  dat <- data[-1]
  x <- dat[which(dat != 0)]
  Q <- function(i, data){length(data[which(data == i)])}
  
  D <- length(x)
  D_infreq <- length(x[which(x <= k)])
  
  if (Q(1, x) > 0 & Q(2, x) > 0){
    A <- 2*Q(2, x)/((t-1)*Q(1, x) + 2*Q(2, x))
  } else if (Q(1, x) > 0 & Q(2, x) == 0){
    A <- 2/((t-1)*(Q(1, x) - 1) + 2)
  } else {
    A <- 1
  }
  C_infreq <- 1 - Q(1, x)/sum(x[which(x <= k)])*(1-A)
  
  j <- c(1:k)
  b1 <- sum(sapply(j, function(j)j*(j-1)*Q(j, x)))
  b2 <- sum(sapply(j, function(j)j*Q(j, x)))
  gamma_infreq_square <- max(D_infreq/C_infreq*t/(t - 1)*b1/b2/(b2) - 1, 0)
  CV_infreq <- sqrt(gamma_infreq_square)
  D_freq <- length(x[which(x > k)])
    
  Q <- function(i, data){length(data[which(data == i)])}
  rsg <- matrix(sapply(j, function(j)f(j, x)), 1, k)
  
  #cat("(1) BASIC DATA INFORMATION  \n ", fill = TRUE)
  cat("             (Number of observed species)  D ", sep = " = ", D, fill = TRUE)    
  cat("             (Number of samples/quadrats)  t ", sep = " = ", t, fill = TRUE)    
  cat("             (Cut-off point)               k ", sep = " = ", k, fill = TRUE)    
  cat("\n", fill = TRUE)    
  cat("   \"Infrequent\" Shared Species Group: (Incidence counts up to cut-off point k)", fill = TRUE)
  cat("          Some Statistics:", fill = TRUE)
  cat("          ---------------------------------------------------------------------", fill = TRUE)
  cat("         ", paste("Q1 =", rsg[1], "; ","Q2 =", rsg[2], "; ",
                         "Q3 =", rsg[3], "; ","Q4 =", rsg[4], "; ",
                         "Q5 =", rsg[5], "; ","Q6 =", rsg[6], "; ",
                         "Q7 =", rsg[7], "; ","Q8 =", rsg[8], "; ",
                         "Q9 =", rsg[9], "; ","Q10 =", rsg[10]), fill = TRUE)
  cat("          ---------------------------------------------------------------------", fill = TRUE)
  cat("               (Number of observed species for infrequent species)  D_infreq " , sep = " = ", D_infreq, fill = TRUE)
  cat("               (Estimated sample converage for infrequent species)  C_infreq " , sep = " = ", C_infreq, fill = TRUE)
  cat("               (Estimated CV for infrequent species)                CV_infreq" , sep = " = ", CV_infreq, fill = TRUE)
  cat("\n", fill = TRUE)    
  cat("   \"Frequent\" Species Group: (Incidence counts over the cut-off point) ", fill = TRUE)
  cat("               (Number of observed species for frequent species)       D_freq",  sep = " = ", D_freq, fill = TRUE)
  cat("\n", fill = TRUE)    
}