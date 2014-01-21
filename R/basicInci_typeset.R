basicIncitype <- function(data, k){
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
    
    basicInfo <- matrix(c(D, t, k), nrow = 1); colnames(basicInfo) <- c("S.obs", "t", "cutpt")
      rownames(basicInfo) <- "BasicInfo"
    infreqInfo <- matrix(D_infreq, nrow = 1); colnames(infreqInfo) <- c("S.obs_infreq")
      rownames(infreqInfo) <- "InfreqInfo"
    infreqInfo2 <- matrix(c(C_infreq, CV_infreq), nrow = 1); colnames(infreqInfo2) <- c("C_infreq", "CV_infeq")
      rownames(infreqInfo2) <- "OtherInfo"
    freqInfo <- matrix(c(D_freq), nrow = 1); colnames(freqInfo) <- c("S.obs_freq")
      rownames(freqInfo) <- "FreqInfo"
    
    return(list(cbind(basicInfo, infreqInfo, freqInfo), infreqInfo2))
  }
