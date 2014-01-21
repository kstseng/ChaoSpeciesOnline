basicIncitype <- function(data, k){
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
    rownames(basicInfo) <- ""
  infreqInfo <- matrix(D_infreq, nrow = 1); colnames(infreqInfo) <- c("S.obs_infreq")
    rownames(infreqInfo) <- "InfreqInfo"
  infreqInfo2 <- matrix(c(C_infreq, CV_infreq), nrow = 1); colnames(infreqInfo2) <- c("C_infreq", "CV_infeq")
    rownames(infreqInfo2) <- ""
  freqInfo <- matrix(c(D_freq), nrow = 1); colnames(freqInfo) <- c("S.obs_freq")
    rownames(freqInfo) <- "FreqInfo"
  
  InfreqSpeciesGroup <- function(data, k){
    data <- data[-1]
    data <- as.numeric(data)
    Q <- function(i, data){length(data[which(data == i)])}
    
    x <- data[which(data != 0)]
    r <- c(1:k)
    Rare.Species.Group <- matrix(sapply(r, function(r)Q(r, x)), 1, k)
    rownames(Rare.Species.Group) <- ""
    colnames(Rare.Species.Group) <- paste("Q", r, sep="")
    return(Rare.Species.Group)
  }
  
  
  return(list(Result1 = cbind(basicInfo, infreqInfo, freqInfo),
              Result2 = infreqInfo2, 
              Result3 = InfreqSpeciesGroup(data, k)))
  }
