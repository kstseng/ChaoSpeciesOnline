SpecInciOut <-
function(data, method = c("all", "Homogeneous", "Chao", "CE", "Jackknife"), k, conf){
    a <- round(SpecInciHomo(data, k, conf), 3)
    b <- round(SpecInciChao2(data, k, conf), 3)
    c <- round(SpecInciChao2bc(data, k, conf), 3)
    d <- round(SpecInciiChao2(data, k, conf), 3)
    e <- round(SpecInciModelh(data, k, conf)[[1]], 3)
    f <- round(SpecInciModelh1(data, k, conf)[[1]], 3)
    g <- round(SpecInciJack1(data, k, conf), 3)
    h <- round(SpecInciJack2(data, k, conf), 3)
    est.cv <- data.frame(c("", "", "", "", round(SpecInciModelh(data, k, conf)[[2]], 3), 
                           round(SpecInciModelh1(data, k, conf)[[2]], 3), "", ""))
    colnames(est.cv) <- "Est.CV(rare)"
  if (method == "all") {
    out <- cbind(rbind(a, b, c, d, e, f, g, h), est.cv)
  }
  
  if (method == "Homogeneous")
    out <- a
  if (method == "Chao")
    out <- rbind(b, c, d)
  if (method == "CE"){
    est.cv <- matrix(c(SpecInciModelh(data, k, conf)[[2]], SpecInciModelh1(data, k, conf)[[2]]), ncol = 1)
    est.cv <- round(est.cv, 3)
    colnames(est.cv) <- "Est.CV(rare)"
    out1 <- rbind(e, f)
    out <- cbind(out1, est.cv)
  }
  if (method == "Jackknife")
    out <- rbind(g, h)
  
  return(out)
}
