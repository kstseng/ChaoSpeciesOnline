SpecAbunOut <-
function(data, 
                       method = c("all", "Homogeneous", "Chao", "CE", "Jackknife"),
                       k, conf){
  data <- as.numeric(data)
  method <- match.arg(method)
  if (method == "all") {
    a <- SpecAbunHomo(data, k, conf)
    b <- SpecAbunHomoMle(data, k, conf)
    c <- SpecAbunChao1(data, k, conf)
    d <- SpecAbunChao1bc(data, k, conf)
    e <- SpecAbuniChao1(data, k, conf)
    f <- SpecAbunAce(data, k, conf)
    g <- SpecAbunAce1(data, k, conf)
    h <- SpecAbunJack1(data, k, conf)
    i <- SpecAbunJack2(data, k, conf) 
    out <- rbind(a, b, c, d, e, f, g, h)
  }
  
  if (method == "Homogeneous")
    out <- rbind(SpecAbunHomo(data, k, conf), SpecAbunHomoMle(data, k, conf))
  if (method == "Chao")
    out <- rbind(SpecAbunChao1(data, k, conf), SpecAbunChao1bc(data, k, conf), SpecAbuniChao1(data, k, conf))
  if (method == "CE")
    out <- rbind(SpecAbunAce(data, k, conf), SpecAbunAce1(data, k, conf))
  if (method == "Jackknife")
    out <- rbind(SpecAbunJack1(data, k, conf), SpecAbunJack2(data, k, conf))
  
  return(out)
}
