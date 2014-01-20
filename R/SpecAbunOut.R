SpecAbunOut <-
function(data, method = c("all", "Homogeneous", "Chao", "CE", "Jackknife"),
         k, conf){
  data <- as.numeric(data)
  method <- match.arg(method)
  a <- SpecAbunHomo(data, k, conf)
  b <- SpecAbunHomoMle(data, k, conf)
  c <- SpecAbunChao1(data, k, conf)
  d <- SpecAbunChao1bc(data, k, conf)
  e <- SpecAbuniChao1(data, k, conf)
  f <- SpecAbunAce(data, k, conf)
  g <- SpecAbunAce1(data, k, conf)
  h <- SpecAbunJack1(data, k, conf)
  i <- SpecAbunJack2(data, k, conf) 
  if (method == "all") {
    out <- rbind(a, b, c, d, e, f, g, h)
  }
  
  if (method == "Homogeneous")
    out <- rbind(a, b)
  if (method == "Chao")
    out <- rbind(c, d, e)
  if (method == "CE")
    out <- rbind(f, g)
  if (method == "Jackknife")
    out <- rbind(h, i)
  
  
  
  return(out)
}
