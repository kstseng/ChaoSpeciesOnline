SpecAbunOut <-
  function(data, method = c("all", "Homogeneous", "Chao", "CE", "Jackknife"),
           k, conf){
    data <- as.numeric(data)
    #method <- match.arg(method, several.ok=T)
    ind <- pmatch(method, c("Homogeneous", "Chao", "CE", "Jackknife"))
    Homo <- round(rbind(SpecAbunHomo(data, k, conf), SpecAbunHomoMle(data, k, conf)),3)
    Chao <- round(rbind(SpecAbunChao1(data, k, conf), SpecAbunChao1bc(data, k, conf), 
                  SpecAbuniChao1(data, k, conf)),3)
    CE <- round(rbind(SpecAbunAce(data, k, conf), SpecAbunAce1(data, k, conf)),3)
    Jackknife <- round(rbind(SpecAbunJack1(data, k, conf), SpecAbunJack2(data, k, conf)),3)
    tmp1 <- rbind(Homo, Chao, CE, Jackknife)
    tmp2 <- list(Homo, Chao, CE, Jackknife)
    
    if (sum(method == "all") != 0){ #if "all" in method
      out <- tmp1
    }else{
      out <- do.call("rbind", lapply(ind, function(ind)tmp2[[ind]]))
    }
    
    return(out)
  }
