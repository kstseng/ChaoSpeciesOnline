SpecInciOut <-
function(data, method = c("all", "Homogeneous", "Chao", "CE", "Jackknife"), k, conf){
    ind <- pmatch(method, c("Homogeneous", "Chao", "CE", "Jackknife"))
    Homo <- round(SpecInciHomo(data, k, conf), 3)
    Chao <- rbind(round(SpecInciChao2(data, k, conf), 3), round(SpecInciChao2bc(data, k, conf), 3),
                  round(SpecInciiChao2(data, k, conf), 3))
    CE <- rbind(round(SpecInciModelh(data, k, conf)[[1]], 3), round(SpecInciModelh1(data, k, conf)[[1]], 3))
    Jackknife <- rbind(round(SpecInciJack1(data, k, conf), 3), round(SpecInciJack2(data, k, conf), 3))

    est.cv <- data.frame(c("", "", "", "", round(SpecInciModelh(data, k, conf)[[2]], 3), 
                           round(SpecInciModelh1(data, k, conf)[[2]], 3), "", ""))
    colnames(est.cv) <- "Est.CV(rare)"
    
#     est.cv_J <- matrix(c(round(SpecInciModelh(data, k, conf)[[2]], 3), round(SpecInciModelh1(data, k, conf)[[2]], 3)),
#                        nrow = 2)
    est.cv_J <- c(round(SpecInciModelh(data, k, conf)[[2]], 3), round(SpecInciModelh1(data, k, conf)[[2]], 3))
                      
    Jackknife.cv <- cbind(Jackknife, est.cv_J)
    tmp1 <- cbind(rbind(Homo, Chao, CE, Jackknife), est.cv)
    tmp2 <- list(Homo, Chao, CE, Jackknife)
    
    if (sum(method == "all") != 0){ #if "all" in method
      out <- tmp1
    }else if (sum(method == "Jackknife") != 0){
      out <- do.call("rbind", lapply(ind, function(ind)tmp2[[ind]]))
      det <- which(rownames(out) == "1st order jackknife" | rownames(out) == "2nd order jackknife")
      value.na <- matrix(rep(NA, nrow(out)), ncol = 1)
      colnames(value.na) <- "Est.CV(rare)"
      
      out <- cbind(out, value.na)   
      out[det, 5] = est.cv_J
      
    }else{
      out <- do.call("rbind", lapply(ind, function(ind)tmp2[[ind]]))
    }
  
  return(out)
}
