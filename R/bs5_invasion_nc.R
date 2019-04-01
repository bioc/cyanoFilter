

bs5_inc <- function(bs4bs5, p1, p2, others, take = "liberal", to_retain = "refined") {
  yel.bhlin_peaks <- getPeaks(bs4bs5, p2)
  red.bhlin_peaks <- getPeaks(bs4bs5, p1)
  
  yel_cuts <- deGate(bs4bs5, p2, all.cuts = T)
  red_cuts <- deGate(bs4bs5, p1, all.cuts = T)
  v_cut <- min(red_cuts)
  abline(v = v_cut, col = 2, lty = 2)
  
  if(take == "liberal") {
    bs5_poten <- bs4bs5[which(bs4bs5@exprs[, p2] > min(yel_cuts) & 
                                bs4bs5@exprs[, p1] > min(red_cuts)), ]
    others_bs5_pot <- others[which(bs4bs5@exprs[, p2] > min(yel_cuts) & 
                                     bs4bs5@exprs[, p1] > min(red_cuts))]
    others_nk <- others[which(!(bs4bs5@exprs[, p2] > min(yel_cuts) & 
                      bs4bs5@exprs[, p1] > min(red_cuts)))]
    h_cut <- min(yel_cuts)
    abline(h = h_cut, col = 2, lty = 2)
  } else if(take == "conservative") {
    bs5_poten <- bs4bs5[which(bs4bs5@exprs[, p2] > max(yel_cuts) & 
                                bs4bs5@exprs[, p1] > min(red_cuts)), ]
    others_bs5_pot <- others[which(bs4bs5@exprs[, p2] > max(yel_cuts) & 
                                     bs4bs5@exprs[, p1] > min(red_cuts))]
    others_nk <- others[which(!(bs4bs5@exprs[, p2] > max(yel_cuts) & 
                      bs4bs5@exprs[, p1] > min(red_cuts)))]
    h_cut <- max(yel_cuts)
    abline(h = h_cut, col = 2, lty = 2)
  } else stop("Take a stand! conservative or liberal? With apologies to the independents")
  
  
  if(to_retain == "refined") {
    bs5s <- flowDensity::flowDensity(bs5_poten, channels = c(p1, p2),
                                     position = c(F, T), 
                                     ellip.gate = T,
                                     n.sd = c(1.5, 2), 
                                     use.percentile = c(T, F),
                                     percentile = c(0.95, NA),
                                     use.upper = c(T, T), 
                                     upper = c(T, F))
    points(bs5s@filter, type = "l", col = "red", lwd = 2, lty = 2)
    text(mean(bs5s@filter[,1]), mean(bs5s@filter[, 2]), "BS5", col = 2)
    
    bs5_reduced <- bs5_poten[which(!is.na(bs5s@flow.frame@exprs[, 1])), ]
    bs5_pos <- others_bs5_pot[which(!is.na(bs5s@flow.frame@exprs[, 1]))]
    others_nk2 <- others_bs5_pot[which(is.na(bs5s@flow.frame@exprs[, 1]))]
  } else if(to_retain == "potential") {
    text(mean(bs5_poten@exprs[, p1]), mean(bs5_poten@exprs[, p2]), "BS5", col = 2)
    bs5_reduced <- bs5_poten
    bs5_pos <- others_bs5_pot
    others_nk2 <- NA
  } else stop("supply to_retian")
  
  return(list(bs5_reduced = bs5_reduced, others_bs5 = others_bs5, others_nk = others_nk, 
              others_nk2 = others_nk2, v_cut = v_cut, h_cut = h_cut))
}