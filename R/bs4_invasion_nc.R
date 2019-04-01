



bs4_inc <- function(bs4bs5, p1, p2, others, take = "conservative", to_retain = "refined") {
  
  yel_cuts <- deGate(bs4bs5, p2, after.peak = T)
  
  if(take == "liberal") {
    bs4_pot <- bs4bs5[which(bs4bs5@exprs[, p2] < max(yel_cuts)), ]
    bs4_pot_pos <- others[which(bs4bs5@exprs[, p2] < max(yel_cuts))]
    others_nk <- others[which(!(bs4bs5@exprs[, p2] < max(yel_cuts)))]
    abline(h = max(yel_cuts), col = 2, lty = 2)
    h_cut <- max(yel_cuts)
  } else if(take == "conservative") {
    bs4_pot <- bs4bs5[which(bs4bs5@exprs[, p2] < min(yel_cuts)), ]
    bs4_pot_pos <- others[which(bs4bs5@exprs[, p2] < min(yel_cuts))]
    others_nk <- others[which(!(bs4bs5@exprs[, p2] < min(yel_cuts)))]
    abline(h = min(yel_cuts), col = 2, lty = 2)
    h_cut <- min(yel_cuts)
  } else stop("Take a stand! conservative or liberal? With apologies to the independents")
  
  red.bhlin_peaks <- getPeaks(bs4_pot, p1)
  p4_peaks <- flowDensity::getPeaks(bs4bs5, p1)
  abline(v = red.bhlin_peaks$Peaks, col = 5, lty = 1)
  red_cuts <- deGate(bs4bs5, p1, all.cuts = T)
  v_cut <- red_cuts
  abline(v = red_cuts, lty = 2)
  
  if(length(p4_peaks$Peaks) == 1) {
    #BS4_pot2 is all cells below the yel line
    bs4_pot2 <- bs4_pot
    b4pot <- bs4_pot_pos
  } else if(length(p4_peaks$Peaks) > 1 ) {
    bs4_pot2 <- bs4_pot
    b4pot <- bs4_pot_pos
  } else {
    abline(v = red_cuts, lty = 2)
  }
  
  text(mean(bs4_pot@exprs[, p1]), mean(bs4_pot@exprs[, p2]), length(p4_peaks$Peaks), col = 2, 
       cex = 3)
  #head(b4pot)
  points(bs4_pot2@exprs[, p1], bs4_pot2@exprs[, p2], col = 3, pch = ".")
   if(to_retain == "refined") {
     
   } else if(to_retain == "potential") {
     
   } else stop("supply to_retian")
  
  
  
  #function to find the closest point to the computed peaks
  # p4_peaks <- flowDensity::getPeaks(bs4bs5, p1)
  # .Peaks_Ind <- function(Peaks = p4_peaks$Peaks, bs4bs5, p1) {
  #     cpoints <- unlist(sapply(1:length(Peaks), function(i) {
  #       #find the closest point from the left
  #       dp_max <- max((bs4bs5@exprs[, p1] - Peaks[i])[(bs4bs5@exprs[, p1] - Peaks[i]) < 0])
  #       cpoint <- which((bs4bs5@exprs[, p1] - Peaks[i]) == dp_max)
  #     return(cpoint)
  #   }))
  # 
  #   return(cpoints)
  # }
  # 
  # 
  # if(length(p4_peaks$Peaks) > 2) {
  #   
  #   
  #   bs4_1 <- flowDensity::flowDensity(bs4bs5, channels = c(p1, p2), position = c(T, F),
  #                                     ellip.gate = T, use.upper = c(T, F), 
  #                                     upper = c(F, NA))
  #   ttext <- 1
  #   
  # } else if(length(p4_peaks$Peaks) == 1) {
  #   
  #   bs4_1 <- flowDensity::flowDensity(bs4bs5, channels = c(p1, p2), position = c(F, F),
  #                                     ellip.gate = T, use.percentile = c(T, F), 
  #                                     percentile = c(0.95, NA))
  #   ttext <- 2
  # } else if(length(p4_peaks$Peaks) == 2) {
  #   P.Ind <- .Peaks_Ind(Peaks = p4_peaks$Peaks, bs4bs5 = bs4bs5, p1 = p1)
  #   
  #   peak_diff <- max(bs4bs5@exprs[P.Ind, p2]) - min(bs4bs5@exprs[P.Ind, p2])
  #   if(floor(peak_diff) > sd(bs4bs5@exprs[, p2])) {
  #     bs4_1 <- flowDensity::flowDensity(bs4bs5, channels = c(p1, p2), position = c(F, F),
  #                                       ellip.gate = T)
  #   } else {
  #     bs4_1 <- flowDensity::flowDensity(bs4bs5, channels = c(p1, p2), position = c(T, F),
  #                                       ellip.gate = T)
  #   }
  #   
  #   ttext <- ifelse(length(peak_diff) > 1, 3.5, 3)
  #   
  # } else {
  #   
  #   bs4_1 <- flowDensity::flowDensity(bs4bs5, channels = c(p1, p2), position = c(F, NA),
  #                                     ellip.gate = T) 
  #   ttext <- 4
  # }
  #   #plotting
  # points(bs4_1@filter, type = "l", col = 2, lwd = 2, lty = 4)
  # text(mean(bs4_1@filter[, 1]), mean(bs4_1@filter[, 2]), paste0("BS4", ",", ttext), col = 2)
  # 
  # BS5_and_others <- others[which(is.na(bs4_1@flow.frame@exprs[, 1]))]
  # bs4_pos <- others[which(!is.na(bs4_1@flow.frame@exprs[, 1]))]
  # 
  # #reduced flowframe for BS4
  # bs4_reduced <- bs4bs5[which(!is.na(bs4_1@flow.frame@exprs[, 1])), ]
  # bs5_rem <- bs4bs5[which(is.na(bs4_1@flow.frame@exprs[, 1])), ]
  return(list(bs4_reduced = NA, others_bs4 = NA, others_nk = NA, 
              others_nk2 = NA, v_cut = v_cut, h_cut = h_cut))
}