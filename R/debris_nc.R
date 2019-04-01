
debris_nc <- function(flowframe, p1, p2) {
  #plotting
  flowDensity::plotDens(flowframe, c(p1, p2), 
                        main = flowCore::identifier(flowframe), frame.plot = F)
  #debris gating
  p1_peaks <- flowDensity::getPeaks(flowframe, p1)
  if(length(p1_peaks$Peaks) == 2) {
    #all is well
    deb_cut <- flowDensity::deGate(flowframe, p1, bimodal = T)
    
  } else if(length(p1_peaks$Peaks) > 2) {
    #invader present
    deb_cut <- flowDensity::deGate(flowframe, p1, bimodal = F)#p1_peaks$Peaks[2]
    
  } else if(length(p1_peaks$Peaks) < 2) {
    # not much debris, hence one peak
    deb_cut <- flowDensity::deGate(flowframe, p1, sd.threshold = T, 
                                   n.sd = 1)
    if(deb_cut < 2) {
      deb_cut <- flowDensity::deGate(flowframe, p1, after.peak = T)
    } else if(deb_cut >= 2 & deb_cut >= p1_peaks$Peaks[1]) {
      deb_cut <- flowDensity::deGate(flowframe, p1, use.upper = T, 
                                     upper = F)
    } else {
      deb_cut <- flowDensity::deGate(flowframe, p1, sd.threshold = T, 
                                     n.sd = 1)
    }
  }
  
  abline(v = deb_cut, lty = 2, col = 2)
  
  bs4bs5 <- flowframe[which(flowCore::exprs(flowframe)[, p1] > deb_cut), ]
  other_pos <- which(flowCore::exprs(flowframe)[, p1] > deb_cut)
  deb_pos <- which(flowCore::exprs(flowframe)[, p1] <= deb_cut)
  
  text(x = mean(flowframe@exprs[which(flowCore::exprs(flowframe)[, p1] <= deb_cut), p1]),
       y = mean(flowframe@exprs[which(flowCore::exprs(flowframe)[, p2] <= deb_cut), p2]), 
       "Deb", col = 2)
  
  return(list(bs4bs5 = bs4bs5, deb_pos = deb_pos, bs4bs5_pos = other_pos))
}
