


bs4_ic <- function(bs4bs5, p1 = "RED.B.HLin", p2 = "YEL.B.HLin", others, b4_control) {
  bs4 <- flowDensity::flowDensity(bs4bs5, channels = c(p1, p2), position = c(NA, F),
                                  use.percentile = c(F, T), ellip.gate = T,
                                  percentile = c(0.975, 0.90),use.control = c(T, T),
                                  control = c(b4_control, b4_control))
  #plotting
  points(bs4@filter, type = "l", lty = 4, col = 2, lwd = 1.5)
  text(mean(bs4@filter[,1]), mean(bs4@filter[,2]),"BS4",col = 2)
  
  bs4_reduced <- flowDensity::getflowFrame(bs4)
  bs4_pos <- others[which(!is.na(bs4@flow.frame@exprs[, 1]))]
  bs5_others <- others[which(is.na(bs4@flow.frame@exprs[, 1]))]
  
  return(list(bs4 = bs4, bs4_reduced = bs4_reduced, bs4_pos = bs4_pos, bs5_others = bs5_others))
}