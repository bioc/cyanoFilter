

bs5_ic <- function(bs4bs5, bs4, p1 = "RED.B.HLin", p2 = "YEL.B.HLin", b5_control, bs5_others) {
  msg <- stringr::str_sub(capture.output(flowDensity::flowDensity(bs4bs5,
                                                                  channels = c(p1, p2),
                                                                  position = c(NA, F)))[1], 6, 27)
  
  if(stringr::str_detect(msg, "peak") == T |
     modes::bimodality_coefficient(flowCore::exprs(bs4bs5)[,p2]) < (5/9)) {
    bs5 <- flowDensity::flowDensity(bs4bs5[which(is.na(flowCore::exprs(bs4@flow.frame)[, 1])),],
                                    channels = c(p1, p2), ellip.gate = T,
                                    position = c(T, T), 
                                    use.percentile = c(T, T), percentile = c(0.50, 0.50),
                                    use.control = c(T, T),
                                    control = c(b5_control, b5_control))
  } else {
    bs5 <- flowDensity::flowDensity(bs4bs5[which(is.na(flowCore::exprs(bs4@flow.frame)[, 1])),],
                                    channels= c(p1, p2), ellip.gate = T,
                                    position = c(T, T),
                                    use.percentile = c(F, T), percentile = c(NA, 0.50),
                                    use.control = c(T, T),
                                    control = c(b5_control, b5_control))
  }
  
  #plotting
  points(bs5@filter,type = "l", lty = 4, col = 2)
  text(mean(bs5@filter[,1]), mean(bs5@filter[, 2]), "BS5", col = 2)
  
  #reduced flowframe
  bs5_reduced <- flowDensity::getflowFrame(bs5)
  bs5_pos <- bs5_others[which(!is.na(bs5@flow.frame@exprs[, 1]))]
  others_nk <- bs5_others[which(is.na(bs5@flow.frame@exprs[, 1]))]
  
  return(list(bs5_reduced = bs5_reduced, bs5_pos = bs5_pos, others_nk = others_nk))
}
