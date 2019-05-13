





bs4_nc <- function(bs4bs5, p1, p2, others, to_retain = "potential") {

  yel.bhlin_peaks <- flowDensity::getPeaks(bs4bs5, p2)
  red.bhlin_peaks <- flowDensity::getPeaks(bs4bs5, p1)
  msg <- capture.output(flowDensity::deGate(bs4bs5, p2, all.cuts = T))[1]

  if(length(yel.bhlin_peaks$Peaks) == 1) {

    yel_cuts <- flowDensity::deGate(bs4bs5, p2, use.percentile = T, percentile = 0.97)
    bs4_pot <- bs4bs5[which(bs4bs5@exprs[, p2] < yel_cuts), ]
    others_pot <- others[which(bs4bs5@exprs[, p2] < yel_cuts)]
    others_nk <- others[which(!(bs4bs5@exprs[, p2] < yel_cuts))]
    ptt <- "1"

  } else if(length(yel.bhlin_peaks$Peaks) == 2) {

    mgroup1 <- yel.bhlin_peaks$Peaks[2] - 0
    if(mgroup1 < 2.9) {
      yel_cuts <- flowDensity::deGate(bs4bs5, p2, use.percentile = T, percentile = 0.99)
      bs4_pot <- bs4bs5[which(bs4bs5@exprs[, p2] < yel_cuts), ]
      others_pot <- others[which(bs4bs5@exprs[, p2] < yel_cuts)]
      others_nk <- others[which(!(bs4bs5@exprs[, p2] < yel_cuts))]
      ptt <- "2a"

    } else {

      yel_cuts <- flowDensity::deGate(bs4bs5, p2, all.cuts = T)
      #BS4 is to the left of the minimum red_cuts
      bs4_pot <- bs4bs5[which(bs4bs5@exprs[, p2] < yel_cuts), ]
      others_pot <- others[which(bs4bs5@exprs[, p2] < yel_cuts)]
      others_nk <- others[which(!(bs4bs5@exprs[, p2] < yel_cuts))]
      ptt <- "2b"

    }
  } else {

    yel_cuts <- flowDensity::deGate(bs4bs5, p2, use.percentile = T, percentile = 0.05)
    #min(yel.bhlin_peaks$Peaks) + 0.50*sd(bs4bs5@exprs[, p2])
    bs4_pot <- bs4bs5[which(bs4bs5@exprs[, p2] > yel_cuts), ]
    others_pot <- others[which(bs4bs5@exprs[, p2] > yel_cuts)]
    others_nk <- others[which(!(bs4bs5@exprs[, p2] > yel_cuts))]
    ptt <- "3"

  }
  #
  abline(h = yel_cuts, lty = 2, col = 2)
  #points(bs4_pot@exprs[, c(p1, p2)], pch = ".", col = 2)
  #text(0, 2.5,
  #     paste("BS4", ptt, sep = ":"),
  #     col = 2)

  if(to_retain == "refined") {

    bs4s <- flowDensity::flowDensity(bs4_pot, channels = c(p1, p2), position = c(F, F),
                                      use.upper = c(T, F),
                                      upper = c(T, NA),
                                      use.percentile = c(F, T),
                                      percentile = c(NA, 0.70),
                                      ellip.gate = T)
    points(bs4s@filter, type = "l", col = 2, lwd = 2, lty = 4)
     text(mean(bs4s@filter[, 1]), mean(bs4s@filter[, 2]), paste("BS4", ptt, sep = "-"),
          col = 2)
     #reduced flowframe for BS4
    bs4_reduced <- bs4_pot[which(!is.na(bs4s@flow.frame@exprs[, 1])), ]
    #positions of BS4s and other unidentified particles
    bs4_pos <- others_pot[which(!is.na(bs4s@flow.frame@exprs[, 1]))]
    others_nk2 <- others_pot[which(is.na(bs4s@flow.frame@exprs[, 1]))]

  } else if(to_retain == "potential") {

    text(mean(bs4_pot@exprs[, p1]), mean(bs4_pot@exprs[, p2]), "BS4",
         col = "red4")
    #reduced flowframe for BS4
    bs4_reduced <- bs4_pot
    #positions of BS4s and other unidentified particles
    others_nk <- others_nk
    bs4_pos <- others_pot
    others_nk2 <- NA

  } else stop("supply to_retain")

  return(list(bs4_reduced = bs4_reduced, others_nk = others_nk, bs4_pos = bs4_pos,
              others_nk2 = others_nk2))

}
