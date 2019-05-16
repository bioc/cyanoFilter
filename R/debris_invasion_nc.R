#' gates out or assign indicators to BS4 cyano cells from a flowframe.
#'
#' @param bs4bs5 flowframe with.
#' @param p1 first flowcytometer channel that can be used to separate BS4 cells from the rest, e.g. "RED.B.HLin".
#' @param p2 second flowcytometer channel that can be used to separate BS4 cells from the rest, e.g. "YEL.B.HLin"
#' @param others row numbers for non-debris events. This is provided by the debris_nc or debris_inc function.
#' @param  retain should potential candidates be retained or further gating be applied to filter out only certain BS4 cells.
#' @return list containing; \itemize{
#' \item \strong{bs4_reduced -}
#' \item \strong{others_nk -}
#' \item \strong{bs4_pos -}
#' \item \strong{others_nk2 -}
#' }
#'
#' @examples
#' \dontrun{
#' celldebris_nc(bs4bs5 = flowfile, p1 = "RED.B.HLin", p2 = "YEL.B.HLin", others = b4b5_others, to_retain = "refined")
#' }
#'
#'
#' @export celldebris_nc


debris_inc <- function(flowframe, p1, p2) {
    # plotting
    flowDensity::plotDens(flowframe, c(p1, p2), main = flowCore::identifier(flowframe), frame.plot = F)
    # debris gating
    p1_peaks <- flowDensity::getPeaks(flowframe, p1)
    if (length(p1_peaks$Peaks) == 2) {
        # all is well
        deb_cut <- flowDensity::deGate(flowframe, p1, bimodal = T)
        if (p1_peaks$Peaks[1] - (1 * sd(flowframe@exprs[, p1])) > 0 & deb_cut >= p1_peaks$Peaks[1]) {
            deb_cut <- flowDensity::deGate(flowframe, p1, use.upper = T, upper = F)
            ptt <- "1a"
        } else {

            deb_cut <- deb_cut
            ptt <- "1b"

        }

    } else if (length(p1_peaks$Peaks) == 3) {
        # invader present
        deb_cut <- flowDensity::deGate(flowframe, p1, all.cuts = T)
        if (p1_peaks$Peaks[1] - (1 * sd(flowframe@exprs[, p1])) <= 0) {

            deb_cut <- flowDensity::deGate(flowframe, p1, all.cuts = T)[1]
            ptt <- "2a"

        } else if (p1_peaks$Peaks[1] - (1 * sd(flowframe@exprs[, p1])) > 0) {

            deb_cut <- flowDensity::deGate(flowframe, p1, all.cuts = T)[1]

            flowframe2 <- flowframe[which(flowCore::exprs(flowframe)[, p1] > deb_cut), ]
            deb_cut <- ifelse(length(flowDensity::getPeaks(flowframe2, p1, tinypeak.removal = 1/5)$Peaks) < 2, flowDensity::deGate(flowframe, p1, use.upper = T,
                upper = F), deb_cut)
            ptt <- "2b"

        } else {

            deb_cut <- flowDensity::deGate(flowframe, p1, bimodal = F)
            ptt <- "2c"

        }


    } else if (length(p1_peaks$Peaks) < 2) {
        # not much debris, hence one peak
        deb_cut <- flowDensity::deGate(flowframe, p1, sd.threshold = T, n.sd = 1)
        if (deb_cut < 2) {

            deb_cut <- flowDensity::deGate(flowframe, p1, after.peak = T)
            ptt <- "3a"

        } else if (p1_peaks$Peaks[1] - (2 * sd(flowframe@exprs[, p1])) > 0 & deb_cut >= p1_peaks$Peaks[1]) {

            deb_cut <- flowDensity::deGate(flowframe, p1, use.upper = T, upper = F)
            ptt <- "3b"

        } else {

            deb_cut <- deb_cut
            ptt <- "3c"
        }
    } else {

        deb_cut <- flowDensity::deGate(flowframe, p1, all.cuts = T)[2]  #p1_peaks$Peaks[2] - (0.5 * sd(flowframe@exprs[, p1]))
        ptt <- "4"
    }

    abline(v = deb_cut, lty = 2, col = 2)
    # abline(v = flowDensity::deGate(flowframe, p1, all.cuts = T), lty = 2, col = 1)

    bs4bs5 <- flowframe[which(flowCore::exprs(flowframe)[, p1] > deb_cut), ]
    other_pos <- which(flowCore::exprs(flowframe)[, p1] > deb_cut)
    deb_pos <- which(flowCore::exprs(flowframe)[, p1] <= deb_cut)

    text(x = mean(flowframe@exprs[which(flowCore::exprs(flowframe)[, p1] <= deb_cut), p1]), y = mean(flowframe@exprs[which(flowCore::exprs(flowframe)[, p2] <=
        deb_cut), p2]), paste("Deb", ptt, sep = "-"), col = 2)

    return(list(bs4bs5 = bs4bs5, deb_pos = deb_pos, bs4bs5_pos = other_pos))
}
