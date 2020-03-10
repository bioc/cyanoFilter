#' gates out or assign indicators to phytoplankton cells based on the expression of
#' Phycocyanin.
#'
#' @param bs4bs5 flowframe with debris (left) removed.
#'
#' @param ch1 first flowcytometer channel measuring first pigmentation.
#'            For example, channel measuring phycocyanin
#'
#' @param ch2 second flowcytometer channel measuring second pigmentation.
#'            For example, channel measuring phycoerytherin
#'
#' @param others row numbers for non-debris events. This is provided by the
#'               \code{\link{debris_nc}} or \code{\link{debris_inc}} function.
#' @param to_retain should potential candidates be retained or further gating be applied
#'                  to filter out only certain phytoplankton cells.
#' @return list containing; \itemize{
#' \item \strong{syn_reduced -} flowframe containing only phytoplankton cells
#' \item \strong{others_nk -} unidentified particle positions
#' \item \strong{syn_pos -} phytoplankton cells positions
#' \item \strong{others_nk2 -} other unidentified particle positions
#' }
#'
#' @description This function takes in a flowframe with debris removed and identifies
#'              phytoplankton cell population in the provided frame.
#'
#' @details The function uses the \code{\link[flowDensity]{getPeaks}} and
#'          \code{\link[flowDensity]{deGate}} functions in the \emph{flowDensity} package to
#'          identify peaks and identify cut-off points between these peaks. This function is
#'          not designed to be called in isolation, if called in isolation an error will be
#'          returned. It is preferably called on the results from \code{\link{debris_nc}} or
#'          \code{\link{debris_inc}} function. A graph with horizontal
#'          and vertical lines used in separating the populations is returned and
#'          if \emph{to_retain = "refined"}, a circle made of dashed lines is drawn around
#'          phytoplankton cell population points.
#'
#' @seealso \code{\link{bs5_nc}}
#'
#'
#' @importFrom utils capture.output
#' @export bs4_nc

bs4_nc <- function(bs4bs5, ch1, ch2, others, ph = 0.1, to_retain = c("refined", "potential")) {

    ch2_peaks <- flowDensity::getPeaks(bs4bs5, ch2, tinypeak.removal = ph)
    ch1_peaks <- flowDensity::getPeaks(bs4bs5, ch1, tinypeak.removal = ph)
    msg <- capture.output(flowDensity::deGate(bs4bs5, ch2, all.cuts = T))[1]

    if (length(ch2_peaks$Peaks) == 1) {

        ch2_cuts <- flowDensity::deGate(bs4bs5, ch2, use.percentile = T, percentile = 0.97)
        bs4_pot <- bs4bs5[which(bs4bs5@exprs[, ch2] < ch2_cuts), ]
        others_pot <- others[which(bs4bs5@exprs[, ch2] < ch2_cuts)]
        others_nk <- others[which(!(bs4bs5@exprs[, ch2] < ch2_cuts))]
        ptt <- "1"

    } else if (length(ch2_peaks$Peaks) == 2) {

        # mgrouch1 <- ch2_peaks$Peaks[2] - 0
        # if (mgrouch1 < 2.9) {
        #     ch2_cuts <- flowDensity::deGate(bs4bs5, ch2, use.percentile = T, percentile = 0.99)
        #     bs4_pot <- bs4bs5[which(bs4bs5@exprs[, ch2] < ch2_cuts), ]
        #     others_pot <- others[which(bs4bs5@exprs[, ch2] < ch2_cuts)]
        #     others_nk <- others[which(!(bs4bs5@exprs[, ch2] < ch2_cuts))]
        #     ptt <- "2a"
        #
        # } else {

            ch2_cuts <- flowDensity::deGate(bs4bs5, ch2, all.cuts = T)
            # BS4 is to the left of the minimum ch2_cuts
            bs4_pot <- bs4bs5[which(bs4bs5@exprs[, ch2] < ch2_cuts), ]
            others_pot <- others[which(bs4bs5@exprs[, ch2] < ch2_cuts)]
            others_nk <- others[which(!(bs4bs5@exprs[, ch2] < ch2_cuts))]
            ptt <- "2b"

        #}
    } else {

        ch2_cuts <- flowDensity::deGate(bs4bs5, ch2, use.percentile = T, percentile = 0.05)
        # min(ch2_peaks$Peaks) + 0.50*sd(bs4bs5@exprs[, ch2])
        bs4_pot <- bs4bs5[which(bs4bs5@exprs[, ch2] > ch2_cuts), ]
        others_pot <- others[which(bs4bs5@exprs[, ch2] > ch2_cuts)]
        others_nk <- others[which(!(bs4bs5@exprs[, ch2] > ch2_cuts))]
        ptt <- "3"

    }
    #
    abline(h = ch2_cuts, lty = 2, col = 2)
    # points(bs4_pot@exprs[, c(ch1, ch2)], pch = '.', col = 2) text(0, 2.5, paste('BS4', ptt, sep = ':'), col = 2)

    if (to_retain == "refined") {

        bs4s <- flowDensity::flowDensity(bs4_pot, channels = c(ch1, ch2), position = c(F, F),
                                         use.upper = c(T, F),
                                         upper = c(T, NA),
                                         use.percentile = c(F, T),
                                         percentile = c(NA, 0.7),
                                         ellip.gate = T)
        points(bs4s@filter, type = "l", col = 2, lwd = 2, lty = 4)
        text(mean(bs4s@filter[, 1]), mean(bs4s@filter[, 2]),
             "Syn", col = 2)
        # reduced flowframe for BS4
        bs4_reduced <- bs4_pot[which(!is.na(bs4s@flow.frame@exprs[, 1])), ]
        # positions of BS4s and other unidentified particles
        bs4_pos <- others_pot[which(!is.na(bs4s@flow.frame@exprs[, 1]))]
        others_nk2 <- others_pot[which(is.na(bs4s@flow.frame@exprs[, 1]))]

    } else if (to_retain == "potential") {

        text(mean(bs4_pot@exprs[, ch1]), mean(bs4_pot@exprs[, ch2]), "Syn", col = "red4")
        # reduced flowframe for BS4
        bs4_reduced <- bs4_pot
        # positions of BS4s and other unidentified particles
        others_nk <- others_nk
        bs4_pos <- others_pot
        others_nk2 <- NA

    } else stop("supply to_retain")

    return(list(syn_reduced = bs4_reduced, others_nk = others_nk,
                syn_pos = bs4_pos, others_nk2 = others_nk2))

}
