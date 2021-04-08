#' gates out or assign indicators to debris particle based on their chlorophyll 
#' expression.
#'
#' @param flowframe flowframe with debris and other cells.
#' @param ch_chlorophyll first flowcytometer channel that can be used to 
#'                       separate debris from the rest, e.g. "RED.B.HLin".
#' @param ch_p2 second flowcytometer channel use for plotting
#'              from the rest, e.g. "YEL.B.HLin"
#' @param ph the minimum peak height that should be considered. 
#'           This aids the removal of tiny peaks. Defaults to 0.1
#' @param n_sd number of standard deviations away from peak should be 
#'             considered to filter out debris
#'
#' @return list containing; \itemize{
#' \item \strong{syn - flowframe containing non-debris particles}
#' \item \strong{deb_pos - position of particles that are debris}
#' \item \strong{syn_pos - position of particles that are not debris}
#' }
#'
#' @description The function takes in a flowframe and identifies debris 
#'              contained in the provided flowframe.
#'
#' @details The function uses the \code{\link[flowDensity]{getPeaks}} and
#'          \code{\link[flowDensity]{deGate}} functions in the flowDensity 
#'          package to
#'          identify peaks in ch_chlorophyll, and identify cut-off points 
#'          #between these peaks.
#'          A plot of both channels supplied with horizontal line separating
#'          debris from other cell populations is also returned.
#'
#' @examples
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#'                   package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::noNA(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noNeg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#' c('SSC.W', 'TIME'))
#' cells_nonmargin <- cellMargin(flowframe = flowfile_logtrans, 
#' Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' debrisNc(flowframe = reducedFlowframe(cells_nonmargin), 
#'           ch_chlorophyll = "RED.B.HLin",
#'           ch_p2 = "YEL.B.HLin",
#'           ph = 0.05)
#'
#' @export debrisNc

debrisNc <- function(flowframe, ch_chlorophyll, 
                     ch_p2, 
                      ph = 0.09, n_sd = 2) {

    # debris gating
    ch_chlorophyll_peaks <- flowDensity::getPeaks(flowframe, ch_chlorophyll)

    #length of the peak vectors
    lp <- length(ch_chlorophyll_peaks$Peaks[ch_chlorophyll_peaks$Peaks >= 0])

    #check if the peak heights are above a certain threshold, i.e. check 
    #for tinypeaks
    #if(sum(ch_chlorophyll_peaks$P.h[ch_chlorophyll_peaks$Peaks >= 0] > ph)
    #   < lp) {

        ch_chlorophyll_peaks2 <- flowDensity::getPeaks(flowframe, 
                                                       ch_chlorophyll,
                                                       tinypeak.removal = ph)


    #} 
    #else {
    # 
    #     ch_chlorophyll_peaks2 <- flowDensity::getPeaks(flowframe, 
    #                                                    ch_chlorophyll,
    #                                                    tinypeak.removal = ph)
    # 
    # }

    #recompute length of the peak vectors
    lp2 <- length(ch_chlorophyll_peaks2$Peaks[ch_chlorophyll_peaks2$Peaks >= 0])
    pks <- ch_chlorophyll_peaks2$Peaks[ch_chlorophyll_peaks2$Peaks >= 0]

    pkks <- ifelse(lp2 >= 2, min(pks), pks)
    pks_ci <- pkks + c(-n_sd, n_sd) * sd(exprs(flowframe[, ch_chlorophyll]))
    
    if (lp2 >= 2 & (min(pks_ci) <= 0 )) {

        # at least 2 peaks with confience interval for smallest 
        # peak containing 0
        deb_cut <- flowDensity::deGate(flowframe, ch_chlorophyll, 
                                       all.cuts = TRUE,
                                       tinypeak.removal = ph)[1]

    } else if(lp2 >= 2 & (min(pks_ci) > 0 )) {

        # at least 2 peaks with confience interval for 
        # smallest peak not containing 0
        deb_cut <- flowDensity::deGate(flowframe, ch_chlorophyll, 
                                       use.percentile = TRUE,
                                       percentile = 0.05,
                                       tinypeak.removal = ph)

    } else if (lp2 < 2 & (min(pks_ci) <= 0 )) {

        #one peak, but confidence interval for the peak contains 0
        deb_cut <- flowDensity::deGate(flowframe, ch_chlorophyll, 
                                       use.percentile = TRUE,
                                       percentile = 0.90,
                                       tinypeak.removal = ph)

    } else {

            deb_cut <- flowDensity::deGate(flowframe, ch_chlorophyll, 
                                           use.percentile = TRUE,
                                           percentile = 0.10,
                                           tinypeak.removal = ph)

    }

    # plotting
    #plott1 <- ggplotDens(flowframe, c(ch_chlorophyll, ch_p2)) +
    #    geom_vline(xintercept = deb_cut, color = "red", linetype = "dashed") +
    #    geom_text(aes(x = mean(flowframe@exprs[which(flowCore::exprs(flowframe)
    #[, ch_chlorophyll] <= deb_cut),
    #                                           ch_chlorophyll]),
    #                  y = mean(flowframe@exprs
    #[which(flowCore::exprs(flowframe)[, ch_p2] <=
    #                 deb_cut), ch_p2])), inherit.aes = FALSE,
    #              label = paste0("Debris"), colour = "blue", size = 6)
    syn <- flowframe[which(flowCore::exprs(flowframe)[, ch_chlorophyll] > 
                               deb_cut), ]
    other_pos <- which(flowCore::exprs(flowframe)[, ch_chlorophyll] > deb_cut)
    deb_pos <- which(flowCore::exprs(flowframe)[, ch_chlorophyll] <= deb_cut)

    ret_result <- DebrisFilter(
                       fullflowframe = flowframe,
                       reducedflowframe = syn,
                       deb_pos = deb_pos,
                       syn_all_pos = other_pos,
                       deb_cut = deb_cut,
                       ch_chlorophyll = ch_chlorophyll,
                       ch_p2 = ch_p2
                      )
    #attr(ret_result, 'class') <- c('cyanoFilter', 'Debris')
    #cyanoFilterDebris(ret_result)
    return(ret_result)
}
