#' Removes or assign indicators to margin events.
#'
#' @param flowframe Flowframe containing margin events to be filtered out
#' @param Channel The channel on which margin events are. Defaults to SSC.W (side scatter width)
#' @param type    The method to be used in gating out the margin cells. Can either be 'manual' where
#'                user supplies a cut off point on the channel, 1 = not margin 0 = margin
#' @param cut sould not be NULL if type = 'manual'
#' @param y_toplot channel on y-axis of plot with \emph{Channel} used to gate out margin events
#'
#' @return an object of cyanoFilter:MarginEvents class containing; \itemize{
#' \item \strong{reducedflowframe -} flowframe without margin events
#' \item \strong{fullflowframe -} flowframe with an Margin.Indicator added as an extra column added to the expression matrix
#' to indicate which particles are margin events. 1 = not margin event, 0 = margin event
#' \item \strong{N_margin -} number of margin events recorded
#' \item \strong{N_cell -} numner of non-margin events
#' \item \strong{N_particle -} is the number of particles in total, i.e. N_cell + N_margin
#' }
#'
#' @description The function identifies margin events, i.e. cells that are too large for the flow cytometer to measure.
#'
#' @details Users can either supply a cut-off point along the channel describing particle width
#'          or allow the function to estimate the cut-off point using the
#'          \code{\link[flowDensity]{deGate}} function from the \emph{flowDensity} package.
#'          A plot of channel against "FSC.HLin" is provided with a vertical line showing the
#'          cut-off point separating margin events from other cells.
#'
#' @examples
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) #FCS file contains only one data object
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- lnTrans(x = flowfile_noneg, c('SSC.W', 'TIME'))
#' cellmargin(flowframe = flowfile_logtrans, Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#'
#'
#' @importFrom methods new
#' @export cellmargin

cellmargin <- function(flowframe, Channel = "SSC.W",
                       type = c("manual", "estimate"),
                       cut = NULL, y_toplot = "FSC,HLin") {

    dvarMetadata <- flowframe@parameters@varMetadata
    ddimnames <- flowframe@parameters@dimLabels
    describe <- flowframe@description

    if (type == "manual" & !is.null(cut)) {
        margin.ind <- ifelse(flowCore::exprs(flowframe)[, Channel] <= cut, T, F)
        n_margin <- sum(margin.ind == F)

        plt_margin <- ggplotDens(flowframe, c(Channel, y_toplot)) +
            geom_vline(xintercept = cut, linetype = 'dashed')

        # constructing full flow frame with both margin and non-margin events
        exx <- as.matrix(cbind(flowCore::exprs(flowframe), as.numeric(margin.ind)))  #1=not margin, 0=margin
        colnames(exx)[ncol(exx)] <- "Margin.Indicator"

        # constructing the annotated data frame for the parameter
        ddata <- data.frame(rbind(flowframe@parameters@data, c("Margin.Indicator", "Margin Indicator", 1, 0, 1)))
        paraa <- Biobase::AnnotatedDataFrame(data = ddata, varMetadata = dvarMetadata, dimLabels = ddimnames)
        row.names(ddata) <- c(row.names(flowframe@parameters@data), paste0("$", "P", ncol(exx)))

        fflowframe <- flowCore::flowFrame(exprs = exx, parameters = paraa, description = describe)

        # constructing flowframe with only the non margin events
        exx2 <- exx[exx[, "Margin.Indicator"] == 1, ]
        rflowframe <- flowCore::flowFrame(exprs = exx2, parameters = paraa, description = describe)

    } else if (type == "estimate") {

        infl_point1 <- flowDensity::deGate(flowframe, Channel, bimodal = T) #upper boundary
        infl_point2 <- flowDensity::deGate(flowframe, Channel, use.upper = T, upper = F) #lower boundary

        # plotting
        plt_margin <- ggplotDens(flowframe, c(Channel, y_toplot)) +
            geom_vline(xintercept = infl_point1, linetype = 'dashed') +
            geom_vline(xintercept = infl_point2,
                       linetype = 'dashed')

        margin.ind <- ifelse(flowCore::exprs(flowframe)[, Channel] > infl_point2 &
                             flowCore::exprs(flowframe)[, Channel] < infl_point1, T, F)  #FALSE = Margin Event
        n_margin <- sum(margin.ind == F)  #part of output

        # constructing full flow frame with both margin and non-margin events
         exx <- as.matrix(cbind(flowCore::exprs(flowframe), as.numeric(margin.ind)))  #1=not margin, 0=margin
         colnames(exx)[ncol(exx)] <- "Margin.Indicator"  #1 = not amrgin, 0 = margin
        #
         # constructing the annotated data frame for the parameter
         ddata <- data.frame(rbind(flowframe@parameters@data, c("Margin.Indicator", "Margin Indicator", 1, 0, 1)))
         paraa <- Biobase::AnnotatedDataFrame(data = ddata, varMetadata = dvarMetadata,
                                              dimLabels = ddimnames)
         row.names(ddata) <- c(row.names(flowframe@parameters@data),
                               paste("$P", length(row.names(flowframe@parameters@data))+1,
                                     sep = "")
                               )

         fflowframe <- flowCore::flowFrame(exprs = exx,
                                           parameters = paraa,
                                           description = describe)

        # constructing flowframe with only the non margin events
        exx2 <- exx[exx[, "Margin.Indicator"] == 1, ]
        rflowframe <- flowCore::flowFrame(exprs = exx2, parameters = paraa, description = describe)
    } else stop("Error: check your inputs")

    ret_result <- list(reducedflowframe = rflowframe,
                       fullflowframe = fflowframe,
                       N_margin = n_margin,
                       N_nonmargin = flowCore::nrow(rflowframe),
                       N_particle = flowCore::nrow(fflowframe),
                       plot = plt_margin)
    #cyanoFilterMargin(ret_result)
    attr(ret_result, 'class') <- c('cyanoFilter', 'marginEvents')

    return(ret_result)
}
