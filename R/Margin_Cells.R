#' Removes or assign indicators to margin events.
#'
#' @param flowframe Flowframe containing margin events to be filtered out
#' @param Channel The channel on which margin events are. Defaults to 
#'                SSC.W (side scatter width)
#' @param type    The method to be used in gating out the margin cells. 
#'                Can either be 'manual' where
#'                user supplies a cut off point on the channel, 
#'                1 = not margin 0 = margin
#' @param cut sould not be NULL if type = 'manual'
#' @param y_toplot channel on y-axis of plot with \emph{Channel} used to gate 
#'                 out margin events
#'
#' @return an object of class MarginEvents class containing slots; \itemize{
#' \item \strong{reducedflowframe -} flowframe without margin events
#' \item \strong{fullflowframe -} flowframe with an Margin.Indicator added as 
#'              an extra column added to the expression matrix
#'              to indicate which particles are margin events. 
#'              1 = not margin event,  0 = margin event
#' \item \strong{N_margin -} number of margin events recorded
#' \item \strong{N_cell -} numner of non-margin events
#' \item \strong{N_particle -} is the number of particles in total, 
#'              i.e. N_cell + N_margin
#' }
#'
#' @description The function identifies margin events, 
#'              i.e. cells that are too large for the flow 
#'              cytometer to measure.
#'
#' @details Users can either supply a cut-off point along the channel 
#'          describing particle width
#'          or allow the function to estimate the cut-off point using the
#'          \code{\link[flowDensity]{deGate}} function from the 
#'          \emph{flowDensity} package.
#'          A plot of channel against "FSC.HLin" is provided with a vertical 
#'          line showing the cut-off point separating margin events 
#'          from other cells.
#'
#' @examples
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#'                  package = "cyanoFilter",
#'                  mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1)
#' flowfile_nona <- cyanoFilter::noNA(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noNeg(x = flowfile_nona)
#' flowfile_logtrans <- lnTrans(x = flowfile_noneg, c('SSC.W', 'TIME'))
#' cellMargin(flowframe = flowfile_logtrans, Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#'
#' @importFrom methods new
#' @export cellMargin

cellMargin <- function(flowframe, Channel = "SSC.W",
                       type = c("manual", "estimate"),
                       cut = NULL, y_toplot = "FSC,HLin") {

    dvarMetadata <- flowCore::varMetadata(flowCore::parameters(flowframe))
    ddimnames <- Biobase::dimLabels(flowCore::parameters(flowframe))
    describe <- flowCore::keyword(flowframe)
    pd <- pData(flowCore::parameters(flowframe))

    if (type == "manual" & !is.null(cut)) {
        margin.ind <- ifelse(flowCore::exprs(flowframe)[, Channel] <= cut, 
                             TRUE, FALSE)
        n_margin <- sum(margin.ind == FALSE)

        #plt_margin <- ggplotDens(flowframe, c(Channel, y_toplot)) +
        #    geom_vline(xintercept = cut, linetype = 'dashed')

        # constructing full flow frame with both margin and non-margin events
        #1=not margin, 0=margin
        exx <- as.matrix(cbind(flowCore::exprs(flowframe), 
                               as.numeric(margin.ind)))  
        colnames(exx)[ncol(exx)] <- "Margin.Indicator"

        # constructing the annotated data frame for the parameter
        ddata <- data.frame(rbind(pd, 
                                  c("Margin.Indicator", "Margin Indicator", 
                                    1, 0, 1)))
        paraa <- Biobase::AnnotatedDataFrame(data = ddata, 
                                             varMetadata = dvarMetadata, 
                                             dimLabels = ddimnames
                                             )
        row.names(ddata) <- c(row.names(pd), 
                              paste0("$", "P", ncol(exx))
                              )

        fflowframe <- flowCore::flowFrame(exprs = exx, parameters = paraa, 
                                          description = describe)

        # constructing flowframe with only the non margin events
        exx2 <- exx[exx[, "Margin.Indicator"] == 1, ]
        rflowframe <- flowCore::flowFrame(exprs = exx2, parameters = paraa, 
                                          description = describe)

    } else if (type == "estimate") {
        #upper boundary
        infl_point1 <- flowDensity::deGate(flowframe, 
                                           Channel, 
                                           bimodal = TRUE) 
        #lower boundary
        infl_point2 <- flowDensity::deGate(flowframe, Channel, 
                                           use.upper = TRUE, 
                                           upper = FALSE) 

        # plotting
        #plt_margin <- ggplotDens(flowframe, c(Channel, y_toplot)) +
        #    geom_vline(xintercept = infl_point1, linetype = 'dashed') +
        #    geom_vline(xintercept = infl_point2,
        #               linetype = 'dashed')

        margin.ind <- ifelse(flowCore::exprs(flowframe)[, Channel] > 
                                 infl_point2 &
                             flowCore::exprs(flowframe)[, Channel] < 
                                 infl_point1, TRUE, FALSE) 
        n_margin <- sum(margin.ind == FALSE)  #part of output

        # constructing full flow frame with both margin and non-margin events
        #1=not margin, 0=margin
         exx <- as.matrix(cbind(flowCore::exprs(flowframe), 
                                as.numeric(margin.ind)))  
         colnames(exx)[ncol(exx)] <- "Margin.Indicator"  
        #
         # constructing the annotated data frame for the parameter
         ddata <- data.frame(rbind(pd, 
                                   c("Margin.Indicator", "Margin Indicator", 
                                     1, 0, 1)))
         paraa <- Biobase::AnnotatedDataFrame(data = ddata, 
                                              varMetadata = dvarMetadata,
                                              dimLabels = ddimnames)
         row.names(ddata) <- c(row.names(pd),
                               paste("$P", 
                               length(row.names(pd))+1,
                                     sep = "")
                               )

         fflowframe <- flowCore::flowFrame(exprs = exx,
                                           parameters = paraa,
                                           description = describe)

        # constructing flowframe with only the non margin events
        exx2 <- exx[exx[, "Margin.Indicator"] == 1, ]
        rflowframe <- flowCore::flowFrame(exprs = exx2, parameters = paraa, 
                                          description = describe)
    } else stop("check your inputs")

    ret_result <- MarginEvents(fullflowframe = fflowframe,
                       reducedflowframe = rflowframe,
                       N_margin = n_margin,
                       N_nonmargin = flowCore::nrow(rflowframe),
                       N_particle = flowCore::nrow(fflowframe),
                       Channel = Channel,
                       y_toplot = y_toplot,
                       cut = ifelse(type == "manual", 
                                    cut, 
                                    infl_point1)
                      )
    
    return(ret_result)
}
