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
#' celldebris_nc(bs4bs5 = flowfile, channel1 = "RED.B.HLin", channel2 = "YEL.B.HLin",
#' interest = "Both", to_retain = "refined")
#' }
#'
#'
#' @export celldebris_nc


celldebris_nc <- function(flowframe, channel1 = "RED.B.HLin", channel2 = "YEL.B.HLin",
                          interest = c("BS4", "BS5", "Both"),
                          to_retain = c("refined", "potential")) {


    if (interest == "Both") {

        ddata <- data.frame(rbind(flowframe@parameters@data, c("BS4BS5.Indicator", "BS4BS5.Indicator", 1, 0, 1)))
        row.names(ddata) <- c(row.names(flowframe@parameters@data), "$P13")

    } else if (interest == "BS4" | interest == "BS5") {

        ddata <- data.frame(rbind(flowframe@parameters@data, c("BS4BS5.Indicator", "BS4BS5.Indicator", 1, 0, 1)))
        row.names(ddata) <- c(row.names(flowframe@parameters@data), "$P13")

    } else stop("incorrect options supplied, check interest")

    dvarMetadata <- flowframe@parameters@varMetadata
    ddimnames <- flowframe@parameters@dimLabels
    paraa <- Biobase::AnnotatedDataFrame(data = ddata, varMetadata = dvarMetadata, dimLabels = ddimnames)
    describe <- flowframe@description
    BS4BS5_ind <- rep(NA, flowCore::nrow(flowframe))

    # gating debris
    if (interest == "Both") {

        debris <- debris_inc(flowframe = flowframe, p1 = channel1, p2 = channel2)

    } else debris <- debris_nc(flowframe = flowframe, p1 = channel1, p2 = channel2)

    # BS4s and BS5s
    bs4bs5 <- debris$bs4bs5
    # BS4BS5 positions in the expression matrix
    bs4bs5_pos <- debris$bs4bs5_pos

    # Fill debris positions with 0
    BS4BS5_ind[debris$deb_pos] <- 0

    if (interest == "BS4") {

        bs4s <- bs4_nc(bs4bs5, p1 = channel1, p2 = channel2, others = bs4bs5_pos, to_retain = to_retain)

        ### indicators for each cell type
        BS4BS5_ind[bs4s$bs4_pos] <- 1  #BS4
        BS4BS5_ind[bs4s$others_nk] <- 2  #others not identified
        BS4BS5_ind[bs4s$others_nk2] <- 2  #others not identified

        # number of BS4
        n_cyano <- sum(BS4BS5_ind == 1)

    } else if (interest == "BS5") {

        bs5s <- bs5_nc(bs4bs5, p1 = channel1, p2 = channel2, others = bs4bs5_pos, to_retain = to_retain)

        ### indicators for each cell type
        BS4BS5_ind[bs5s$bs5_pos] <- 1  #BS5
        BS4BS5_ind[bs5s$others_nk] <- 2  #others not identified
        BS4BS5_ind[bs5s$others_nk2] <- 2  #others not identified

        # number of BS5
        n_cyano <- sum(BS4BS5_ind == 1)

    } else if (interest == "Both") {

        # BS4
        bs4s <- bs4_nc(bs4bs5, p1 = channel1, p2 = channel2, others = bs4bs5_pos, to_retain = to_retain)
        # BS5
        bs5s <- bs5_nc(bs4bs5, p1 = channel1, p2 = channel2, others = bs4bs5_pos, to_retain = to_retain)

        ### indicators for each cell type
        BS4BS5_ind[bs4s$bs4_pos] <- 1  #BS4
        BS4BS5_ind[bs5s$bs5_pos] <- 2  #BS5
        BS4BS5_ind[is.na(BS4BS5_ind)] <- 3  #others not identified

        # number of BS4, number of BS5
        n_cyano <- c(BS4 = sum(BS4BS5_ind == 1), BS5 = sum(BS4BS5_ind == 2))


    } else stop("Supply interest")

    ### Full flowframe Forming a new expression matrix for the full flowframe with indicator added for BS4 or BS5
    nexp_mat <- as.matrix(cbind(flowCore::exprs(flowframe), as.numeric(BS4BS5_ind)))
    # giving a name to the newly added column to the expression matrix
    colnames(nexp_mat)[length(colnames(nexp_mat))] <- "BS4BS5.Indicator"
    # full flow frame with indicator for particly type
    fflowframe <- methods::new("flowFrame", exprs = nexp_mat, parameters = paraa, description = describe)

    ### reduced flowframe
    if (interest == "Both") {

        rframe <- fflowframe[fflowframe@exprs[, 13] %in% c(1, 2), ]

    } else {

        rframe <- fflowframe[fflowframe@exprs[, 13] == 1, ]
    }

    # number of debris
    n_debris <- sum(BS4BS5_ind == 0)

    return(list(fullframe = fflowframe, reducedframe = rframe, Cell_count = n_cyano, Debris_Count = n_debris))

}
