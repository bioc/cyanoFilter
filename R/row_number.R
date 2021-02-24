#' returns the position of the cells below, above or between estimated gates
#'
#' @param flowframe after debris are removed.
#' @param gates cut point between the identified clusters
#' @param ch gated channel
#'
#' @return a numeric vector
#' 
#' @examples
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#'                   package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#' c('SSC.W', 'TIME'))
#' oneDgate(flowfile, 'RED.B.HLin')
#'
#' @export row_numbers

row_numbers <- function(flowframe, gates, ch) {

  if(length(gates) == 1) {

    # 2 populations detected
    # the row_number of the guys
    others_pot <- which(flowframe@exprs[, ch] <= gates) 
    #setdiff(others, others_pot) # row number of the other guys
    others_pot_compl <- which(flowframe@exprs[, ch] > gates) 
    return(list(row_num_pop1 = others_pot, row_num_pop2 = others_pot_compl))

  } else {

    #more than 2 populations detected
    others_pot <- vector("list", length = length(gates))
    n <- length(gates) + 1

    for(i in seq_len(n)) {

      if(i == 1) {

        oth1 <- which(flowframe@exprs[, ch] <= gates[i])
        others_pot[[i]] <- oth1

      } else if(i == n) {

        oth1 <- which(flowframe@exprs[, ch] > gates[i-1])
        others_pot[[i]] <- oth1

      } else {

        oth2 <- which((flowframe@exprs[, ch] > gates[i-1]) & 
                       (flowframe@exprs[, ch] <= gates[i]))
        others_pot[[i]] <- oth2

      }

    }

    names(others_pot) <- paste0("row_num_pop", seq_along(others_pot))
    return(others_pot)

  }

}
