#' returns the position of the cells below, above or between estimated gates
#'
#' @param flowframe after debris are removed.
#' @param gates cut point between the identified clusters
#' @param ch gated channel
#'
#'
#'
#' @export row_numbers

<<<<<<< HEAD
row_numbers <- function(flowframe, gates, ch) {
=======
row_numbers <- function(flowfile, gates, ch) {
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a

  if(length(gates) == 1) {

    # 2 populations detected
<<<<<<< HEAD
    others_pot <- which(flowframe@exprs[, ch] <= gates) # the row_number of the guys
    others_pot_compl <- which(flowframe@exprs[, ch] > gates) #setdiff(others, others_pot) # row number of the other guys
=======
    others_pot <- which(flowfile@exprs[, ch] <= gates) # the row_number of the guys
    others_pot_compl <- which(flowfile@exprs[, ch] > gates) #setdiff(others, others_pot) # row number of the other guys
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a
    return(list(row_num_pop1 = others_pot, row_num_pop2 = others_pot_compl))

  } else {

    #more than 2 populations detected
    others_pot <- vector("list", length = length(gates))
    n <- length(gates) + 1

    for(i in 1:n) {

      if(i == 1) {

<<<<<<< HEAD
        oth1 <- which(flowframe@exprs[, ch] <= gates[i])
=======
        oth1 <- which(flowfile@exprs[, ch] <= gates[i])
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a
        others_pot[[i]] <- oth1

      } else if(i == n) {

<<<<<<< HEAD
        oth1 <- which(flowframe@exprs[, ch] > gates[i-1])
=======
        oth1 <- which(flowfile@exprs[, ch] > gates[i-1])
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a
        others_pot[[i]] <- oth1

      } else {

<<<<<<< HEAD
        oth2 <- which((flowframe@exprs[, ch] > gates[i-1]) & (flowframe@exprs[, ch] <= gates[i]))
=======
        oth2 <- which((flowfile@exprs[, ch] > gates[i-1]) & (flowfile@exprs[, ch] <= gates[i]))
>>>>>>> 74ebbc1f0b57e1f0410570dfc73e0a04247a1d3a
        others_pot[[i]] <- oth2

      }

    }

    names(others_pot) <- paste0("row_num_pop", 1:length(others_pot))
    return(others_pot)

  }


}
