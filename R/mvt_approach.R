#' function to cluster BS4, BS5 and Debris in a flowfile using an EM style algorithm.
#' are used and the data is clustered into 5 clusters (automatically reduces this number)
#'
#' @param flowfile flowframe with.
#' @param channels channels to use for the clustering
#' @param mu pre-specified mean matrix for the clusters. Number of rows should equal ncluster and number
#'           of columns should equal length(channels). Defaults to NULL and can be computed from the data internally.
#' @param sigma pre-specified list of variance-covariace matrix for the clusters. Each element of the list should contain a square matrix
#'        of variance-covariance matrix with length equal ncluster.
#' @param  ncluster number of cluster to desired
#' @return list containing; \itemize{
#' \item \strong{percentages -}
#' \item \strong{mus -}
#' \item \strong{sigmas -}
#' \item \strong{result -}
#' }
#'
#' @examples
#' \dontrun{
#' cyano_emclustering(flowfile = flowfile, channel1 = c("RED.B.HLin", "YEL.B.HLin", "FSC.HLin", "RED.R.HLin"),
#' ncluster = 5, min.itera = 20)
#' }
#'
#'
#' @export cyano_emclustering




cyano_emclustering <- function(flowfile, channels, mu = NULL, sigma = NULL, ncluster = 5,
                       min.itera = 20) {

      data <- flowCore::exprs(flowfile)


      if(is.null(mu)) {
          #sapply(, )
          #mu <- matrix(, nrow = ncluster, ncol = length(channels))

      } else {
          mu <- mu
      }

      if(is.null(sigma)) {
          sigma <- vector("list", length = ncluster)
          sigma <- lapply(sigma, diag(var(data), ncol = length(channels),
                                      nrow = length(channels)))
      } else {
        sigma <- sigma
      }

      # probability of each cell to belong to each cluster
      lambda <- matrix(NA, nrow = nrow(data), ncol = ncluster)
      tau <- rep(1/ncluster, ncluster) # wheight of each cluster

    i <- 0

    while(TRUE) {
      i <- i + 1
      # to check progress
      tau.old <- tau
      mu.old <- mu
      # compute the probability of each cell for each cluster (E step)
      for (p in 1:ncluster) {
        lambda[,p] <- tau[p] * mvnorm(data, mu[,p], sig[[p]])
      }
      if(anyNA(lambda)) { # one cluster has size 0, repeat algorithm
        #clusters.pres <- clusters.pres[!is.na(lambda[1,clusters.pres])]
        return( clustering(data, ncluster = ncluster - 1))
    }

    # each cell should have probability 1
    lambda <- lambda/rowSums(lambda)

    # probabilities/sizes of the different clusters
    tau <- colMeans(lambda)

    # adapt the mean and the std for the different clusters (M - step)
    for(p in 1:ncluster) {
        mu[,p] <- colSums(lambda[, p] * data) / sum(lambda[, p])
        sig[[p]] <- t(lambda[,p] * sweep(data, 2, mu[,p])) %*% sweep(data, 2, mu[, p]) /
                    sum(lambda[, p])
    }
    # check whether local maxima has been achieved
    rel.diff.mu <- max((abs(mu - mu.old) / mu))
    rel.diff.tau <- max((abs(tau - tau.old) / tau))
    #if (min(c(rel.diff.mu, rel.diff.tau)) < 0) print(c(i,rel.diff.mu, rel.diff.tau))
    if (i > min.itera & rel.diff.mu < 1e-3 & rel.diff.tau < 1e-3) {
      assigned <- proportions(lambda)
      print(c(i, rel.diff.mu, rel.diff.tau))
      return(list(percentages = tau, mus = mu, sigmas = sig,
                  assigned = assigned))
      }
    }

    ddata <- data.frame()
    for(i in 1:ncol(lambda)) {

      ddata <- data.frame(rbind(ddata, flowfile@parameters@data,
                       c(paste("Cluster_Prob", i, "_"), paste("Cluster_Prob", i, sep = "_"), 1, 0, 1))
      )
    }


    paraa <- Biobase::AnnotatedDataFrame(data = ddata, varMetadata = dvarMetadata, dimLabels = ddimnames)
    describe <- flowframe@description
    row.names(ddata) <- c(row.names(flowframe@parameters@data), "$P13", "$P14", "$P15", "$P16", "$P17")

    ### Full flowframe Forming a new expression matrix for the full flowframe with indicator added for BS4 or BS5
    nexp_mat <- as.matrix(cbind(flowCore::exprs(flowframe), as.numeric(lambda)))
    # giving a name to the newly added column to the expression matrix
    colnames(nexp_mat)[13:length(colnames(nexp_mat))] <- paste("Cluster_Prob", 13:17, sep = "_")

    # full flow frame with indicator for particly type
    fflowframe <- methods::new("flowFrame", exprs = nexp_mat, parameters = paraa, description = describe)

    return(list(percentages = tau, mus = mu, sigmas = sig,
                  result = fflowframe))
}
