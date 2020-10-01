
#' produces a scatter plot of the expression matrix of a flowframe. Note that, it takes some time to display the plot.
#'
#' @param x flowframe to be plotted
#' @param notToPlot column in expression matrix not to be plotted
#' @param ... other arguments. Not used at the moment
#' @import ggplot2 GGally flowCore
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics plot.default pairs.default abline points panel.smooth pairs smoothScatter text
#' @export pairs_plot

pairs_plot <- function(x, notToPlot = c("TIME"), ...) {

  toplot <- setdiff(flowCore::colnames(x), notToPlot)
  col.palette <- colorRampPalette(c("white", "blue", "cyan", "green", "orange", "red"),
                                  space = "Lab")

  if(class(x) == 'flowFrame') {


    pairs.default(flowCore::exprs(x)[, toplot], pch = ".",
                  panel = function(...) smoothScatter(..., nrpoints = 0,
                                                      colramp = col.palette,
                                                      add = TRUE), gap = 0.2,
                  main = flowCore::identifier(x))
  } else if(!(class(x)[2] == 'Debris')) {


    pairs.default(flowCore::exprs(x$fullflowframe)[, toplot], pch = ".",
                  panel = function(...) smoothScatter(..., nrpoints = 0,
                                                      colramp = col.palette,
                                                      add = TRUE), gap = 0.2,
                  main = flowCore::identifier(x))

  } else if(class(x)[2] == 'Debris') {


    pairs.default(flowCore::exprs(x$syn)[, toplot], pch = ".",
                  panel = function(...) smoothScatter(..., nrpoints = 0,
                                                      colramp = col.palette,
                                                      add = TRUE), gap = 0.2,
                  main = flowCore::identifier(x))

  } else stop('Error: object supplied not supported')

}



#' plots two channels of a flowframe.
#'
#' @param flowfile flowframe to be plotted
#' @param channels a character vector of length 2, must contain channel names in the flowfile.
#' @param ... not used at the moment
#'
#' @examples
#' \donttest{
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) #FCS file contains only one data object
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, c('SSC.W', 'TIME'))
#' ggplotDens(flowfile_logtrans,
#'            channels = c("FSC.HLin", "RED.R.HLin"))
#' }
#' @return a ggplot object
#'
#' @export ggplotDens

ggplotDens <- function(flowfile, channels, ...) {

  plotdata <- as.data.frame(flowfile@exprs[, channels])
  d <- densCols(plotdata, colramp = colorRampPalette(c("white", "blue", "cyan",
                                                       "green", "orange", "red"),
                                                     space = "Lab"))

  plotdata <- cbind(plotdata, d = d)
  ggplot(data = plotdata, aes(x = get(channels[1], plotdata), y = get(channels[2], plotdata), color = d)) +
    geom_point(pch = ".", size = 25) +
    theme_minimal() +
    scale_color_identity() +
    theme(axis.line.x = element_line(color = "black", size = 0.8),
          axis.line.y = element_line(color = "black", size = 0.8),
          plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          axis.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          axis.text = element_text(hjust = 0.5, size = 15, face = "bold")) +
    labs(x = channels[1], y = channels[2]) +
    ggtitle(identifier(flowfile))

}


#' plots two channels of a flowframe with different colors for clusters identified.
#'
#' @param flowfile flowframe to be plotted
#' @param channels a character vector of length 2, must contain channel names in the flowfile.
#' @param group cluster groups. must be equal to the number of particles in the flow cytometer.
#' @param ... not used at the moment
#'
#' @return a ggplot object
#'
#' @export ggplotDens2

ggplotDens2 <- function(flowfile, channels, group, ...) {

  plotdata <- as.data.frame(flowfile@exprs[, c(channels, group)])

  plotdata_mean <- summary.cyanoFilter(flowfile, channels = channels,
                                       cluster_var = group, summary = "mean")
  Group <- unique(flowfile@exprs[, group])
  hulls <- lapply(Group, function(x) {

    pd <- plotdata[plotdata[, group] == x, ]
    pd[chull(pd[, 1:2]), ]

  })

  hulldata <- do.call(rbind.data.frame, hulls)

  d <- densCols(plotdata, colramp = colorRampPalette(c("white", "blue", "cyan",
                                                       "green", "orange", "red"),
                                                     space = "Lab"))
  plotdata <- cbind(plotdata, d = d)
  ggplot(data = plotdata, aes(x = get(channels[1], plotdata), y = get(channels[2], plotdata), color = d)) +
    geom_point(pch = ".", size = 35) +
    theme_minimal() +
    scale_color_identity() +
    theme(axis.line.x = element_line(color = "black", size = 0.8),
          axis.line.y = element_line(color = "black", size = 0.8),
          plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          axis.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          axis.text = element_text(hjust = 0.5, size = 15, face = "bold")) +
    labs(x = channels[1], y = channels[2]) +
    ggtitle(identifier(flowfile)) +
    geom_polygon(data = hulldata, aes(x = get(channels[1], hulldata),
                                      y = get(channels[2], hulldata)),
                 linetype = "dashed", color = get(group, hulldata),
                 group = get(group, hulldata),
                 inherit.aes = FALSE, fill = NA,
                 size = 1.5) +
    geom_text(data = plotdata_mean, aes(x = get(channels[1], plotdata_mean),
                                        y = get(channels[2], plotdata_mean),
                                        label = paste("Cluster", plotdata_mean[, 1], sep = "-")
    ),
    size = 5,
    color = get(group, plotdata_mean),
    inherit.aes = FALSE)

}


#' produces a scatter plot of the expression matrix of the flowframe. If a cluster variable is given,
#' it assigns different colors to the clusters.
#'
#' @param flowfile flowframe to be plotted
#' @param channels a character vector of length 2 or more. It must contain channel names in the flowfile.
#' @param group cluster groups. It must be equal to the number of particles in the flowfile. If group is
#'              null cluster boundaries are not drawn.
#' @param notToPlot columns not to plot. This is especially useful for for plotting all columns in a
#' @param ... not used at the moment
#'
#' @examples \donttest{
#'
#'  # example without clustering
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) #FCS file contains only one data object
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- lnTrans(x = flowfile_noneg, c('SSC.W', 'TIME'))
#' ggpairsDens(flowfile = flowfile_logtrans,
#'             channels = c("FSC.HLin", "RED.R.HLin", "RED.B.HLin", "NIR.R.HLin"))
#'
#'  # example plot after clustering
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", package = "cyanoFilter",
#'                              mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) #FCS file contains only one data object
#' flowfile_nona <- cyanoFilter::nona(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, c('SSC.W', 'TIME'))
#' cells_nonmargin <- cellmargin(flowframe = flowfile_logtrans, Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' cells_nodebris <- debris_nc(flowframe = cells_nonmargin$reducedflowframe,
#'                             ch_chlorophyll = "RED.B.HLin",
#'                             ch_p2 = "YEL.B.HLin",
#'                             ph = 0.05)
#' fin <- phyto_filter(flowfile = cells_nodebris$syn,
#'               pig_channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin"),
#'               com_channels = c("FSC.HLin", "SSC.HLin"))
#'
#' ggpairsDens(flowfile = fin$flowframe_proportion,
#'                     channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin", "FSC.HLin"),
#'                     cluster_var = "Clusters")
#'
#'  }
#'
#' @return a ggplot object
#' @export ggpairsDens

ggpairsDens <- function(flowfile, channels = NULL, group = NULL, notToPlot = NULL,
                        ...) {

  if(is.null(notToPlot) & !is.null(channels)) {

    tsd1 <- expand.grid(channels,
                        channels,
                        stringsAsFactors = FALSE)

  } else if(!is.null(notToPlot) & is.null(channels)) {

    channels <- setdiff(flowCore::colnames(flowfile), notToPlot)
    tsd1 <- expand.grid(channels,
                        channels,
                        stringsAsFactors = FALSE)
  } else stop('channels and notToPlot cannot be empty')


  plts2 <- vector("list", nrow(tsd1))
  for(i in 1:nrow(tsd1)) {
    if(as.character(tsd1$Var1[i]) == as.character(tsd1$Var2[i])) {

      plts2[[i]] <- ggally_text(as.character(tsd1$Var1[i])) +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()
        )

    } else {

      if(is.null(group)) {

        plts2[[i]] <- ggplotDens(flowfile, as.character(tsd1[i, ]))

      } else {

        #plts2[[i]] <- ggplotDens2(flowfile, as.character(tsd1[i, ]),
        #                          group = group)
        chs <- as.character(tsd1[i, ])
        plotdata <- as.data.frame(flowfile@exprs[, c(chs, group)])
        plotdata$Clusters <- factor(plotdata$Clusters)
        plotdata_mean <- aggregate(flowfile@exprs[, channels],
                                   list(flowfile@exprs[, group]),
                                   median, na.rm = TRUE)

        plts2[[i]] <- ggplot(data = plotdata, aes_string(x = chs[1], y = chs[2], color = group)) +
          geom_point(pch = 19, size = 0.5, alpha = 0.3) +
          theme_minimal() +
          scale_color_brewer(type = "qual", palette = "Set1") +
          theme(axis.line.x = element_line( size = 0.2),
                axis.line.y = element_line( size = 0.2),
                plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
                axis.title = element_text(hjust = 0.5, size = 25, face = "bold"),
                axis.text = element_text(hjust = 0.5, size = 7, face = "bold")) +
          labs(x = chs[1], y = chs[2]) +
          geom_text(data = plotdata_mean,
                    aes_string(x = chs[1], y = chs[2]),
                    label = paste("Cluster", plotdata_mean[, 1], sep = "-"),
                    size = 6,
                    color = "grey0",
                    inherit.aes = FALSE)

      }


    }

  }

  ggmatrix(plts2,
           nrow = nrow(tsd1) / length(channels),
           ncol = nrow(tsd1) / length(channels),
           xAxisLabels = unique(tsd1$Var1),
           yAxisLabels = unique(tsd1$Var1),
           title = identifier(flowfile)
  ) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"))

}

#' plot method for cyanoFilter objects
#'
#' @param x object of class cyanoFilter to plotted
#' @param ch channels to be plotted. This must be supplied if x is a flowFrame.
#' @param ... other arguments. Not used at the moment
#'
#' @export
plot <- function(x, ch = NULL, ...) {
UseMethod('plot')
}

#' @export
plot.cyanoFilter <- function(x, ch = NULL, ...) {

  #chs <- match.arg(ch)

  if(class(x)[2] == 'phytoFilter' | class(x)[2] == 'Debris' | class(x)[2] == 'marginEvents') {

    return(x$plot)

  } else {
    stop('Error: incorrect input')
  } #else if(class(x) == 'flowFrame' & !is.null(ch) & length(ch) > 2) {

  #return(ggpairsDens(x, channels = ch))

  #} else if(class(x) == 'flowFrame' & !is.null(ch) & length(ch) == 2) {

  # return(ggplotDens(x, channels = ch))

  #}
}

#' function to assign particles to cluster based on a supplied classifier. It is useful for the EM clustering approach.
#'
#' @param flowfile flowframe to be plotted
#' @param classifier minimum desired probability
#'
#'@export cluster_assign

cluster_assign <- function(flowfile, classifier) {

  ddata <- flowCore::exprs(flowfile)

  ddata2 <- ddata[, stringr::str_detect(colnames(ddata), "Prob") ]

  cluster_assign <- apply(ddata2, 1, function(x) {

    rest <- which(x >= classifier & x == max(x))
    # maximum = not sure in other words NA
    frest <- ifelse(length(rest) == 0, 0, rest)

    return(frest)

  })

  return(cluster_assign)


}


#' plots the expression matrix of a flowframe analysed with celldebris_emclustering.
#'
#' @param flowfile flowframe to be plotted
#' @param channels channels used in gating
#' @param mus matrix of means obtained from celldebris_emclustering
#' @param tau vector of cluster weights obtained from celldebris_emclustering
#' @param classifier cells will be assigned to a cluster if belongs to that cluster
#'                    by at least this probability. Only for plotting purposes.
#'
#' @importFrom grDevices chull densCols
#' @importFrom graphics plot polygon
#' @export cluster_plot


cluster_plot <- function(flowfile, channels,
                         mus = NULL, tau = NULL,
                         classifier) {

  if(is.null(mus) | is.null(tau)) stop("supply the matrix of mean and
                                        vector of percentages obtained
                                       from the celldebris_emclustering function")

  ddata <- flowCore::exprs(flowfile)

  ddata2 <- ddata[, stringr::str_detect(colnames(ddata), "Prob") ]

  color_code <- apply(ddata2, 1, function(x) {

    rest <- which(x >= classifier & x == max(x))
    # maximum = not sure in other words NA
    frest <- ifelse(length(rest) == 0, 0, rest)

    return(frest)

  })

  #plotting
  color_seq <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "seq", ]
  panel.xy <- function(x, y, ...) {

    for(i in unique(color_code)) {
      if(i != 0) {

        cols.pal <- RColorBrewer::brewer.pal(n = 9, row.names(color_seq)[i])
        col.palette <- colorRampPalette(c("white", cols.pal))
        pal <- densCols(x[color_code == i], y[color_code == i],
                        colramp =  col.palette)
        points(x[color_code == i], y[color_code == i],
               col = pal, ...)

        #drawign cluster boundaries
        if(i == which(tau == max(tau)) ) {
          plot_data <- cbind(x[color_code == i], y[color_code == i])
          convhull <- chull(plot_data)
          polygon(plot_data[convhull, ], lty = 4, lwd = 2,
                  border = "red"
          ) }


      } else {

        cols.pal <- RColorBrewer::brewer.pal(n = 9, "Greys")
        col.palette <- colorRampPalette(c("white", cols.pal))
        pal <- densCols(x[color_code == i], y[color_code == i],
                        colramp =  col.palette)
        points(x[color_code == i], y[color_code == i],
               col = pal, ...)

      }

    }

  }

  pairs(ddata[, channels], pch = ".",
        panel = panel.xy, gap = 0.2,
        main = flowCore::identifier(flowfile)
  )

  msg <- paste("points are assigned to a cluster if they have at least",
               classifier, "of belonging to that cluster.
                 The boundary of the largest cluster is drawn but the
                 colors differentiate each cluster", sep = " ")
  message(msg)


}
