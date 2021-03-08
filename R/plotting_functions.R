
#' produces a scatter plot of the expression matrix of a flowframe. Note that, 
#' it takes some time to display the plot.
#'
#' @param x flowframe to be plotted
#' @param notToPlot column in expression matrix not to be plotted
#' @param ... other arguments. Not used at the moment
#' @return a plot object
#' @import ggplot2 GGally flowCore Biobase
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics plot.default pairs.default abline points panel.smooth 
#'             pairs smoothScatter text
#' @examples
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#'               package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::noNA(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noNeg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#'                       c('SSC.W', 'TIME'))
#' pairsPlot(flowfile_logtrans,
#'            notToPlot = c("TIME", "SSC.W",
#'            "SSC.HLin", "NIR.R.HLin", 
#'            "FSC.HLin"))  
#'            
#' @importFrom stats var cov quantile runif
#' @importFrom grDevices chull densCols
#' @export pairsPlot

pairsPlot <- function(x, notToPlot = c("TIME"), ...) {

  toplot <- setdiff(flowCore::colnames(x), notToPlot)
  col.palette <- colorRampPalette(c("white", "blue", "cyan", "green", "orange", 
                                    "red"),
                                  space = "Lab")

  if(!is.DebrisFilter(x) | !is.MarginEvents(x) | !is.PhytoFilter(x)) {


    pairs.default(flowCore::exprs(x)[, toplot], pch = ".",
                  panel = function(...) smoothScatter(..., nrpoints = 0,
                                                      colramp = col.palette,
                                                      add = TRUE), gap = 0.2,
                  main = flowCore::identifier(x))
  } else if(is.DebrisFilter(x) | is.MarginEvents(x) | is.PhytoFilter(x)) {


    pairs.default(flowCore::exprs(fullFlowframe(x))[, toplot], pch = ".",
                  panel = function(...) smoothScatter(..., nrpoints = 0,
                                                      colramp = col.palette,
                                                      add = TRUE), gap = 0.2,
                  main = flowCore::identifier(x))

  } else stop('Error: object supplied not supported')

}



#' plots two channels of a flowframe.
#'
#' @param flowfile flowframe to be plotted
#' @param channels a character vector of length 2, must contain channel names 
#'                 in the flowfile.
#' @param ... not used at the moment
#'
#' @examples
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#'               package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::noNA(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noNeg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#'                       c('SSC.W', 'TIME'))
#' ggplotDens(flowfile_logtrans,
#'            channels = c("FSC.HLin", "RED.R.HLin"))
#' @return a ggplot object
#'
#' @export ggplotDens

ggplotDens <- function(flowfile, channels, ...) {

  plotdata <- as.data.frame(flowfile@exprs[, channels])
  d <- densCols(plotdata, colramp = colorRampPalette(c("white", "blue", "cyan",
                                                       "green", "orange", 
                                                       "red"),
                                                     space = "Lab"))

  plotdata <- cbind(plotdata, d = d)
  ggplot(data = plotdata, aes(x = get(channels[1], plotdata), y = 
                                get(channels[2], plotdata), color = d)) +
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


#' plots two channels of a flowframe with different colors for clusters 
#'  identified.
#'
#' @param flowfile flowframe to be plotted
#' @param channels a character vector of length 2, must contain channel names 
#'                 in the flowfile.
#' @param group cluster groups. must be equal to the number of particles in the 
#'             flow cytometer.
#' @param ... not used at the moment
#' 
#' @examples 
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#' package = "cyanoFilter",
#'               mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, 
#'                                emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::noNA(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noNeg(x = flowfile_nona)
#' flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
#'                       c('SSC.W', 'TIME'))
#' cells_nonmargin <- cellMargin(flowframe = flowfile_logtrans, 
#'                               Channel = 'SSC.W',
#'            type = 'estimate', y_toplot = "FSC.HLin")
#' cells_nodebris <- debrisNc(flowframe = reducedFlowframe(cells_nonmargin),
#'                             ch_chlorophyll = "RED.B.HLin",
#'                             ch_p2 = "YEL.B.HLin",
#'                             ph = 0.05)
#' cct <- phytoFilter(flowfile = reducedFlowframe(cells_nodebris),
#'               pig_channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin"),
#'               com_channels = c("FSC.HLin", "SSC.HLin"))
#' ggplotDens2(reducedFlowframe(cct), 
#' c("RED.B.HLin", "YEL.B.HLin"),
#' group = "Clusters")
#' @return a ggplot object
#'
#' @export ggplotDens2

ggplotDens2 <- function(flowfile, channels, group, ...) {

  plotdata <- as.data.frame(flowfile@exprs[, c(channels, group)])

  plotdata_mean <- aggregate(flowfile@exprs[, channels],
                             list(flowfile@exprs[, group]),
                             mean, na.rm = TRUE)
  names(plotdata_mean)[1] <- group
  Group <- unique(flowfile@exprs[, group])
  hulls <- lapply(Group, function(x) {

    pd <- plotdata[plotdata[, group] == x, ]
    pd[chull(pd[, seq_len(2)]), ]

  })

  hulldata <- do.call(rbind.data.frame, hulls)

  d <- densCols(plotdata, colramp = colorRampPalette(c("white", "blue", "cyan",
                                                       "green", "orange", 
                                                       "red"),
                                                     space = "Lab"))
  plotdata <- cbind(plotdata, d = d)
  ggplot(data = plotdata, aes(x = get(channels[1], plotdata), 
                              y = get(channels[2], plotdata), 
                              color = d)) +
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
                                        label = paste("Cluster", 
                                                      plotdata_mean[, 1], 
                                                      sep = "-")
    ),
    size = 5,
    color = get(group, plotdata_mean),
    inherit.aes = FALSE)

}


#' produces a scatter plot of the expression matrix of the flowframe. 
#' If a cluster variable is given,
#' it assigns different colors to the clusters.
#'
#' @param flowfile flowframe to be plotted
#' @param channels a character vector of length 2 or more. It must contain 
#'                channel names in the flowfile.
#' @param group cluster groups. It must be equal to the number of particles 
#'              in the flowfile. If group is null cluster boundaries are not 
#'              drawn.
#' @param notToPlot columns not to plot. This is especially useful for for 
#'                  plotting all columns in a
#' @param ... not used at the moment
#' 
#'  @return a ggplot object
#'
#' @examples
#'
#'  # example without clustering
#' flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
#'                  package = "cyanoFilter",
#'                  mustWork = TRUE)
#' flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
#'                                transformation = FALSE, emptyValue = FALSE,
#'                                dataset = 1) 
#' flowfile_nona <- cyanoFilter::noNA(x = flowfile)
#' flowfile_noneg <- cyanoFilter::noNeg(x = flowfile_nona)
#' flowfile_logtrans <- lnTrans(x = flowfile_noneg, c('SSC.W', 'TIME'))
#' ggpairsDens(flowfile = flowfile_logtrans,
#'             channels = c("FSC.HLin", "RED.R.HLin", "RED.B.HLin", 
#'             "NIR.R.HLin"))
#' @return a ggplot object
#' @export ggpairsDens

ggpairsDens <- function(flowfile, channels = NULL, group = NULL, 
                        notToPlot = NULL,
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
  for(i in seq_len(nrow(tsd1))) {
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

        plts2[[i]] <- ggplot(data = plotdata, aes_string(x = chs[1], y = chs[2],
                                                         color = group)) +
          geom_point(pch = 19, size = 0.5, alpha = 0.3) +
          theme_minimal() +
          scale_color_brewer(type = "qual", palette = "Set1") +
          theme(axis.line.x = element_line( size = 0.2),
                axis.line.y = element_line( size = 0.2),
                plot.title = element_text(hjust = 0.5, size = 25, 
                                          face = "bold"),
                axis.title = element_text(hjust = 0.5, size = 25, 
                                          face = "bold"),
                axis.text = element_text(hjust = 0.5, size = 7, 
                                         face = "bold")) +
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
