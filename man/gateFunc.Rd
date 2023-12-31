% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accTest.R
\name{gateFunc}
\alias{gateFunc}
\title{tests the accuracy of several automated gating functions on monoculture
flow cytometry experiments.}
\usage{
gateFunc(
  flowfile,
  sfts = c("phytoFilter", "flowClust", "cytometree"),
  channels,
  funargs_list
)
}
\arguments{
\item{flowfile}{flowSet with each flowFrame being a phytoplankton 
monoculture FCM experiment}

\item{sfts}{character vector of gating function to test.}

\item{channels}{channels to be used for gating}

\item{funargs_list}{additional options for the chosen gating function}
}
\value{
a flowFrame with cluster indicator generated by the software used
        added to the expression matrix.
}
\description{
This function gates all flowFrames in the supplied flowSet
             to attach cluster labels. Then it mixes up the flowSet into
             one giant flowFrame and re-gates this to attach another label.
             These labels are used to examine if the gating algorithms can
             reproduce the earlier clusters before the mixing.
}
\examples{
 flowfile_path <- system.file("extdata", "B4_18_1.fcs", 
package = "cyanoFilter",
              mustWork = TRUE)
flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
                               transformation = FALSE, 
                               emptyValue = FALSE,
                               dataset = 1) 
flowfile_nona <- cyanoFilter::noNA(x = flowfile)
flowfile_noneg <- cyanoFilter::noNeg(x = flowfile_nona)
flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, 
                      c('SSC.W', 'TIME'))
cells_nonmargin <- cellMargin(flowframe = flowfile_logtrans, 
                              Channel = 'SSC.W',
           type = 'estimate', y_toplot = "FSC.HLin")
cells_nodebris <- debrisNc(flowframe = reducedFlowframe(cells_nonmargin),
                            ch_chlorophyll = "RED.B.HLin",
                            ch_p2 = "YEL.B.HLin",
                            ph = 0.05)
#phytoFilter specification
gateFunc(flowfile = reducedFlowframe(cells_nodebris),
              channels = c("RED.B.HLin", "YEL.B.HLin", 
              "RED.R.HLin", "FSC.HLin", "SSC.HLin"),
              sfts = "phytoFilter", 
              list(ph = 0.1, proportion = 0.90)
              )
#flowClust specification
gateFunc(flowfile = reducedFlowframe(cells_nodebris),
              channels = c("RED.B.HLin", "YEL.B.HLin", 
              "RED.R.HLin", "FSC.HLin", "SSC.HLin"),
              sfts = "flowClust", 
              list(K = 1:4, B = 100)
              )
#cytometree specification
gateFunc(flowfile = reducedFlowframe(cells_nodebris),
              channels = c("RED.B.HLin", "YEL.B.HLin", 
              "RED.R.HLin", "FSC.HLin", "SSC.HLin"),
              sfts = "cytometree", 
              list(minleaf = 1, t = 0.10)
              )

}
