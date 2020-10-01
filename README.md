
# cyanoFilter

<!---[![Travis-CI Build Status](https://travis-ci.org/fomotis/cyanoFilter.svg?branch=master)](https://travis-ci.org/fomotis/cyanoFilter) -->

# Motivation and Background

Phytoplankton are primary producers responsible for about 50% of global
primary production through photosynthesis. Studies seeking a better
understanding of the ecology of phytoplankton have relied on flow
cytometry (FCM) to measure phytoplankton abundance and traits. FCM is a
technique which involves the suspension of cells or particles within a
fluid stream which is made to pass through one or more laser beams. A
crucial step in FCM application is to separate signal from noise, a
process termed gating. This process can be done manually, termed manual
gating, or with the use of automated algorithms in software packages,
termed automated gating. Gating is typically done manually or using
model-based tools such as *flowClust* and *flowEMMi* or machine
learning. While manual gating is often not fully reproducible,
model-based tools and machine learning require tuning of global
parameters that are not related to the biological properties of the
measured cells and cannot readily integrate experts’ knowledge.
*cyanoFilter* for semi-automated gating of phytoplankton FCM data. The
package uses pigment and complexity information to identify each
phytoplankton population contained in a sample. Aside from identifying
phytoplankton populations, *cyanoFilter* can also assist in the
identification of previously unknown light channels useful for
differentiating phytoplankton cells.

# Installation and Dependencies

Run the `code` below to install the package and all its dependencies.

``` r
remotes::install_github("https://github.com/fomotis/cyanoFilter")
```

All dependencies both on **CRAN** and **bioconductor** should be
installed when you install the package itself. However, do install the
following needed **bioconductor** packages should you run into errors
while attempting to use the functions in this package.

``` r
install.packages("BiocManager")
library(BiocManager)
install(c("Biobase", "flowCore", "flowDensity", 'flowViz))
```

# Usage

The package comes with 2 internal datasets that we use for demonstrating
the usage of the functions contained in the package. The **meta data**
file contains *BS4* and *BS5* samples measured with a guava easyCyte HT
series at 3 dilution levels (*2000*, *10000* and *20000*) each. The
**FCS** file contains the flow cytometer channel measurements for one of
these sample.

# Monoculture experiment

## Meta File Preprocessing

### The Good Measurements

The **goodfcs()** is deigned to check the \(cell/\mu\)L of the meta file
(normally csv) obtained from the flow cytometer and decide if the
measurements in the FCS file can be trusted. This function is especially
useful for flow cytometers that are not equipped to perform automated
dilution.

``` r
library(flowCore)
library(cyanoFilter)
#internally contained datafile in cyanoFilter
metadata <- system.file("extdata", "2019-03-25_Rstarted.csv", 
  package = "cyanoFilter", 
  mustWork = TRUE)
metafile <- read.csv(metadata, skip = 7, stringsAsFactors = FALSE, 
  check.names = TRUE)
#columns containing dilution, $\mu l$ and id information
metafile <- metafile[, c(1:3, 6:8)] 
knitr::kable(metafile) 
```

| Sample.Number | Sample.ID  | Number.of.Events | Dilution.Factor | Original.Volume |   Cells.L |
| ------------: | :--------- | ---------------: | --------------: | --------------: | --------: |
|             1 | BS4\_20000 |             6918 |           20000 |              10 |  62.02270 |
|             2 | BS4\_10000 |             6591 |           10000 |              10 | 116.76311 |
|             3 | BS4\_2000  |             6508 |            2000 |              10 | 517.90008 |
|             4 | BS5\_20000 |             5976 |           20000 |              10 |  48.31036 |
|             5 | BS5\_10000 |             5844 |           10000 |              10 |  90.51666 |
|             6 | BS5\_2000  |             5829 |            2000 |              10 | 400.72498 |

Each row in the csv file corresponds to a measurement from two types of
cyanobacteria cells carried out at one of three dilution levels. The
columns contain information about the dilution level, the number of
cells per micro-litre (\(cell/\mu l\)), number of particles measured and
a unique identification code for each measurement. The *Sample.ID*
column is structured in the format cyanobacteria\_dilution. We extract
the cyanobacteria part of this column into a new column and also rename
the \(cell/\mu l\) column with the following code:

``` r
#extract the part of the Sample.ID that corresponds to BS4 or BS5
metafile <- metafile %>% dplyr::mutate(Sample.ID2 = 
                                         stringr::str_extract(metafile$Sample.ID, "BS*[4-5]")
                                       )
#clean up the Cells.muL column
names(metafile)[which(stringr::str_detect(names(metafile), "Cells."))] <- "CellspML"
```

### Good Measurements

To determine the appropriate data file to read from a FCM datafile, the
desired minimum, maximum and column containing the \(cell\mu l\) values
are supplied to the **goodfcs()** function. The code below demonstrates
the use of this function for a situation where the desired minimum and
maximum for \(cell/\mu l\) is 50 and 1000 respectively.

``` r
metafile <- metafile %>% mutate(Status = cyanoFilter::goodfcs(metafile = metafile, col_cpml = "CellspML", 
                                        mxd_cellpML = 1000, mnd_cellpML = 50)
                                )
knitr::kable(metafile)
```

| Sample.Number | Sample.ID  | Number.of.Events | Dilution.Factor | Original.Volume |  CellspML | Sample.ID2 | Status |
| ------------: | :--------- | ---------------: | --------------: | --------------: | --------: | :--------- | :----- |
|             1 | BS4\_20000 |             6918 |           20000 |              10 |  62.02270 | BS4        | good   |
|             2 | BS4\_10000 |             6591 |           10000 |              10 | 116.76311 | BS4        | good   |
|             3 | BS4\_2000  |             6508 |            2000 |              10 | 517.90008 | BS4        | good   |
|             4 | BS5\_20000 |             5976 |           20000 |              10 |  48.31036 | BS5        | bad    |
|             5 | BS5\_10000 |             5844 |           10000 |              10 |  90.51666 | BS5        | good   |
|             6 | BS5\_2000  |             5829 |            2000 |              10 | 400.72498 | BS5        | good   |

The function adds an extra column, *Status*, with entries *good* or
*bad* to the metafile. Rows containing \(cell/\mu l\) values outside the
desired minimum and maximum are labelled *bad*. Note that the *Status*
column for the fourth row is labelled *bad*, because it has a
\(cell/\mu l\) value outside the desired range.

### Files to Retain

Although any of the files labelled good can be read from the FCM file,
the **retain()** function can help select either the file with the
highest \(cell/\mu l\) or that with the smallest \(cell/\mu l\) value.
To do this, one supplies the function with the status column,
\(cell/\mu l\) column and the desired decision. The code below
demonstrates this action for a case where we want to select the file
with the maximum \(cell/\mu l\) from the good measurements for each
unique sample ID.

``` r
broken <- metafile %>% group_by(Sample.ID2) %>% nest()
metafile$Retained <- unlist(map(broken$data, function(.x) {
  retain(meta_files = .x, make_decision = "maxi",
  Status = "Status",
  CellspML = "CellspML")
 })
)
knitr::kable(metafile)
```

| Sample.Number | Sample.ID  | Number.of.Events | Dilution.Factor | Original.Volume |  CellspML | Sample.ID2 | Status | Retained |
| ------------: | :--------- | ---------------: | --------------: | --------------: | --------: | :--------- | :----- | :------- |
|             1 | BS4\_20000 |             6918 |           20000 |              10 |  62.02270 | BS4        | good   | No\!     |
|             2 | BS4\_10000 |             6591 |           10000 |              10 | 116.76311 | BS4        | good   | No\!     |
|             3 | BS4\_2000  |             6508 |            2000 |              10 | 517.90008 | BS4        | good   | Retain   |
|             4 | BS5\_20000 |             5976 |           20000 |              10 |  48.31036 | BS5        | bad    | No\!     |
|             5 | BS5\_10000 |             5844 |           10000 |              10 |  90.51666 | BS5        | good   | No\!     |
|             6 | BS5\_2000  |             5829 |            2000 |              10 | 400.72498 | BS5        | good   | Retain   |

This function adds another column, *Retained*, to the metafile. The
third and sixth row in the metadata are with the highest \(cell/\mu l\)
values, thus one can proceed to read the fourth and sixth file from the
corresponding FCS file for *BS4* and *BS5* respectively. This implies
that we are reading in only two FCS files rather than the six measured
files.

## Flow Cytometer File Processing

To read **B4\_18\_1.fcs** file into **R**, we use the **read.FCS()**
function from the **flowCore** package. The *dataset* option enables the
specification of the precise file to be read. Since this datafile
contains one file only, we set this option to 1. If this option is set
to 2, it gives an error since **text.fcs** contains only one datafile.

``` r
flowfile_path <- system.file("extdata", "B4_18_1.fcs", package = "cyanoFilter",
  mustWork = TRUE)
flowfile <- read.FCS(flowfile_path, alter.names = TRUE,
  transformation = FALSE, emptyValue = FALSE,
  dataset = 1)
flowfile
> flowFrame object ' B4_18_1'
> with 8729 cells and 11 observables:
>            name                                desc range    minRange maxRange
> $P1    FSC.HLin          Forward Scatter (FSC-HLin) 1e+05    0.000000    99999
> $P2    SSC.HLin             Side Scatter (SSC-HLin) 1e+05  -34.479282    99999
> $P3  GRN.B.HLin   Green-B Fluorescence (GRN-B-HLin) 1e+05  -21.194536    99999
> $P4  YEL.B.HLin  Yellow-B Fluorescence (YEL-B-HLin) 1e+05  -10.327441    99999
> $P5  RED.B.HLin     Red-B Fluorescence (RED-B-HLin) 1e+05   -5.347203    99999
> $P6  NIR.B.HLin Near IR-B Fluorescence (NIR-B-HLin) 1e+05   -4.307983    99999
> $P7  RED.R.HLin     Red-R Fluorescence (RED-R-HLin) 1e+05  -25.490185    99999
> $P8  NIR.R.HLin Near IR-R Fluorescence (NIR-R-HLin) 1e+05  -16.020023    99999
> $P9    SSC.ALin        Side Scatter Area (SSC-ALin) 1e+05    0.000000    99999
> $P10      SSC.W          Side Scatter Width (SSC-W) 1e+05 -111.000000    99999
> $P11       TIME                                Time 1e+05    0.000000    99999
> 368 keywords are stored in the 'description' slot
```

The **R** object *flowfile* contains measurements about cells across 10
channels since the time channel does not contain any information about
the properties of the measured cells.

### Transformation and visualisation

To examine the need for transformation, a visual representation of the
information in the expression matrix is of great use. The
**ggpairsDens()** function produces a panel plot of all measured
channels. Each plot is also smoothed to show the cell density at every
part of the plot.

``` r
flowfile_nona <- nona(x = flowfile)
ggpairsDens(flowfile_nona, notToPlot = "TIME")
```

![**Panel plot for all channels measured in flowfile\_nona. A bivariate
kernel smoothed color density is used to indicate the cell
density.**](man/figures/README-remove_na-1.png)

We obtain Figure above by using the **ggpairsDens()** function after
removing all `NA` values from the expression matrix with the **nona()**
function. There is a version of the function, **pairs\_plot()** that
produces standard base scatter plots also smoothed to indicate cell
density.

``` r

flowfile_noneg <- noneg(x = flowfile_nona)
flowfile_logtrans <- lnTrans(x = flowfile_noneg, 
  notToTransform = c("SSC.W", "TIME"))
ggpairsDens(flowfile_logtrans, notToPlot = "TIME")
```

![Panel plot for log-transformed channels for flowfile\_logtrans. A
bivariate kernel smoothed color density is used to indicate the cell
density.](man/figures/README-logtrans-1.png)

The second figure is the result of performing a logarithmic
transformation in addition to the previous actions taken. The
logarithmic transformation appears satisfactory in this case, as it
allow a better examination of the information contained in each panel of
the figure. Moreover, the clusters are clearly visible in this figure
compared to the former figure. Other possible transformation (linear,
bi-exponential and arcsinh) can be pursued if the logarithm
transformation is not satisfactory. Functions for these transformations
are provided in the **flowCore** package.

## Gating

Flow cytometry outcomes can be divided into 3 and they are not entirely
mutually exclusive but this is not a problem as scientists are often
interested in a pre-defined outcome.

![Flow Cytometry
Outcomes](README_files/figure-gfm/flowcytometryOutcome.PNG)

  - Margin Events are particles too big to be measured
  - Doublets/Multiplets are cells with disproportionate Area, Height
    relationship
  - Singlets are the ‘normal cells’ but these could either be dead
    cells/particles (debris) or living cells (good cells).

The set of functions below identifies margin events and singlets.
Doublets are normally pre-filtered during the event acquiring phase when
running the flow cytometer.

The set of functions below identifies margin events and singlets.
Doublets are normally pre-filtered during the event

### Gating margin events

To remove margin events, the **cellmargin()** function takes the column
in the expression matrix corresponding to measurements about the width
of each cell. The code below demonstrates the removal of margin events
using the SSC.W column with the option to estimate the cut point between
the margin events and the good cells.

``` r
flowfile_marginout <- cellmargin(flowframe = flowfile_logtrans,
                                 Channel = 'SSC.W', type = 'estimate', 
                                 y_toplot = "FSC.HLin")
plot(flowfile_marginout)
```

![Smoothed Scatterplot of measured width (SSC.W) and height (FSC.HLin).
The red line is the estimated cut point by flowDensity, and every
particle below the red line has their width properly
measured.](man/figures/README-marginEvents-1.png)

``` r

summary(flowfile_marginout, 
       channels = c('FSC.HLin', 'SSC.HLin', 
                    'SSC.W'))
> $means
>     FSC.HLin     SSC.HLin        SSC.W 
>     4.981515     4.732556 37223.493869 
> 
> $medians
>     FSC.HLin     SSC.HLin        SSC.W 
>     4.780560     4.543738 36168.914062 
> 
> $variances
>            FSC.HLin     SSC.HLin        SSC.W
> FSC.HLin  0.8243962    0.7505259 8.770121e+01
> SSC.HLin  0.7505259    1.1430019 4.388930e+03
> SSC.W    87.7012064 4388.9298025 3.418691e+08
> 
> $marginCount
> [1] 3092
> 
> $nonMarginCount
> [1] 3831
```

*flowfile\_marginout* is an S3 object of class `cyanoFilter` with
**summary()** and **plot()** method. This object can be subsetted like a
normal list. Running **plot()** on *flowfile\_marginout* produces a plot
of the width channel against the channel supplied in *y\_toplot*. This
action returns the figure @ref(fig:marginEvents). flowfile\_marginout
contains the following sub-objects:

The function returns a figure (Figure @ref(fig:marginEvents)) in this
case) and a list containing:

  - *fullflowframe*, flowframe with indicator for margin and non-margin
    events in the expression matrix,
  - *reducedflowframe*, flowframe containing only non-margin events
  - *N\_margin*, number of margin events contained in the input
    flowframe
  - *N\_nonmargin*, number of non-margin events
  - *N\_particle*, number of particles in the input flowframe

Running **plot()** on *flowfile\_marginout* gives you the number of
margin and non-margin particles as well as descriptives on channels
supplied. These descriptives are computed on the flowfile after the
margin events have been removed.

### Gating Debris

To identify debris, we leverage on the presence of chlorophyll *a*

``` r

cells_nodebris <-  debris_nc(flowframe = flowfile_marginout$reducedflowframe, 
                             ch_chlorophyll = "RED.B.HLin", ch_p2 = "YEL.B.HLin",
                             ph = 0.05)
plot(cells_nodebris)
```

![Smoothed Scatterplot of measured chlorophyll *a* channel (RED.B.HLin)
and phycoerythrin channel (YEL.B.HLin). The red lines are the estimated
minimum intersection points between the detected
peaks.](man/figures/README-Debris-1.png)

### Gating cyanobacteria

We conceptualized the division of cells into clusters in two ways in
cyanoFilter and this is reflected in two main functions that perform the
clustering exercise; **phyto\_filter()** and
**celldebris\_emclustering()**. The **phyto\_filter()** function employs
the following algorithm to separate particles into different clusters;

1.  Search for peaks along the supplied pigment and cell complexity
    channels.
2.  Idneify the minimum intersection point between the peaks observed
    these channels.
3.  Divide particles into groups based on the minimum intersection
    points identified in 1 and label each group.
4.  Formulate all possible combinations of labels in step 2.
5.  Assign a new label to the combinations in 3.
6.  Retain clusters that make up a desired proportion of the total
    number of particles clustered.

<!-- end list -->

``` r

bs4_gate1 <- phyto_filter(flowfile = cells_nodebris$syn,
               pig_channels = c("RED.B.HLin", "YEL.B.HLin", "RED.R.HLin"),
               com_channels = c("FSC.HLin", "SSC.HLin"))

plot(bs4_gate1)
```

![Smoothed Scatterplot of all channels used in the gating
process.](man/figures/README-kdapproach-1.png)

The resulting object is a figure (Figure @ref(fig:kdapproach)) and a
list containing the following:

  - *reducedframe*, a flowFrame with all debris removed
  - *fullframe*, flowFrame with all measured particles and indicator for
    debris and cyanobacteria cells
  - *Cell\_count*, the number of BS4 cells counted
  - *Debris\_Count*, the number of debris particles.

# Gating a bi-culture experiment

# License

This is a free to use package for anyone who has the need. However,
users must adhere to the licensing agreement of some dependencies that
require that their packages be used only for educational and research
purposes.
