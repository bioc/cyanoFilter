
# cyanoFilter

[![Travis-CI Build
Status](https://travis-ci.org/fomotis/cyanoFilter.svg?branch=master)](https://travis-ci.org/fomotis/cyanoFilter)

cyaoFilter is a package designed to identify, assign indicators and/or
filter out synechoccus type cyanobacteria from a water sample examined
with flowcytometry.

# Installation and Dependencies

Run the `code` below to install the package and all its dependencies.

``` r
install.packages("cyanoFilter")
```

All dependencies both on **CRAN** and **bioconductor** should be
installed when you install the package itself. However, do install the
following needed **bioconductor** packages should you run into errors
while attempting to use the functions in this package.

``` r
install.packages("BiocManager")
library(BiocManager)
install(c("Biobase", "flowCore", "flowDensity"))
```

# Motivation and Background

Flow cytometry is a well-known technique for identifying cell
populations in fluids. It is largely applied in biological and medical
sciences for cell sorting, counting, biomarker detections and protein
engineering. Identifying cell populations in flow cytometry data often
rely on manual gating, a subjective and generally irreproducible method
based on expert knowledge. To address this issue, two filtering
frameworks were developed in **R**, identifying and filtering out two
strains of Synechococcus type cyanobacteria (*BS4* and *BS5*) from flow
cytometry data.

# Usage

The package comes with 2 internal datasets that is used for
demonstrating the usage of the functions contained in the package. The
**meta data** file contains *BS4* and *BS5* samples examined with a
GUAVAVA flow cytometer at 3 dilution levels (*2000*, *10000* and
*20000*) each. The **FCS** file contains the flow cytometer channel
measurements for one of these sample.

## Meta File Preprocessing

### The Good Measurements

`goodfcs()` is deigned to check the cells/\(\mu\)L of the meta file
(normally csv) obtained from the flow cytometer and decide if the
measurements in the FCS file can be trusted.

``` r
#internally contained datafile in cyanoFilter
metadata <- system.file("extdata", "2019-03-25_Rstarted.csv", package = "cyanoFilter",
               mustWork = TRUE)
metafile <- read.csv(metadata, skip = 7, stringsAsFactors = FALSE, check.names = TRUE)
metafile <- metafile[, 1:65] #first 65 columns contains useful information
#extract the part of the Sample.ID that corresponds to BS4 or BS5
metafile$Sample.ID2 <- stringr::str_extract(metafile$Sample.ID, "BS*[4-5]")
#clean up the Cells.muL column
names(metafile)[which(stringr::str_detect(names(metafile), "Cells."))] <- "CellspML"
metafile$Status <- cyanoFilter::goodfcs(metafile = metafile, col_cpml = "CellspML", 
                                        mxd_cellpML = 1000, mnd_cellpML = 50)
#should work fine with tidyverse setup
metafile <- metafile %>% mutate(Status = cyanoFilter::goodfcs(metafile = metafile, col_cpml = 
                                                                "CellspML", mxd_cellpML = 1000,
                                                              mnd_cellpML = 50))
#the whole metadata file
knitr::kable(metafile)
```

| Sample.Number | Sample.ID  | Number.of.Events | Termination.Count | Count.Gate | Dilution.Factor | Original.Volume |  CellspML | Total.Volume | Acquisition.Time..s. | Date.of.Acquisition | Time.of.Acquisition | FSC.Gain | SSC |    GRN.B |    YEL.B | RED.B | NIR.B | RED.R | NIR.R | SSC.1 | GRN.B.1 | YEL.B.1 | RED.B.1 | NIR.B.1 | RED.R.1 | NIR.R.1 | YEL.B….GRN.B | RED.B….GRN.B | NIR.B….GRN.B | RED.R….GRN.B | NIR.R….GRN.B | GRN.B….YEL.B | RED.B….YEL.B | NIR.B….YEL.B | RED.R….YEL.B | NIR.R….YEL.B | GRN.B….RED.B | YEL.B….RED.B | NIR.B….RED.B | RED.R….RED.B | NIR.R….RED.B | GRN.B….NIR.B | YEL.B….NIR.B | RED.B….NIR.B | RED.R….NIR.B | NIR.R….NIR.B | GRN.B….RED.R | YEL.B….RED.R | RED.B….RED.R | NIR.B….RED.R | NIR.R….RED.R | GRN.B….NIR.R | YEL.B….NIR.R | RED.B….NIR.R | NIR.B….NIR.R | RED.R….NIR.R | Threshold.Parameter | Threshold.Value | Flow.Rate | High.Concentration.Warning.Trigger | X..of.Errors…Warnings | User.Login.Name | User.Full.Name | WorkList                                                           | Sample.ID2 | Status |
| ------------: | :--------- | ---------------: | ----------------: | :--------- | --------------: | --------------: | --------: | -----------: | -------------------: | :------------------ | :------------------ | -------: | --: | -------: | -------: | ----: | ----: | ----: | ----: | :---- | :------ | :------ | :------ | :------ | :------ | :------ | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :------------------ | --------------: | --------: | ---------------------------------: | --------------------: | :-------------- | :------------- | :----------------------------------------------------------------- | :--------- | :----- |
|             1 | BS4\_20000 |             6918 |              5000 | Cyanos     |           20000 |              10 |  62.02270 |    111.53980 |            189.05051 | 25-MAR-2019         | 12:42:13            | 10.37472 |   1 | 9.934862 | 3.084422 |     1 |     1 |     1 |     1 | Yes   | Yes     | Yes     | Yes     | Yes     | Yes     | Yes     | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | FSC                 |              50 |      0.59 |                                500 |                     1 | DESKTOP-9MNKV60 | Unknown        | C:/Users/User/Documents/Backup data/2019-03-25\_at\_12-37-08pm.xml | BS4        | good   |
|             2 | BS4\_10000 |             6591 |              5000 | Cyanos     |           10000 |              10 | 116.76311 |     56.44762 |             95.67394 | 25-MAR-2019         | 12:44:25            | 10.37472 |   1 | 9.934862 | 3.084422 |     1 |     1 |     1 |     1 | Yes   | Yes     | Yes     | Yes     | Yes     | Yes     | Yes     | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | FSC                 |              50 |      0.59 |                                500 |                     1 | DESKTOP-9MNKV60 | Unknown        | C:/Users/User/Documents/Backup data/2019-03-25\_at\_12-37-08pm.xml | BS4        | good   |
|             3 | BS4\_2000  |             6508 |              5000 | Cyanos     |            2000 |              10 | 517.90008 |     12.56613 |             21.29853 | 25-MAR-2019         | 12:45:19            | 10.37472 |   1 | 9.934862 | 3.084422 |     1 |     1 |     1 |     1 | Yes   | Yes     | Yes     | Yes     | Yes     | Yes     | Yes     | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | FSC                 |              50 |      0.59 |                                500 |                     2 | DESKTOP-9MNKV60 | Unknown        | C:/Users/User/Documents/Backup data/2019-03-25\_at\_12-37-08pm.xml | BS4        | good   |
|             4 | BS5\_20000 |             5976 |              5000 | Cyanos     |           20000 |              10 |  48.31036 |    123.70018 |            209.66132 | 25-MAR-2019         | 12:49:22            | 10.37472 |   1 | 9.934862 | 3.084422 |     1 |     1 |     1 |     1 | Yes   | Yes     | Yes     | Yes     | Yes     | Yes     | Yes     | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | FSC                 |              50 |      0.59 |                                500 |                     2 | DESKTOP-9MNKV60 | Unknown        | C:/Users/User/Documents/Backup data/2019-03-25\_at\_12-37-08pm.xml | BS5        | bad    |
|             5 | BS5\_10000 |             5844 |              5000 | Cyanos     |           10000 |              10 |  90.51666 |     64.56270 |            109.42831 | 25-MAR-2019         | 12:51:49            | 10.37472 |   1 | 9.934862 | 3.084422 |     1 |     1 |     1 |     1 | Yes   | Yes     | Yes     | Yes     | Yes     | Yes     | Yes     | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | FSC                 |              50 |      0.59 |                                500 |                     1 | DESKTOP-9MNKV60 | Unknown        | C:/Users/User/Documents/Backup data/2019-03-25\_at\_12-37-08pm.xml | BS5        | good   |
|             6 | BS5\_2000  |             5829 |              5000 | Cyanos     |            2000 |              10 | 400.72498 |     14.54614 |             24.65447 | 25-MAR-2019         | 12:52:47            | 10.37472 |   1 | 9.934862 | 3.084422 |     1 |     1 |     1 |     1 | Yes   | Yes     | Yes     | Yes     | Yes     | Yes     | Yes     | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | FSC                 |              50 |      0.59 |                                500 |                     1 | DESKTOP-9MNKV60 | Unknown        | C:/Users/User/Documents/Backup data/2019-03-25\_at\_12-37-08pm.xml | BS5        | good   |

``` r
#some interesting columns and the newly added column
knitr::kable(metafile %>% 
              dplyr::select(Sample.ID, Sample.ID2, Number.of.Events, Dilution.Factor, CellspML, Status) )
```

| Sample.ID      | Sample.ID2      |    Number.of.Events |   Dilution.Factor |     CellspML | Status   |
| :------------- | :-------------- | ------------------: | ----------------: | -----------: | :------- |
| BS4\_20000     | BS4             |                6918 |             20000 |     62.02270 | good     |
| BS4\_10000     | BS4             |                6591 |             10000 |    116.76311 | good     |
| BS4\_2000      | BS4             |                6508 |              2000 |    517.90008 | good     |
| BS5\_20000     | BS5             |                5976 |             20000 |     48.31036 | bad      |
| BS5\_10000     | BS5             |                5844 |             10000 |     90.51666 | good     |
| BS5\_2000      | BS5             |                5829 |              2000 |    400.72498 | good     |
| The \*\*Status | \*\* columns in | dicates if the file | at the current di | lution level | is good. |

### Files to Retain

Generally, reading **FCS** files with the `read.FCS()` function from the
`flowCore` package takes time, hence it can save you considerable amount
of time to read only the good files. An **FCS** normally contains data
from many experiment or same experiment measured at different dilution
levels. Some of these would be determined bad by the `goodfcs()`
function and should be avoided. The `retain()` function is especially
designed for this.

Since each sample, i.e. *BS4* and *BS5* were measured at 3 dilution
levels it means that the 3 rows in *metafile* containing the
measurements of *BS4* are to be examined together by the `retain()`
function and the same should be done for the *BS5* dilution levels. This
can easily be achieved by using some `tidyverse` function to break the
*metafile* into two based on *Sample.ID2* and apply `retain()` on the
broken dataset.

``` r
#break csv file into groups (2 in this case) based on sample ID2
broken <- metafile %>% group_by(Sample.ID2) %>% nest() 
#this is how broken looks like
broken
> # A tibble: 2 x 2
>   Sample.ID2 data             
>   <chr>      <list>           
> 1 BS4        <tibble [3 x 66]>
> 2 BS5        <tibble [3 x 66]>

# Let's apply the function
metafile$Retained <- unlist(map(broken$data, function(.x) {
    cyanoFilter::retain(meta_files = .x, make_decision = "maxi",
                      Status = "Status", 
                      CellspML = "CellspML")
  })
)  
#the whole metadata file
knitr::kable(metafile)
```

| Sample.Number | Sample.ID  | Number.of.Events | Termination.Count | Count.Gate | Dilution.Factor | Original.Volume |  CellspML | Total.Volume | Acquisition.Time..s. | Date.of.Acquisition | Time.of.Acquisition | FSC.Gain | SSC |    GRN.B |    YEL.B | RED.B | NIR.B | RED.R | NIR.R | SSC.1 | GRN.B.1 | YEL.B.1 | RED.B.1 | NIR.B.1 | RED.R.1 | NIR.R.1 | YEL.B….GRN.B | RED.B….GRN.B | NIR.B….GRN.B | RED.R….GRN.B | NIR.R….GRN.B | GRN.B….YEL.B | RED.B….YEL.B | NIR.B….YEL.B | RED.R….YEL.B | NIR.R….YEL.B | GRN.B….RED.B | YEL.B….RED.B | NIR.B….RED.B | RED.R….RED.B | NIR.R….RED.B | GRN.B….NIR.B | YEL.B….NIR.B | RED.B….NIR.B | RED.R….NIR.B | NIR.R….NIR.B | GRN.B….RED.R | YEL.B….RED.R | RED.B….RED.R | NIR.B….RED.R | NIR.R….RED.R | GRN.B….NIR.R | YEL.B….NIR.R | RED.B….NIR.R | NIR.B….NIR.R | RED.R….NIR.R | Threshold.Parameter | Threshold.Value | Flow.Rate | High.Concentration.Warning.Trigger | X..of.Errors…Warnings | User.Login.Name | User.Full.Name | WorkList                                                           | Sample.ID2 | Status | Retained |
| ------------: | :--------- | ---------------: | ----------------: | :--------- | --------------: | --------------: | --------: | -----------: | -------------------: | :------------------ | :------------------ | -------: | --: | -------: | -------: | ----: | ----: | ----: | ----: | :---- | :------ | :------ | :------ | :------ | :------ | :------ | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :------------------ | --------------: | --------: | ---------------------------------: | --------------------: | :-------------- | :------------- | :----------------------------------------------------------------- | :--------- | :----- | :------- |
|             1 | BS4\_20000 |             6918 |              5000 | Cyanos     |           20000 |              10 |  62.02270 |    111.53980 |            189.05051 | 25-MAR-2019         | 12:42:13            | 10.37472 |   1 | 9.934862 | 3.084422 |     1 |     1 |     1 |     1 | Yes   | Yes     | Yes     | Yes     | Yes     | Yes     | Yes     | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | FSC                 |              50 |      0.59 |                                500 |                     1 | DESKTOP-9MNKV60 | Unknown        | C:/Users/User/Documents/Backup data/2019-03-25\_at\_12-37-08pm.xml | BS4        | good   | No\!     |
|             2 | BS4\_10000 |             6591 |              5000 | Cyanos     |           10000 |              10 | 116.76311 |     56.44762 |             95.67394 | 25-MAR-2019         | 12:44:25            | 10.37472 |   1 | 9.934862 | 3.084422 |     1 |     1 |     1 |     1 | Yes   | Yes     | Yes     | Yes     | Yes     | Yes     | Yes     | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | FSC                 |              50 |      0.59 |                                500 |                     1 | DESKTOP-9MNKV60 | Unknown        | C:/Users/User/Documents/Backup data/2019-03-25\_at\_12-37-08pm.xml | BS4        | good   | No\!     |
|             3 | BS4\_2000  |             6508 |              5000 | Cyanos     |            2000 |              10 | 517.90008 |     12.56613 |             21.29853 | 25-MAR-2019         | 12:45:19            | 10.37472 |   1 | 9.934862 | 3.084422 |     1 |     1 |     1 |     1 | Yes   | Yes     | Yes     | Yes     | Yes     | Yes     | Yes     | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | FSC                 |              50 |      0.59 |                                500 |                     2 | DESKTOP-9MNKV60 | Unknown        | C:/Users/User/Documents/Backup data/2019-03-25\_at\_12-37-08pm.xml | BS4        | good   | Retain   |
|             4 | BS5\_20000 |             5976 |              5000 | Cyanos     |           20000 |              10 |  48.31036 |    123.70018 |            209.66132 | 25-MAR-2019         | 12:49:22            | 10.37472 |   1 | 9.934862 | 3.084422 |     1 |     1 |     1 |     1 | Yes   | Yes     | Yes     | Yes     | Yes     | Yes     | Yes     | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | FSC                 |              50 |      0.59 |                                500 |                     2 | DESKTOP-9MNKV60 | Unknown        | C:/Users/User/Documents/Backup data/2019-03-25\_at\_12-37-08pm.xml | BS5        | bad    | No\!     |
|             5 | BS5\_10000 |             5844 |              5000 | Cyanos     |           10000 |              10 |  90.51666 |     64.56270 |            109.42831 | 25-MAR-2019         | 12:51:49            | 10.37472 |   1 | 9.934862 | 3.084422 |     1 |     1 |     1 |     1 | Yes   | Yes     | Yes     | Yes     | Yes     | Yes     | Yes     | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | FSC                 |              50 |      0.59 |                                500 |                     1 | DESKTOP-9MNKV60 | Unknown        | C:/Users/User/Documents/Backup data/2019-03-25\_at\_12-37-08pm.xml | BS5        | good   | No\!     |
|             6 | BS5\_2000  |             5829 |              5000 | Cyanos     |            2000 |              10 | 400.72498 |     14.54614 |             24.65447 | 25-MAR-2019         | 12:52:47            | 10.37472 |   1 | 9.934862 | 3.084422 |     1 |     1 |     1 |     1 | Yes   | Yes     | Yes     | Yes     | Yes     | Yes     | Yes     | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | —            | FSC                 |              50 |      0.59 |                                500 |                     1 | DESKTOP-9MNKV60 | Unknown        | C:/Users/User/Documents/Backup data/2019-03-25\_at\_12-37-08pm.xml | BS5        | good   | Retain   |

``` r
#some interesting column
knitr::kable(metafile %>% 
              dplyr::select(Sample.ID, Sample.ID2, Number.of.Events, Dilution.Factor, CellspML, Status, Retained) )
```

| Sample.ID  | Sample.ID2 | Number.of.Events | Dilution.Factor |  CellspML | Status | Retained |
| :--------- | :--------- | ---------------: | --------------: | --------: | :----- | :------- |
| BS4\_20000 | BS4        |             6918 |           20000 |  62.02270 | good   | No\!     |
| BS4\_10000 | BS4        |             6591 |           10000 | 116.76311 | good   | No\!     |
| BS4\_2000  | BS4        |             6508 |            2000 | 517.90008 | good   | Retain   |
| BS5\_20000 | BS5        |             5976 |           20000 |  48.31036 | bad    | No\!     |
| BS5\_10000 | BS5        |             5844 |           10000 |  90.51666 | good   | No\!     |
| BS5\_2000  | BS5        |             5829 |            2000 | 400.72498 | good   | Retain   |

Notice that the function suggests you retain only the measurements
associated with dilution *2000* since *make\_decision = “maxi”* and
diltion *2000* happens to have the highest cells/\(\mu\)L measurement
among the dilution levels for both *BS4* and *BS5*. Furthermore, rather
than reading in 6 files, we have narrowed down to reading only the 2
needed files.

## Flow Cytometer File Processing

### Removing NAs in Expression Matrix

`nona()` functions removes `NA` values from the inputed
flowframe.

``` r
flowfile_path <- system.file("extdata", "text.fcs", package = "cyanoFilter",
                              mustWork = TRUE)
flowfile <- flowCore::read.FCS(flowfile_path, alter.names = TRUE,
                               transformation = FALSE, emptyValue = FALSE,
                               dataset = 1) #FCS file contains only one data object
flowfile_nona <- cyanoFilter::nona(x = flowfile)
```

### Removing Negative Values in Expression Matrix

Typically negative values are removed from the expression matrix as they
are deemed measurement error (there are some arguments against this) so
the `noneg()` rids a flowframe of all negative values in its expression
matrix.

``` r
flowfile_noneg <- cyanoFilter::noneg(x = flowfile_nona)
```

### Log-Transforming Expression Matrix

`lnTrans()` transforms all values (except those specified in the
*notToTransform* option) in an expression matrix to the log scale. This
function has a counterpart in the `flowCore()` package but we made
things simpler and also give the opportunity for users to specify
columns in the flowframe that are not to be
transformed.

``` r
flowfile_logtrans <- cyanoFilter::lnTrans(x = flowfile_noneg, notToTransform = c("SSC.W", "TIME"))
```

### Plotting

pair\_plot() gives a scatter plot of all columns in the flowframe,
except that specified in the *notToPlot*
option.

``` r
cyanoFilter::pair_plot(flowfile_noneg, notToPlot = "TIME") ##untransfrmed
```

![](man/figures/README-plotting-1.png)<!-- -->

``` r
cyanoFilter::pair_plot(flowfile_logtrans, notToPlot = "TIME") ##logtransformed
```

![](man/figures/README-plotting-2.png)<!-- -->

## Clustering and Gating

### Margin Events

Margin Events are particles that are too large for teh flow cytometer to
measure. It is desired to eliminate this particles or assign indicators
to them in the flow frame to allow its identification. The
`cellmargin()` function achieves this. It returns a list containing;

  - *fullflowframe* with indicator for margin and non-margin events in
    th eexpression matrix,
  - *reducedflowframe* containing only non-margin events
  - *N\_margin* number of margin events contained in the input flowframe
  - *N\_nonmargin* number of non-margin events
  - *N\_particle* number of particles in the input
flowframe

<!-- end list -->

``` r
flowfile_marginout <- cyanoFilter::cellmargin(flow.frame = flowfile_logtrans, Channel = 'SSC.W',
                                              type = 'estimate', y_toplot = "FSC.HLin")
```

![](man/figures/README-marginEvents-1.png)<!-- -->

### cyanobacteria Population Identification (Kernel Density Approach using flowDensity)

``` r

cyanoFilter::celldebris_nc(flowfile_marginout$reducedflowframe, channel1 = "RED.B.HLin",
                    channel2 = "YEL.B.HLin", interest = "BS4", to_retain = "refined" )
```

![](man/figures/README-unnamed-chunk-1-1.png)<!-- -->

    > $fullframe
    > flowFrame object ' B4_18_1'
    > with 3831 cells and 13 observables:
    >                  name                                desc range
    > $P1          FSC.HLin          Forward Scatter (FSC-HLin) 1e+05
    > $P2          SSC.HLin             Side Scatter (SSC-HLin) 1e+05
    > $P3        GRN.B.HLin   Green-B Fluorescence (GRN-B-HLin) 1e+05
    > $P4        YEL.B.HLin  Yellow-B Fluorescence (YEL-B-HLin) 1e+05
    > $P5        RED.B.HLin     Red-B Fluorescence (RED-B-HLin) 1e+05
    > $P6        NIR.B.HLin Near IR-B Fluorescence (NIR-B-HLin) 1e+05
    > $P7        RED.R.HLin     Red-R Fluorescence (RED-R-HLin) 1e+05
    > $P8        NIR.R.HLin Near IR-R Fluorescence (NIR-R-HLin) 1e+05
    > $P9          SSC.ALin        Side Scatter Area (SSC-ALin) 1e+05
    > $P10            SSC.W          Side Scatter Width (SSC-W) 1e+05
    > $P11             TIME                                Time 1e+05
    > 12   Margin.Indicator                    Margin Indicator     1
    > $P13 BS4BS5.Indicator                    BS4BS5.Indicator     1
    >               minRange maxRange
    > $P1                  0    99999
    > $P2  -34.4792823791504    99999
    > $P3  -21.1945362091064    99999
    > $P4  -10.3274412155151    99999
    > $P5  -5.34720277786255    99999
    > $P6  -4.30798292160034    99999
    > $P7  -25.4901847839355    99999
    > $P8  -16.0200233459473    99999
    > $P9                  0    99999
    > $P10              -111    99999
    > $P11                 0    99999
    > 12                   0        1
    > $P13                 0        1
    > 368 keywords are stored in the 'description' slot
    > 
    > $reducedframe
    > flowFrame object ' B4_18_1'
    > with 3066 cells and 13 observables:
    >                  name                                desc range
    > $P1          FSC.HLin          Forward Scatter (FSC-HLin) 1e+05
    > $P2          SSC.HLin             Side Scatter (SSC-HLin) 1e+05
    > $P3        GRN.B.HLin   Green-B Fluorescence (GRN-B-HLin) 1e+05
    > $P4        YEL.B.HLin  Yellow-B Fluorescence (YEL-B-HLin) 1e+05
    > $P5        RED.B.HLin     Red-B Fluorescence (RED-B-HLin) 1e+05
    > $P6        NIR.B.HLin Near IR-B Fluorescence (NIR-B-HLin) 1e+05
    > $P7        RED.R.HLin     Red-R Fluorescence (RED-R-HLin) 1e+05
    > $P8        NIR.R.HLin Near IR-R Fluorescence (NIR-R-HLin) 1e+05
    > $P9          SSC.ALin        Side Scatter Area (SSC-ALin) 1e+05
    > $P10            SSC.W          Side Scatter Width (SSC-W) 1e+05
    > $P11             TIME                                Time 1e+05
    > 12   Margin.Indicator                    Margin Indicator     1
    > $P13 BS4BS5.Indicator                    BS4BS5.Indicator     1
    >               minRange maxRange
    > $P1                  0    99999
    > $P2  -34.4792823791504    99999
    > $P3  -21.1945362091064    99999
    > $P4  -10.3274412155151    99999
    > $P5  -5.34720277786255    99999
    > $P6  -4.30798292160034    99999
    > $P7  -25.4901847839355    99999
    > $P8  -16.0200233459473    99999
    > $P9                  0    99999
    > $P10              -111    99999
    > $P11                 0    99999
    > 12                   0        1
    > $P13                 0        1
    > 368 keywords are stored in the 'description' slot
    > 
    > $Cell_count
    > [1] 3066
    > 
    > $Debris_Count
    > [1] 348

### cyanobacteria Population Identification (EM Approach)

``` r

cyanoFilter::celldebris_emclustering(flowfile_marginout$reducedflowframe, channels =  c("RED.B.HLin",
                    "YEL.B.HLin", "FSC.HLin", "RED.R.HLin") )
```

![](man/figures/README-emapproach-1.png)<!-- -->

    > $percentages
    > [1] 0.16141975 0.05674587 0.03567496 0.10603963 0.64011979
    > 
    > $mus
    >                 [,1]      [,2]      [,3]     [,4]     [,5]
    > RED.B.HLin 3.7667113 0.6918045 0.8377841 4.672208 3.940687
    > YEL.B.HLin 0.8290632 2.0195027 1.5998735 2.189597 1.847443
    > FSC.HLin   4.7099883 5.5553504 5.1949117 6.799144 4.686122
    > RED.R.HLin 7.5222036 2.1122156 3.0892964 8.460369 7.739509
    > 
    > $sigmas
    > $sigmas[[1]]
    >             RED.B.HLin  YEL.B.HLin   FSC.HLin  RED.R.HLin
    > RED.B.HLin  0.08650772 -0.04124413 0.07685445  0.08766450
    > YEL.B.HLin -0.04124413  1.49207964 0.16748962 -0.06058480
    > FSC.HLin    0.07685445  0.16748962 0.20253901  0.08075963
    > RED.R.HLin  0.08766450 -0.06058480 0.08075963  0.12208235
    > 
    > $sigmas[[2]]
    >            RED.B.HLin   YEL.B.HLin   FSC.HLin   RED.R.HLin
    > RED.B.HLin 0.24587507  0.064379340 0.15697100  0.039666562
    > YEL.B.HLin 0.06437934  0.158786398 0.13088514 -0.006178874
    > FSC.HLin   0.15697100  0.130885137 1.07100822  0.020786312
    > RED.R.HLin 0.03966656 -0.006178874 0.02078631  0.362570873
    > 
    > $sigmas[[3]]
    >            RED.B.HLin YEL.B.HLin  FSC.HLin RED.R.HLin
    > RED.B.HLin  1.3504202  0.1361328 0.1780014  1.0482982
    > YEL.B.HLin  0.1361328  0.8100661 0.1622922  0.4575789
    > FSC.HLin    0.1780014  0.1622922 0.9449451  0.4280339
    > RED.R.HLin  1.0482982  0.4575789 0.4280339  5.3530708
    > 
    > $sigmas[[4]]
    >            RED.B.HLin YEL.B.HLin  FSC.HLin RED.R.HLin
    > RED.B.HLin  0.6641398  0.2581150 0.6986730  0.6627531
    > YEL.B.HLin  0.2581150  0.3525741 0.4785388  0.2572977
    > FSC.HLin    0.6986730  0.4785388 1.5192332  0.6987817
    > RED.R.HLin  0.6627531  0.2572977 0.6987817  0.6758216
    > 
    > $sigmas[[5]]
    >            RED.B.HLin  YEL.B.HLin   FSC.HLin  RED.R.HLin
    > RED.B.HLin 0.08377324 0.015955755 0.09307965 0.075609542
    > YEL.B.HLin 0.01595576 0.245660497 0.02724529 0.008981062
    > FSC.HLin   0.09307965 0.027245289 0.15232267 0.094490988
    > RED.R.HLin 0.07560954 0.008981062 0.09449099 0.091592072
    > 
    > 
    > $result
    > flowFrame object ' B4_18_1'
    > with 3831 cells and 17 observables:
    >                  name                                desc range
    > $P1          FSC.HLin          Forward Scatter (FSC-HLin) 1e+05
    > $P2          SSC.HLin             Side Scatter (SSC-HLin) 1e+05
    > $P3        GRN.B.HLin   Green-B Fluorescence (GRN-B-HLin) 1e+05
    > $P4        YEL.B.HLin  Yellow-B Fluorescence (YEL-B-HLin) 1e+05
    > $P5        RED.B.HLin     Red-B Fluorescence (RED-B-HLin) 1e+05
    > $P6        NIR.B.HLin Near IR-B Fluorescence (NIR-B-HLin) 1e+05
    > $P7        RED.R.HLin     Red-R Fluorescence (RED-R-HLin) 1e+05
    > $P8        NIR.R.HLin Near IR-R Fluorescence (NIR-R-HLin) 1e+05
    > $P9          SSC.ALin        Side Scatter Area (SSC-ALin) 1e+05
    > $P10            SSC.W          Side Scatter Width (SSC-W) 1e+05
    > $P11             TIME                                Time 1e+05
    > 12   Margin.Indicator                    Margin Indicator     1
    > 1      Cluster_Prob_1                      Cluster_Prob_1     1
    > 2      Cluster_Prob_2                      Cluster_Prob_2     1
    > 3      Cluster_Prob_3                      Cluster_Prob_3     1
    > 4      Cluster_Prob_4                      Cluster_Prob_4     1
    > 5      Cluster_Prob_5                      Cluster_Prob_5     1
    >                   minRange          maxRange
    > $P1                      0             99999
    > $P2      -34.4792823791504             99999
    > $P3      -21.1945362091064             99999
    > $P4      -10.3274412155151             99999
    > $P5      -5.34720277786255             99999
    > $P6      -4.30798292160034             99999
    > $P7      -25.4901847839355             99999
    > $P8      -16.0200233459473             99999
    > $P9                      0             99999
    > $P10                  -111             99999
    > $P11                     0             99999
    > 12                       0                 1
    > 1                        0 0.999999999999675
    > 2    2.00292214577102e-109 0.996307497800057
    > 3     3.24771377094787e-13                 1
    > 4                        0 0.999999999769265
    > 5                        0 0.994124337680755
    > 368 keywords are stored in the 'description' slot

# License

This is a free to use package for anyone who has the need.
