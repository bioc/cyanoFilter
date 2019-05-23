
# cyanoFilter

cyaoFilter is a package designed to identify, assign indicators and/or
filter out synechoccus type cyanobacteria from a water sample examined
via flowcytometry.

# Installation and Dependencies

To use the package, please intall the following needed **bioconductor**
packages.

``` r
install.packages("BiocManager")
library(BiocManager)
install(c("Biobase", "flowCore", "flowDensity"))
```

All other needed dependencies are avialble on **CRAN** and will be
installed when you install the package itself.

``` r
install.packages("cyanoFilter")
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

# License

This is a free to use package for anyone who has the need.
