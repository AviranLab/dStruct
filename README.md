# dStruct

Method for identifying differential reactive regions from RNA structurome profiling data.

dStruct is a statistical method for identifying regions that display altered reactivity patterns between two groups of samples. It can perform *de novo* discovery or identify regions from a list provided by the user. The latter case is called *guided discovery*.

## Getting started

dStruct is written in R version 3.4.1. To start with, download and install the latest versions of [R](https://cran.r-project.org/) and [RStudio](https://www.rstudio.com/products/rstudio/).

dStruct can be installed directly from source. Download the distribution tarball. Next, in RStudio type

`install.packages(path_to_tarball, repos = NULL, type="source")`.

This should install the package in R. Check by executing the following command.

`library(dStruct)`

In the future, we plan to make dStruct available through Bioconductor.

## Usage

dStruct takes reactivities of groups of samples for all transcripts under consideration.

## Citation

Choudhary K., Lai Y.H., Tran E., Aviran S. (2019) "dStruct: identifying differentially reactive regions from RNA structurome profiling data."

## Contributors

* [Krishna Choudhary](https://www.linkedin.com/in/datamaster/) - _Initial implementation and developer_
* [Sharon Aviran](https://aviranlab.bme.ucdavis.edu/) - _Supervisor_


