
<!-- README.md is generated from README.Rmd. Please edit that file -->

# appendMCP

<!-- badges: start -->

<!-- badges: end -->

The goal of appendMCP is to provide routines for describing graphical
multiple comparison procedures (MCPs) in group sequential design (GSD)
studies. Its primary objective is to generate an R Markdown file that
defines the testing procedure and interim analyses, which can be
appended to a statistical analysis plan (SAP) document. The package
includes predefined templates that users can customize, ideally by
modifying only a small set of high-level input parameters.

## Installation

You can install the development version of appendMCP from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("ytymofyeyev/appendMCP")
```

## Example

To load a template in RStudio, go to **File** $\rightarrow$ **New File**
$\rightarrow$ **R Markdown …** , select **From Template** (on left
panel), then choose “**MCP in GSD**” from the list. You can also load a
template from the command line:

``` r
# install.packages("rmarkdown")
rmarkdown::draft("my_report", template = "mcpgen", package = "appendMCP")

# To open created file in Rstudio editor use
rstudioapi::navigateToFile("my_report/my_report.Rmd")
```

After that, you can knit the R Markdown file to generate a report.

(`library(appendMCP)` is called automatically within the file).

The report describes the trial design and assumptions, interim analyses,
hypothesis testing, and other relevant information.
