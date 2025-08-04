#' Reporting graphical comparison procedures in GSD studies
#'
#' appendMCP provides routines for describing graphical multiple comparison
#' procedures (MCP) in group sequential design (GSD) studies. The objective is
#' to generate an R Markdown file that defines the MTP and interim analyses,
#' which can be appended to a statistical analysis plan document.
#'
#' The package provides predefined templates that can be modified by a user,
#' ideally by changing only high-level input parameters.
#'
#' To load a template in Rstudio, use "File" -> "New File" -> "R Markdown ..."
#' Select "From Template" (left panel), then choose "MCP in GSD"
#'
#' @examples
#' # Load template from command line:
#' # rmarkdown::draft("my_report", template = "mcpgen", package = "appendMCP", edit = FALSE)
#' # To open created file in Rstudio editor use
#' # rstudioapi::navigateToFile("my_report/my_report.Rmd")
#'
#' @name appendMCP
#' @docType package
"_PACKAGE"

