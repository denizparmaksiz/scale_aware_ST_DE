# Install the packages required by tutorial_ALDEx3.Rmd.
#
# Run this script once from the repository root:
# source("documentation_code/aldex_main_crossstudy_repo/tutorial/install_tutorial_dependencies.R")

cran_packages <- c(
  "ggplot2",
  "knitr",
  "lme4",
  "reformulas",
  "remotes",
  "rmarkdown",
  "writexl"
)

missing_packages <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Install the exact standard ALDEx3 commit used by the manuscript. Installing
# by commit avoids silently substituting a later branch or release.
remotes::install_github(
  "jsilve24/ALDEx3",
  ref = "8c05ad40c41279dffa05dc808167ffcd53207740",
  upgrade = "never"
)

message("Tutorial dependencies installed.")
