required_packages <- c(
  "devtools", "ggplot2"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

invisible(lapply(required_packages, install_if_missing))
