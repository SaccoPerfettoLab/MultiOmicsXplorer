# To run MultiOmicsXplorer is necessary to install all the following packages.
# Use this code to have all the packages installed! When you'll run the application the libraries 
# will be uploaded automatically

# Write to eleonorameo.hp@gmail.com for any problem

required_packages <- c(
  "shiny",
  "shinydashboard",
  "shinyBS",
  "dplyr",
  "tidyverse",
  "ggplot2",
  "ggpubr",
  "fs",
  "ggsignif",
  "shinycssloaders",
  "png",
  "plotly",
  "gridExtra",
  "shinyjs",
  "DT",
  "openxlsx",
  "plumber"
)

# install missing CRAN packages
installed <- installed.packages()[, "Package"]
for (pkg in required_packages) {
  if (!pkg %in% installed) {
    message("Installing ", pkg, "...")
    install.packages(pkg, dependencies = TRUE)
  }
}

# install SignalingProfiler from GitHub if not already installed
if (!"SignalingProfiler" %in% installed) {
  if (!"devtools" %in% installed) install.packages("devtools")
  devtools::install_github('https://github.com/SaccoPerfettoLab/SignalingProfiler/')
  }
