## Plot line plots between 2 more time points using ggplot2

# Libraries
requiredPackages = c("tidyverse",)
if (length(setdiff(requiredPackages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(requiredPackages, rownames(installed.packages())))  
} else { 
  lapply( requiredPackages , require, character.only = TRUE) 
}
