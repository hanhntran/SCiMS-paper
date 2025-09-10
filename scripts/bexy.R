if (!requireNamespace("bexy", quietly = TRUE)) {
    install.packages("bexy", repos = "http://cran.us.r-project.org")
}
library(bexy)
args <- commandArgs(TRUE) 
files_dir <- args[[1]]

bex <- bexy(files_dir)
getPosteriorModeSexKaryotypes(bex)

writePosteriorModeSexKaryotypes(bex, file.path(files_dir, 'bexy_output_0.95.txt'), threshold_certainty = 0.95)
