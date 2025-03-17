if (!requireNamespace("bexy", quietly = TRUE)) {
    install.packages("bexy")
}
library(bexy)

bex <- bexy("./bexy/")

getPosteriorModeSexKaryotypes(bex)

result_dir <- "bexy_out"

dir.create(result_dir), showWarnings = FALSE)

writePosteriorModeSexKaryotypes(bex, './bexy_output', threshold_certainty = 0.95)

