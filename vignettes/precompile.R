# Vignettes that are computationally intensive have been precompiled:
# Must manually move image files from BRCore/figure to BRCore/vignettes/ after knit

library(knitr)
library(here)
library(devtools)

knit("vignettes/BRCore-vignette.Rmd.orig", "vignettes/BRCore-vignette.Rmd")

# Moving files
figure_files <- list.files(here::here("images"), full.names = TRUE)
figure_dir <- here::here("images")

# Move files to vignettes and doc directories
file.copy(figure_files, here::here("vignettes/images"), recursive = TRUE)
file.copy(figure_files, here::here("doc/images"), recursive = TRUE)

# Cleanup
file.remove(figure_files, figure_dir)

# Build static vignette
build_vignettes()
