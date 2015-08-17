


pkg <- devtools::as.package("~/work/faskally/ef")

## tools:
devtools::document(pkg)
devtools::check(pkg)
devtools::install(pkg)
devtools::load_all(pkg)


# -------------------------------------------
# ef package development

#devtools::use_vignette("Getting-started", pkg)
devtools::install(pkg, build_vignettes = TRUE)

library(ef)
vignette("Getting-started", package = "ef")

## tests

ef_data <- data.frame(n     = c(100, 53, 24, 50, 26, 12),
                   pass  = c(  1,  2,  3,  1,  2,  3),
                   stage = c(  1,  1,  1,  2,  2,  2))
ef_data

# Fit a simple model
m1 <- efp(n ~ 1 + factor(stage), data = ef_data)
