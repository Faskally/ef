


pkg <- devtools::as.package("~/work/faskally/ef")
pkg <- devtools::as.package("c:/work/repos/faskally/ef")


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

### first test

ef_data <- data.frame(n     = c(100, 53, 24, 50, 26, 12),
                      pass  = c(  1,  2,  3,  1,  2,  3),
                      stage = c(  1,  1,  1,  2,  2,  2))
ef_data

# Fit a simple model
data <- ef_data
formula <- n ~ 1 + factor(stage)
hessian = TRUE 
verbose = TRUE
init = "0"

m1 <- efp(n ~ 1 + factor(stage), data = ef_data)
m1


### second test

ef_data <- data.frame(n     = c(100, 53, 24, 50, 26, 12),
                      pass  = c(  1,  2,  4,  1,  2,  4),
                      stage = c(  1,  1,  1,  2,  2,  2))
ef_data

# Fit a simple model
devtools::load_all(pkg)

m1 <- efp(n ~ 1, data = ef_data, pass = pass)
m1



### test 3

ef_data <- data.frame(n     = c(100, 53, 24, 50, 26, 12),
                      pass  = c(  1,  2,  3,  1,  2,  3),
                      stage = c(  1,  1,  1,  2,  2,  2))
ef_data

# Fit a simple model
devtools::load_all(pkg)

m1 <- efp(n ~ 1, data = ef_data, pass = pass, groups = stage)
m1

### test 4
N <- 200
ef_data <- data.frame(n     = round(N * c(0.5, 0.5^2, 0.5^3, 0.7, 0.3 * 0.7, 0.3^2 * 0.7)),
                      pass  = c(  1,  2,  3,  1,  2,  3),
                      stage = c(  1,  1,  1,  2,  2,  2))
ef_data

# Fit a simple model
devtools::load_all(pkg)

efp(n ~ 1, data = ef_data, pass = pass, groups = stage)
efp(n ~ factor(stage), data = ef_data, pass = pass, groups = stage)
# linear pass model
efp(n ~ pass, data = ef_data, pass = pass, groups = stage)
# pass 2 == pass 3
ef_data $ pass23 <- with(ef_data, replace(pass, pass>2, 2)) 
efp(n ~ factor(pass23), data = ef_data, pass = pass, groups = stage)


### test 5

# get data from CLdata
library(CLdata)
library(tidyr)
library(dplyr)
library(magrittr)
data(ef)

ef %<>% as_data_frame() %>%
        filter(!is.na(Site_OBJECTID) & 
               !is.na(Area)          &
               Runs >= 2             &
               Stocked == "No")

# restructure ef
ef <- do.call(rbind, 
        lapply(c("S0", "SP", "T0", "TP"), 
            function(x) {
              out <- ef[c("Site_OBJECTID", "Site.Name", "Dataset", "Width", "Date", "Runs", "Area", "Trust", paste0(x, "_R", 1:6))]
              names(out) <- gsub(x, "n", names(out))
              out $ Species <- if (substring(x, 1, 1) == "S") "Salmon" else "Trout"
              out $ LifeStage <- if (substring(x, 2, 2) == "0") "Fry" else "Parr"
              out
            }
        ))

# drop rows with all NAs
# keep only salmon:
ef <- filter(ef, Species == "Salmon")

# replace all NA where there should be an observation with a zero observation
for (i in 1:nrow(ef)) {
  fill <- 1:6 %in% 1:ef $ Runs[i]
  vals <- unlist(ef[i, paste0("n_R", 1:6)])
  vals[!fill] <- NA
  vals[fill] <- ifelse(is.na(vals[fill]), 0, vals[fill])
  ef[i, paste0("n_R", 1:6)] <- vals
}

# keep only 3 pass:
ef <- filter(ef, Runs > 2)
ef <- ef[!names(ef) %in% c("n_R4", "n_R5", "n_R6")]

# tag on HMA data
data(hma)
hma <- hma[!(hma $ HAName %in% c("Shetlands", "Orkneys")),]
hma $ hmidx <- 1:nrow(hma)
gis <- CLdata::gis
gis @ data <- cbind(gis @ data, sp::over(gis, hma))

# tag on gis data
ef <- cbind(ef, gis[ef $ Site_OBJECTID,])
# fix missing value
ef $ n_R3[is.na(ef $ n_R3) & ef $ Runs == 3] <- 0

# gather
ef <- gather(ef, pass, n, n_R1:n_R3)

# add in some dates
ef $ pDate <- as.POSIXlt(ef $ Date, tz = "GMT", format = "%d/%m/%Y")
# try using lubridate
ef $ year <- ef $ pDate $ year + 1900
ef $ doy <- ef $ pDate $ yday

ef $ pass <- as.numeric(gsub("n_R", "", ef $ pass))


# fit a model

efp(n ~ 1, data = ef, pass = pass, groups = LifeStage)

efp(n ~ factor(LifeStage), data = ef, pass = pass, groups = LifeStage)
# linear pass model
efp(n ~ pass, data = ef, pass = pass, groups = LifeStage)
# pass 2 == pass 3
ef $ pass23 <- with(ef, replace(pass, pass>2, 2)) 
efp(n ~ factor(pass23), data = ef, pass = pass, groups = LifeStage)



BIC(efp(n ~ LifeStage + factor(pass23), data = ef, pass = pass, groups = LifeStage, verbose = FALSE))
BIC(efp(n ~ LifeStage * factor(pass23), data = ef, pass = pass, groups = LifeStage, verbose = FALSE))



