
# ------------------------------------------------------------------------------
# this file is supposed to test the ability to predict from a
# capture probability model object when new data is available
# ------------------------------------------------------------------------------

# load ef package
#httr::set_config(httr::use_proxy(url="192.168.41.8", port=80))
pkg <- devtools::as.package("C:/work/repos/Faskally/ef")
devtools::document(pkg)
devtools::load_all(pkg)

# load other libs
library(tidyr)
library(dplyr)
library(magrittr)


# ------------------------------------------------------------------------------
# get some data
# ------------------------------------------------------------------------------

# connect to database
db <- "B:/Conservation_Limits/CL_Juvenile_Density/CollateData/db/scottishEF_CLROAME_v01.sqlite3"
con <- DBI::dbConnect(RSQLite::SQLite(), db)
# get data
ef <- DBI::dbReadTable(con, "ef"); DBI::dbDisconnect(con)

# reformat - fobs sites with no site ID are welsh sites etc.  Need to check other sites
ef %<>% #as_data_frame() %>%
  filter(!is.na(Site_OBJECTID) &
           !is.na(Area)          &
           Runs >= 2             &
           Stocked == "No")

# restructure ef
ef <- do.call(rbind,
              lapply(c("S0", "SP", "T0", "TP"),
                     function(x) {
                       out <- ef[c("Site_OBJECTID", "Site.Name", "Dataset", "Width", "BedWidth", "Date", "Runs", "Area", "Trust", paste0(x, "_R", 1:6))]
                       names(out) <- gsub(x, "n", names(out))
                       out $ Species <- if (substring(x, 1, 1) == "S") "Salmon" else "Trout"
                       out $ LifeStage <- if (substring(x, 2, 2) == "0") "Fry" else "Parr"
                       out
                     }
              ))

# drop rows with all NAs
# keep only salmon fry:
ef <- filter(ef, Species == "Salmon" & LifeStage == "Fry")

# replace all NA where there should be an observation with a zero observation
for (i in 1:nrow(ef)) {
  fill <- 1:6 %in% 1:ef $ Runs[i]
  vals <- unlist(ef[i, paste0("n_R", 1:6)])
  vals[!fill] <- NA
  vals[fill] <- ifelse(is.na(vals[fill]), 0, vals[fill])
  ef[i, paste0("n_R", 1:6)] <- vals
}
rm(vals)

# get site info
con <- DBI::dbConnect(RSQLite::SQLite(), db)
gis <- DBI::dbReadTable(con, "gis"); DBI::dbDisconnect(con)
# define points as spatial in BNG
gis <- sp::SpatialPointsDataFrame(
  coords = cbind(gis $ NEAR_X, gis $ NEAR_Y),
  data = gis,
  coords.nrs = c(1,2),
  proj4string = sp::CRS("+init=epsg:27700"))

# tag on HMA data to sites
data(hma, package = "CLdata")
raster::crs(hma) <- raster::crs(gis)
hma <- hma[!(hma $ HAName %in% c("Shetlands", "Orkneys")),]
hma $ hmidx <- 1:nrow(hma)
gis @ data <- cbind(gis @ data, sp::over(gis, hma))

# add in catchment info to sites
data(redctm, package = "CLdata")
raster::crs(redctm) <- raster::crs(gis)
gis @ data <- cbind(gis @ data, sp::over(gis, redctm))

# tag on site info data
ef <- cbind(ef, gis[ef $ Site_OBJECTID,])
ef $ pDate <- as.POSIXlt(ef $ Date, tz = "GMT", format = "%d/%m/%Y")

# add in barrier info.... on the fly by taking from CLdata gis data
warning("Barrer info is being taken from CLdata package.  Any new data is assumed to be below barrier.")
gis2 <- CLdata::gis @ data
row.names(gis2) <- paste(gis2 $ NEAR_X, gis2 $ NEAR_Y)
gis $ barrier <- gis2[paste(gis $ NEAR_X, gis $ NEAR_Y),] $ barrier
# assume new data (not in CLdata) is below barriers
gis $ barrier[gis $ CATCH_ == -1] <- FALSE


# tidy covariates
ef $ year <- lubridate::year(ef $ pDate)
ef $ doy <- lubridate::yday(ef $ pDate)
ef $ CATCH_ID <- factor(ef $ CATCH_ID)
ef $ fyear <- factor(ef $ year)
ef $ sinSlope <- sin(ef $ Slope_deg/180*pi)
landuse <- c("CTrees", "Urban", "NCTrees", "Mixed", "Marsh", "Other")
ef $ totlanduse <- rowSums(ef[landuse])
ef[landuse] <- ef[landuse] / ef $ totlanduse
ef $ Distance_s <- ef $ Distance_s / 1000

## remove covariate outliers
ef $ keep <- with(ef, doy > 150 & doy < 325 &
                    year >= 1997 & year <= 2014 &
                    Area < 5000 &
                    Species == "Salmon" & LifeStage == "Fry" &
                    Water_W < 35 &
                    Elevation_ < 500 &
                    sinSlope < 0.4 &
                    Upcatch_km < 600)
ef <- subset(ef, keep)
ef <- ef[names(ef) != "pDate"]

# organise by pass
ef_byp <- tidyr::gather(ef, pass, n, n_R1:n_R6)
ef_byp$pass <- as.numeric(gsub("n_R", "", as.character(ef_byp$pass)))
ef_byp <- subset(ef_byp, pass < 3)

rm(gis, gis2, fill, hma, landuse, redctm, con, db)

# ------------------------------------------------------------------------------
# Fit a model
# ------------------------------------------------------------------------------
ef_byp $ id <- with(ef_byp, factor(paste(Site_OBJECTID, Date, Species, LifeStage)))

mod1 <- efp(n ~ 1, pass = pass, id = id, data = ef_byp)



