## ---- sim-design
set.seed(2024)

nsite <- 30
years <- 2006:2020

# one visit to each site in each year
visits <- expand.grid(siteID = paste0("Site", 1:nsite), year = years)
visits$siteID <- factor(visits$siteID, levels = paste0("Site", 1:nsite))
nvisit <- nrow(visits)
visits$visitID <- factor(paste(visits$siteID, visits$year))

# sites are fished between June and October
visits$doy <- sample(152:304, nvisit, replace = TRUE)

# the area fished (m^2) is a property of the site
visits$area <- round(runif(nsite, 80, 250))[visits$siteID]

## ---- sim-effects
# smooth effects on logit capture probability
f_year <- function(year) 0.5 * sin(2 * pi * (year - 2006) / 12)
f_doy <- function(doy) 0.9 * exp(-((doy - 200) / 35)^2) - 0.45

## ---- sim-samples
# two samples per visit: fry and parr
samples <- merge(visits, data.frame(lifestage = c("Fry", "Parr")))
samples <- samples[order(samples$visitID, samples$lifestage), ]
samples$lifestage <- factor(samples$lifestage)
samples$sampleID <- factor(paste(samples$visitID, samples$lifestage))
nsample <- nrow(samples)

# random effects on logit capture probability
sigma_visit <- 0.4
sigma_sample <- 0.3
visit_effect <- rnorm(nvisit, 0, sigma_visit)
names(visit_effect) <- levels(visits$visitID)
samples$visit_effect <- visit_effect[as.character(samples$visitID)]
samples$sample_effect <- rnorm(nsample, 0, sigma_sample)

# the number of fish in the area fished
density <- c(Fry = 0.8, Parr = 0.4) # fish per m^2
site_effect <- rnorm(nsite, 0, 0.5)
samples$N <-
  rpois(
    nsample,
    density[as.character(samples$lifestage)] * samples$area *
      exp(site_effect[samples$siteID])
  )

## ---- sim-passes
# three passes per sample
ef_sim <- merge(samples, data.frame(pass = 1:3))
ef_sim <- ef_sim[order(ef_sim$sampleID, ef_sim$pass), ]
rownames(ef_sim) <- NULL

# capture probability is lower after the first pass, and higher for parr
ef_sim$pass12 <- factor(pmin(ef_sim$pass, 2))
ef_sim$eta <-
  0.3 +
  0.5 * (ef_sim$lifestage == "Parr") -
  0.4 * (ef_sim$pass12 == "2") +
  f_year(ef_sim$year) +
  f_doy(ef_sim$doy) +
  ef_sim$visit_effect +
  ef_sim$sample_effect
ef_sim$p <- 1 / (1 + exp(-ef_sim$eta))

## ---- sim-removal
# remove fish pass by pass
ef_sim$count <- NA
remaining <- samples$N
names(remaining) <- samples$sampleID
for (i in 1:3) {
  ipass <- which(ef_sim$pass == i)
  id <- as.character(ef_sim$sampleID[ipass])
  ef_sim$count[ipass] <- rbinom(length(ipass), remaining[id], ef_sim$p[ipass])
  remaining[id] <- remaining[id] - ef_sim$count[ipass]
}

## ---- sim-factors
# factors for random walk smooths: every step of the walk is a level
ef_sim$fyear <- factor(ef_sim$year, levels = years)
days <- min(ef_sim$doy):max(ef_sim$doy)
ef_sim$fdoy <- factor(ef_sim$doy, levels = days)
