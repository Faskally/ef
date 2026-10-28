#' Estimating overdispersion
#'
#' Complete function for returning overdispersion estimates
#'
#'
#' @param data dataframe containing EF data
#' @param visitID a number identifying each unique site visit
#' @param count the number of fish caught for a particular combination of site visit, 
#'              species, lifestage and pass (defaults to "count")
#' @param pass EF pass number e.g. 1,2,3,4 (defaults to "pass")
#' @param sampleID sample ID i.e. unique combinations of site visit, species & lifestage
#' @param largemodel a large model that captures most of the systematic
#'                variation in the data - this is specified before
#'               running the overdispersion function
#' @param control passes control information to optimiser              
#' @return a data.frame summarising overdispersion
#'
#' @note ensure column names in function call are in inverted commas
#'
#' @importFrom dplyr left_join
#' @importFrom dplyr n_distinct
#' @importFrom stats aggregate pchisq qlogis fitted
#'
#' @export
overdispersion <- function(data, visitID = "visitID", count = "count",
                           pass = "pass", sampleID = "sampleID", largemodel, control = "control") {

  # organise data set with standard names
  # varID are the supplied names
  # varID2 are the standard names we will use in the code 
  
  varID <- c(visitID, count, pass, sampleID)
  varID2 <- c("visitID", "count", "pass", "sampleID")
  
  # check mandatory names are present
  if (!all(varID %in% names(data)))
    stop("some mandatory variables are missing")
  
  # subset data
  data <- data[varID]
  
  # rename columns
  names(data) <- varID2
  
  # set visitID to a factor
  data$visitID <- as.factor(data$visitID)
  data$sampleID <- as.factor(data$sampleID)
  
  # return the maximum number of EF passes for each visit
  # checks for missing passes
  # checks for visits with only one pass
  
  # first check that pass is numeric
  if (!is.numeric(data$pass))
    stop("pass should be numeric")
  
  aggregatepass <- aggregate(pass ~ visitID, data = data, FUN = function(x) {
    if (max(x) != length(unique(x))) 
      stop("some visits have missing pass information")
    if (max(x) == 1)
      stop("some visits have only one pass")
    max(x)
  })
  
  # rename the columns
  names(aggregatepass) <- c("visitID", "maxpass")
  
  # set the visitID to be a factor
  aggregatepass$visitID <- as.factor(aggregatepass$visitID)
  
  # join the data and the maximum number of passes together. Use
  # left_join to maintain the dataframe order
  data <- dplyr::left_join(data, aggregatepass, by = "visitID")
  
  ################
  #
  # SATURATED MODEL
  #
  ################
  
  # Produce a saturated model for each visit. 
  # '-1' prevents the model estimate of intercept, and avoids problems if the 
  # intercept for the first sample is poorly estimated. This results in a 
  # separate estimate of intercept for every sample. 
  # Hessian = FALSE prevents the return of the hessian matrix (used to calculate 
  # standard errors) and speeds up model fitting.
  
  
  # Record and print the system time
  t1 <- Sys.time()
  message("Saturated model start time ", t1)
  
  saturated <- by(data, data$visitID, function(x) {
    
    # drop unused factor levels in data, in particular sampleID which is used
    # in both the model formula and directly by efp
    x <- droplevels(x)
    
    # check there is exactly one observation for each sampleID and pass
    if (!all(table(x$sampleID, x$pass) == 1)) 
      stop(
        "unbalanced data: number of passes differs betwen samples within visits"
      )	
    
    # drop samples with no fish; these have no information but their inclusion
    # would bias the dispersion estimate downwards because the effective number 
    # of parameters is inflated
    
    n_fish <- tapply(x$count, x$sampleID, sum)
    
    if (any(n_fish == 0)) {
      keep_id <- names(n_fish)[n_fish > 0]
      x <- x[x$sampleID %in% keep_id, ]
      x <- droplevels(x)
    }
    
    # create pass_reduce, a variable which combines the last two 
    # passes into a single category: 
    # e.g. if there are three passes, then pass_reduce is the same as pass12
    x$pass_reduce <- pmin(x$pass, x$maxpass - 1)
    
    # combine this with sampleID to get a saturated model formula
    x$groupID <- paste(x$sampleID, x$pass_reduce)
    
    
    # fit saturated model
    
    if(dplyr::n_distinct(x$groupID)>1) {
      
      formula <- count ~ -1 + groupID 
      
    } else {
      
      formula <- count ~ 1
      
    }
    
    m <- efp(
      formula, 
      pass = pass, id = sampleID, data = x, 
      control = control,
      hessian = FALSE, verbose = FALSE)
    
    data.frame(llik = m$llik, nparam = length(m$coefficients))
  })      
  
  # combine into a single data frame
  saturated <- do.call("rbind", saturated)
  
  # sum over all visits
  saturated <- data.frame(
    llik = sum(saturated$llik), 
    nparam = sum(saturated$nparam)
  )
  
  
  # store and print the time that it took for the model to run
  t2 <- Sys.time()
  etime <- difftime(t2, t1, units = "secs")
  message("Saturated model duration = ", round(etime, 0), " secs")
  

  ################
  #
  # SITE-VISIT MODEL
  #
  ################


  # The site-visit model has an extra parameter for each site-visit (over the large model). 
  # Due to the number of  parameters, this model is fitted by conditioning on the fitted
  # values from the large model and estimating the site-visit effect for each
  # site-visit in turn.

  # The site-visit model includes estimates of all of the parameters from the full model
  # as the offset, thus the additional variability due to the inclusion of the site-visits
  # gives you an estimate of the between visit overdispersion.


  # Record and print system time
  t1 <- Sys.time()
  message("Site visit model start time ", t1)
  
  # create a column which has the transformed p values from the large model
  data$p <- qlogis(fitted(largemodel))
  
  # estimate an intercept for each site visit
  sitevisit <- by(data, data$visitID, function(x) {
    x <- droplevels(x)
    m <-
      efp(
        count ~ 1,
        pass = pass, id = sampleID, data = x,
        control = control,
        offset = x$p, hessian = FALSE, verbose = FALSE
      )
    
    data.frame(llik = m$llik, nparam = length(m$coefficients))
  })
  
  # combine into a single data frame 
  sitevisit <- do.call("rbind", sitevisit)
  
  # sum over all visits and add on number of parameters from large model
  sitevisit <- data.frame(
    llik = sum(sitevisit$llik), 
    nparam = sum(sitevisit$nparam) + length(largemodel$coefficients)
  )
  
  # store and print the system time and calculate the running time for the site-visit model
  t2 <- Sys.time()
  etime <- difftime(t2, t1, units = "secs")
  message("Site visit model duration = ", round(etime, 0), " secs")
  
  ################
  #
  # ESTIMATING OVERDISPERSION
  #
  ################

  # extract the log-likelihood and number of parameters from the large model
  largeout <-
    data.frame(
      llik = largemodel$llik,
      nparam = length(largemodel$coefficients)
    )

  # bind together the log-likelihoods and number of parameters from each
  # of the three modelling processes
  anova <- rbind(saturated, sitevisit, largeout)
  row.names(anova) <- c("saturated", "sitevisit", "large")
  
  # calculate the difference in the deviance and degrees of freedom between successive
  # model outputs and estimate the within visit (saturated to site-visit model) and between
  # visit (site-visit to large) overdispersion
  anova$deviance <- c(NA, -2 * diff(anova$llik))
  anova$df <- c(NA, -diff(anova$nparam))
  anova$disp <- anova$deviance / anova$df
  anova$p <- pchisq(anova$deviance, anova$df, lower.tail = FALSE)
  
  # return the overdispersion estimates and p values
  # the values for "disp" are then included in the adjusted BIC function
  return(anova)
}
