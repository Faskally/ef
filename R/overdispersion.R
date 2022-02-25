#' Estimating overdispersion
#'
#' Complete function for returning overdispersion estimates
#'
#'
#' @param data dataframe containing EF data
#' @param siteID site name or unique ID
#' @param visitID a number identifying each unique visit
#' @param count the number of fish caught for a particular combination of site visit, 
#'              species, lifestage and pass (defaults to "count")
#' @param pass EF pass number e.g. 1,2,3,4 (defaults to "pass")
#' @param species fish species, e.g. salmon, trout (defaults to "species")
#' @param lifestage lifestage e.g. fry, parr (defaults to "lifestage")
#' @param id sample ID i.e. unique combinations of site visit, species & lifestage
#' @param largemodel a large model that captures most of the systematic
#'                variation in the data - this is specified before
#'               running the overdispersion function
#' @return a data.frame summarising overdispersion
#'
#' @note ensure column names in function call are in inverted commas
#'
#' @importFrom dplyr left_join
#' @importFrom stats aggregate pchisq qlogis fitted
#'
#' @export
overdispersion <- function(data, siteID = "siteID", visitID = "visitID", count = "count",
                           pass = "pass", id = "sampleID", species = "species", lifestage = "lifestage", 
                           largemodel) {

  # organise data set with standard names
  # varID are the supplied names
  # varID2 are the standard names we will use in the code 
  
  varID <- c(siteID, visitID, count, pass, id, species, lifestage)
  varID2 <- c("siteID", "visitID", "count", "pass", "id", "species", "lifestage")
  
  # check mandatory names are present
  if (!all(c(siteID, visitID, count, pass, id) %in% names(data)))
    stop("some mandatory variables are missing")
  
  # subset names in case species and lifestage are not in the dataframe
  ok <- varID %in% names(data)
  varID <- varID[ok]
  varID2 <- varID2[ok]
  
  # subset data
  data <- data[varID]
  
  # rename columns
  names(data) <- varID2
  
  # set visitID to a factor
  data$visitID <- as.factor(data$visitID)
  
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
  
  
  # get species:lifestage grouping variable for saturated model
  # variable has a separate value for each combination of species and lifestage
  # deals with the standard case where species have the same lifestages
  # deals also with the case e.g. where there are juvenile salmon and elver eels
  # fitting the saturate model will also be slightly more efficient
  
  group_var <- c("species", "lifestage")
  ok <- group_var %in% names(data)
  group_var <- group_var[ok]
  
  if (!any(ok)) 
    data$groupID <- "only_value" 
  else 
    data$groupID <- do.call("paste", data[group_var])
  
  
  ################
  #
  # SATURATED MODEL
  #
  ################
  
  # Produce a saturated model for each visit. '-1' prevents
  # the model estimate of intercept, and avoids problems if the intercept for the
  # first sample is poorly estimated. This results in a separate estimate of
  # intercept for every sample. Hessian = FALSE prevents the return of the hessian
  # matrix for returning standard errors, and speeds up model fitting.
  
  
  # Record and print the system time
  t1 <- Sys.time()
  message("Saturated model start time ", t1)
  
  saturated <- by(data, data$visitID, function(x) {
    
    # drop unused factor levels in data - e.g. id which gets passed into efp
    x <- droplevels(x)
    
    # check there is exactly one observation for each groupID and pass
    if (!all(table(x$groupID, x$pass) == 1)) 
      stop(
        "unbalanced data: some visits don't have exactly one observation ", 
        "for each combination of species, lifestage and pass"
      )	
    
    # create pass_reduce, a variable which combines the last two 
    # passes into a single category: 
    # e.g. if there are three passes, then pass_reduce is the same as pass12
    
    x$pass_reduce <- pmin(x$pass, x$maxpass-1)
    x$pass_reduce[x$pass_reduce == x$maxpass] <- x$maxpass - 1
    
    # combine this with groupID to get a saturated model formula
    x$groupID <- do.call("paste", x[c("groupID", "pass_reduce")])
    
    # fit saturated model
    m <- efp(
      count ~ -1 + groupID, 
      pass = pass, id = id, data = x, 
      hessian = FALSE, verbose = FALSE
    )
    
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
  etime <- t2 - t1
  message("Saturated model duration =", round(etime, 3), "s")

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

  # subset the dataframe to select only the required fields
  df <- subset(data, select = c("count", "pass", "id", "visitID"))

  # create a column which has the transformed p values from the large model
  df$p <- qlogis(fitted(largemodel))

  # split into a list of dataframes, one for each unique site visit
  v <- split(df, df$visitID)

  # create an empty dataframe to house the log likelihood and number of parameters
  # extracted from each loop of the model
  svis <- data.frame(llik = as.numeric(), params = as.numeric())

  # loop over the list of dataframes and estimate an intercept for each site visit
  for (i in 1:length(v)) {
    m <-
      efp(
        count ~ 1,
        pass = pass, id = id, data = v[[i]],
        offset = v$p, hessian = FALSE, verbose = FALSE
      )
    mod <-
      data.frame(
        llik = as.numeric(m$llik),
        params = as.numeric(length(m$coefficients))
      )

    # bind the model outputs to the empty dataframe and build up the dataframe during each loop
    svis <- rbind(svis, mod)
  }

  # store and print the system time and calculate the running time for the site-visit model
  t2 <- Sys.time()
  etime <- t2 - t1
  message("Site visit model duration = ", round(etime, 3), "s")

  # sum the log likelihoods and number of parameters from the models
  sitevisit <-
    data.frame(
      llik = sum(svis$llik),
      nparam = (sum(svis$params)) + length(largemodel$coefficients)
    )


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
  wk.anova <- rbind(saturated, sitevisit, largeout)
  row.names(wk.anova) <- c("saturated", "sitevisit", "large")

  # calculate the difference in the deviance and degrees of freedom between successive
  # model outputs and estimate the within visit (saturated to site-visit model) and between
  # visit (site-visit to large) overdispersion
  wk.anova$deviance <- c(NA, -2 * diff(wk.anova$llik))
  wk.anova$df <- c(NA, -diff(wk.anova$nparam))
  wk.anova$disp <- wk.anova$deviance / wk.anova$df
  wk.anova$p <- pchisq(wk.anova$deviance, wk.anova$df, lower.tail = FALSE)

  # return the overdispersion estimates and p values
  # the values for "disp" are then included in the adjusted BIC function
  out <-
    data.frame(
      wk.anova[c("llik", "nparam", "deviance", "df", "disp", "p")]
    )

  # the end result of the function is to return the overdispersion output
  return(out)
}
