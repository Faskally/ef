
#' Estimating overdispersion
#'
#' Complete function for returning overdispersion estimates
#'
#'
#' @param data dataframe containing EF data
#' @param siteID site name or unique ID
#' @param visitID a number identifying each unique visit
#' @param count count of fish (defaults to "count")
#' @param pass the EF pass number (defaults to "pass")
#' @param lifestage the lifestage (defaults to "lifestage")
#' @param pass12 categorical variable with 2 levels where the 1st pass
#'            and subsequent passes are treated separately (defaults to "pass12")
#' @param id sample ID
#' @param largemodel the large model that captures most of the systematic
#'                variation in the data - this is specified before
#'               running the overdispersion function
#' @return a data.frame summarising overdispersion
#'
#' @note ensure column names in function call are in inverted commas
#'
#' @export
overdispersion <- function(data, siteID, visitID, count = "count",
                           pass = "pass", lifestage = "lifestage",
                           pass12 = "pass12", id, largemodel){

  # Subset the dataframe to select only the columns that will be used in the function
  varID <- c(siteID, visitID, count, pass, lifestage, pass12, id)

  # Create a vector of names to standardise the naming convention in the function
  varID2 <- c("siteID", "visitID", "count", "pass", "lifestage", "pass12", "id")

  # subset data
  data <- data[varID]

  # rename columns
  names(data) <- varID2

  # set visitID to a factor
  data$visitID <- as.factor(data$visitID)

  # return the maximum number of EF passes for each visit
  aggregatepass<-aggregate(pass~visitID, FUN=max, data=data)

  # rename the columns
  names(aggregatepass)<-c("visitID", "maxpass")

  # set the visitID to be a factor
  aggregatepass$visitID <- as.factor(aggregatepass$visitID)

  # join the data and the maximum number of passes together. Use
  # left_join to maintain the dataframe order
  data<-left_join(data,aggregatepass, by="visitID")

  # split the dataframe by the maximum number of passes to get a list
  # of dataframes where the total number of passes is either 2, 3 or 4
    x <- split(data[c("pass", "id", "count", "maxpass", "pass12",
                    "lifestage", "visitID")],
             data$maxpass)

  # add warnings if the list of dataframes includes site visits with
  # less than 2 passes or more than 4 passes
  if(min(as.numeric(names(x))) <2 | (max(as.numeric(names(x))) > 4))
      stop('pass number not between 2 and 4')


    ################
    #
    # SATURATED MODEL
    #
    ################

  # create an empty dataframe to store the log likelihood and number
  # of parameters for each model
  Saturated<-data.frame(llik=as.numeric(),params=as.numeric())

  # Produce a saturated model for 2, 3, and 4-pass data if present. '-1' prevents
  # the model estimate of intercept, and avoids problems if the intercept for the
  # first sample is poorly estimated. This results in a separate estimate of
  # intercept for every sample. Hessian = FALSE prevents the return of the hessian
  # matrix for returning standard errors, and speeds up model fitting.


  # NOTE: This next step repeats for 2, 3 and 4-pass data

  # go to list of dataframes and check to see if 4 pass data recorded as "4" in list
  if(max(as.numeric(names(x))) == 4){

    # Create a categorical variable for pass123 which
    # treats passes 1 and 2 separately from 3 and 4
    x[["4"]] <- within(x[["4"]], {
      pass123 <- factor(pass)
      levels(pass123)[levels(pass123) == "4"] <- "3"
    })

    # Record and print the system time
    t1=Sys.time()
    cat(paste("4-pass model start time",t1), "\n")
    cat("Modelling 4-pass data...", "\n")

     # create a dataframe from the 4-pass data
    d <- x[["4"]]

    # drop levels to remove unused factor levels
    d <- droplevels(d)

    # split the dataframe by visitID to get a list of dataframes
    # for each individual visit to speed up model fitting
    v <- split(d[c("pass", "id", "count", "maxpass", "pass123", "lifestage")],
               d$visitID)
    # create a blank dataframe to sore model outputs
    out <- data.frame(llik=as.numeric(),params=as.numeric())

    # loop over each dataframe in the list and extract the log likelihood and
    # number of parameters (length of coefficients) from the saturated model for
    # each site visit
    for(i in 1:length(v)){
      mod <- data.frame(llik=as.numeric(efp(count ~ -1 + lifestage:pass123, pass = pass, id = id,
                                        data = v[[i]], hessian = FALSE, verbose = FALSE)$llik),
                        params=as.numeric(length(efp(count ~ -1 + lifestage:pass123, pass = pass, id = id,
                                          data = v[[i]], hessian = FALSE, verbose = FALSE)$coefficients)))

      # rbind the model outputs with the blank 'out' dataframe - this will
      # build up a list of model outputs during the looping process
      out <- rbind(out, mod)
    }

    # sum all the log likelihoods and number of parameters from the saturated
    # models for 4-pass data
    sat4pass <- data.frame(
      llik = sum(out$llik),
      nparam = sum(out$params))

    # rbind the empty dataframe create outside the saturated function with the
    # summed 4-pass saturated model data
    Saturated <- rbind(Saturated, sat4pass)

    # store and print the time that it took for the model to run
    t2=Sys.time()
    cat(paste("4-pass model end time", t2), "\n")
    etime=t2-t1
    cat(paste("4-pass model duration =", etime), "\n")
  }

  if (3 %in% (as.numeric(names(x)))){
    t1=Sys.time()
    cat(paste("3-pass model start time",t1), "\n")
    cat("Modelling 3-pass data...", "\n")

    d <- x[["3"]]
    d <- droplevels(d)
    v <- split(d[c("pass", "id", "count", "maxpass", "pass12", "lifestage")],
               d$visitID)
    out <- data.frame(llik=as.numeric(),params=as.numeric())
    for(i in 1:length(v)){
      mod <- data.frame(llik=as.numeric(efp(count ~ -1 + lifestage:pass12, pass = pass, id = id,
                                        data = v[[i]], hessian = FALSE, verbose = FALSE)$llik),
                        params=as.numeric(length(efp(count ~ -1 + lifestage:pass12, pass = pass, id = id,
                                          data = v[[i]], hessian = FALSE, verbose = FALSE)$coefficients)))
      out <- rbind(out, mod)
    }

    sat3pass <- data.frame(
      llik = sum(out$llik),
      nparam = sum(out$params))

    Saturated <- rbind(Saturated, sat3pass)
    t2=Sys.time()
    cat(paste("3-pass model end time", t2), "\n")
    etime=t2-t1
    cat(paste("3-pass model duration =", etime), "\n")
  }

  if(min(as.numeric(names(x))) == 2){
    t1=Sys.time()
    cat(paste("2-pass model start time",t1), "\n")
    cat("Modelling 2-pass data...", "\n")

    d <- x[["2"]]
    d <- droplevels(d)
    v <- split(d[c("pass", "id", "count", "maxpass", "lifestage")],
               d$visitID)
    out <- data.frame(llik=as.numeric(),params=as.numeric())
    for(i in 1:length(v)){
      mod <- data.frame(llik=as.numeric(efp(count ~ -1 + lifestage, pass = pass, id = id,
                                        data = v[[i]], hessian = FALSE, verbose = FALSE)$llik),
                        params=as.numeric(length(efp(count ~ -1 + lifestage, pass = pass, id = id,
                                          data = v[[i]], hessian = FALSE, verbose = FALSE)$coefficients)))
      out <- rbind(out, mod)
    }

    sat2pass <- data.frame(
      llik = sum(out$llik),
      nparam = sum(out$params))

    Saturated <- rbind(Saturated, sat2pass)

    t2=Sys.time()
    cat(paste("2-pass model end time", t2), "\n")
    etime=t2-t1
    cat(paste("2-pass model duration =", etime), "\n")
  }


  # Create dataframe with the sum of the logliks and nparams from the 2, 3 and
  # 4 pass data if you have it
    saturated <- data.frame(
        llik = sum(Saturated$llik),
        nparam = sum(Saturated$nparam))


  ################
  #
  # SITE-VISIT MODEL
  #
  ################


  # The site-visit model has an extra parameter for each site-visit. Due to the
  # number of  parameters, this model is fitted by conditioning on the fitted
  # values from the large model and estimating the site-visit effect for each
  # site-visit in turn.

  # The site-visit model includes estimates of all of the parameters from the full model
  # as the offset, thus the additional variability due to the inclusion of the site-visits
  # gives you an estimate of the between visit overdispersion.


  # Record and print system time
  t1=Sys.time()
  cat(paste("Site visit model start time",t1), "\n")
  cat("Running Site Visit model...", "\n")

  # subset the dataframe to select only the required fields
  df <- subset(data, select=c("count", "pass", "id", "visitID"))

  # create a column which has the transformed p values from the large model
  df$p <- qlogis(largemodel$p)

  # split into a list of dataframes, one for each unique site visit
  v <- split(df, df$visitID)

  # create an empty dataframe to house the log likelihood and number of parameters
  # extracted from each loop of the model
  svis <- data.frame(llik=as.numeric(),params=as.numeric())

  # loop over the list of dataframes and estimate an intercept for each site visit
  for(i in 1:length(v)){
    mod <- data.frame(llik=as.numeric(efp(count ~ 1, pass = pass, id = id, data = v[[i]],
                                      offset = v$p, hessian = FALSE, verbose = FALSE)$llik),
                      params=as.numeric(length(efp(count ~ 1, pass = pass, id = id, data = v[[i]],
                                        offset = v$p, hessian = FALSE, verbose = FALSE)$coefficients)))

    # bind the model outputs to the empty dataframe and build up the dataframe during each loop
    svis <- rbind(svis, mod)
  }

  # store and print the system time and calculate the running time for the site-visit model
  t2=Sys.time()
  cat(paste("Site visit model end time", t2), "\n")
  etime=t2-t1
  cat(paste("Site visit model duration =", etime), "\n")

  # sum the log likelihoods and number of parameters from the models
  sitevisit <- data.frame(
        llik = sum(svis$llik),
        nparam = (sum(svis$params)) + length(largemodel$coefficients))


  ################
  #
  # ESTIMATING OVERDISPERSION
  #
  ################

  # extract the log-likelihood and number of parameters from the large model
  largeout <- data.frame(
        llik = largemodel$llik,
        nparam = length(largemodel$coefficients))

  # bind together the log-likelihoods and number of parameters from each
  # of the three modelling processes
  wk.anova <- rbind(saturated, sitevisit, largeout)
  row.names(wk.anova) = c("saturated", "sitevisit", "large")

  # calculate the difference in the deviance and degrees of freedom between successive
  # model outputs and estimate the within visit (saturated to site-visit model) and between
  # visit (site-visit to large) overdispersion
  wk.anova <- within(wk.anova, {
    deviance <- c(NA, - 2 * diff(llik))
    df <- c(NA, - diff(nparam))
    disp <- deviance / df
    p <- pchisq(deviance, df, lower.tail = FALSE)
  })

  # return the overdispersion estimates and p values
  # the values for "disp" are then included in the adjusted BIC function
    overdispersion.output <- data.frame(
    wk.anova[c("llik", "nparam", "deviance", "df", "disp", "p")]
  )
  # the end result of the function is to return the overdispersion output
  return(overdispersion.output)

}




#' Adjusted BIC function
#'
#' Complete function for returning overdispersion estimates
#'
#'
#' @param model a fitted ef model
#' @param data the data used to fit the model
#' @param overdispersion.output output form the function \link{overdispersion}
#' @return the adjusted BIC
#'
#' @note When calling the function, you need to specify the data source for
#'   the model so that the number of site visits can be determined.
#'   You also need to specify the output from the overdispersion model
#'   to get the measure
#'
#' @export
BICadj <- function(model, data, overdispersion.output) {

    # The adjusted BIC value is calculated in a similar way to BIC, however the
  # loglikelihood of the model is divided by the measure of between sample
  # overdispersion, which makes it more stringent when adding or removing terms
  # from the model.

  BICval <- data.frame(

    Adjusted.BIC =

      -((2*model$llik) / overdispersion.output[3,5])
    + ((log(length(unique(paste(data$siteID, data$date)))))
       * (length(model$coefficients))),

    DF =

      # The degrees of freedom are the number of coefficients minus the intercept
      length((model$coefficients)-1))

  # This displays the adjusted BIC value along with the degrees of freedom
  # NOTE: nothing is stored here, only displayed on-screen
  BICval[c("Adjusted.BIC", "DF")]

}


#################################################################################################
#################################################################################################


