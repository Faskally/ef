# edit history

# 21/12/16 adjust for gamm models
# 22/12/16 correct bug when picking up AIC for gamm4 models
# 01/01/17 by argument allows the by terms to be expanded when the by variable is continuous 
# 01/01/17 essential argument allows terms to always be present
# 14/02/17 fix bug when smoothTerms is missing; also deal with random effects in the formula (lme4)
# 29/03/17 allows more general likelihood calculation
# 10/09/19 allows for 2D smoothers - see notes for explanation of how to use
# 10/09/19 fitModels = FALSE returns the candidate models without fitting them 
# 11/09/19 catch for errors in model fits
# 12/09/19 retainModels = FALSE only returns the summary table (useful when models are large)

# v6_9
# 22/06/20 drop ... from argument list to force criterion to pick up information from model object

drop1.efp <- function(
  model, scope, smoothTerms, notLinear, essential, criterion, cluster = NULL, 
  byExpand = FALSE, calcLogLik, fitModels = TRUE, retainModels = TRUE) {

  # scope is the maximal model which, for simplicity, is expressed only using linear terms.
  # So, if 
  #   smoothTerms = list(
  #     altitude = "s(altitude, k = 3)", 
  #     dayOfYear = "s(dayOfYear, k = 3)",
  #     network = "s(network, k = 4)", 
  #     hydroArea = "s(hydroArea, k = 4, bs = \"tp\")", 
  #     longitude = t2(longitude, latitude)
  #   )
  # then 
  #   scope = pass * lifestage * (altitude + dayOfYear + longitude)
  # denotes a maximal model
  #   pass + lifestage + pass:lifestage + s(altitude, k = 3, by = pass:lifestage) + 
  #      s(dayOfYear, k = 3, by = pass:lifestage) + t2(longitude, latitude, by = pass:lifestage)
  # NB Don't use : in scope and don't use 1 unless it is the null model (which is pretty pointless)
  
  # NB Still need to sort out code when by is used with more that one variable
  # Currently, the term s(x, by = y:z) is actually fitted as s(x, by = y) + s(x, by = z) + s(x, by = y:z)
  # If you get models like this turning up in the model selection, interpret them with great care!
  
  # smoothTerms gives the form of the smoother for each term in the model
  # NB for a 2 or 3 dimensional smoother, just give the first term - these smoothers must also 
  #   be identified in the notLinear argument
  
  # notLinear specifies the smoothers which are not reduced to a linear term before being dropped 
  # e.g. notLinear = "hydroArea" would consider dropping the hydroArea term completely, rather than
  # introducing hydroArea as a factor (which doesn't make sense because that is a more complicated model)
  
  # essential specifies terms that must be in every model
  # these must not be the same as any terms in the scope (or even related to them in any way)
  # leave any offset out of this - dealt with separately
  
  # criterion (default = AIC) is a function which allows the user to define the model selection criterion
  
  # cluster is the identity of a cluster to allow parallel processing
  # drop1.efp and other terms (model, the data set used, scope and smoothTerms) must previously be 
  # exported to cluster using clusterExport
  
  # byExpand = FALSE (default) is designed for use when the by variables are factors - fits a separate
  # smoother for each level of the factor, but with no main effect (of the covariate)
  # byExpand = TRUE is for when there is a single by variable and it is continuous - fits both the  
  # main effect of the covariate (as a smoother) and the product of the by variable and a different smooth
  # of the covariate
  
  # calcLogLik (default = logLik) allows for a more general log likelihood calculation, for example, in 
  # gam models with random effects, where you would want to use the marginal (integrated) likelihood in 
  # object$gcv.ubre.  Only argument is model so 
  # calcLogLik = function(model) -model$gcv.ubre is the obvious choice
  
  # fitModels = TRUE (default) fits all the candidate model
  # fitModels = FALSE just returns the list of candidate models without fitting - useful if you want
  # to just check which terms come in and drop out
  
  # retainModels = TRUE (default) returns the updated model
  # retainModels = FALSE just returns the updated model formula - saves memory
 
 
   require(pbapply)
  
  
  # utility functions - throughout id is the identity of a covariate in smoothTerms
  # and rhs is a character vector of model terms

  # converts rhs to single formula-like string
  
  pasteRHS <- function(rhs) paste(rhs, collapse = " + ")

  # turns rhs into a proper formula
    
  rhsToFormula <- function(rhs) formula(paste("~", pasteRHS(rhs)))
  
  # expands rhs to show all main effects and interactions
  
  expandTerms <- function(x) {
    if (is.character(x))
      x <- rhsToFormula(x)
    rhsTerms <- terms(x)
    rhs <- attr(rhsTerms, "term.labels")
    if ("offset" %in% names(attributes(rhsTerms))) {
      offset <- as.character(attr(rhsTerms, "variables")[attr(rhsTerms, "offset") + 1])
      rhs <- c(rhs, offset)
    }
    if (length(rhs)) gsub(" ", "", rhs, fixed = TRUE) else as.character(1)
  }


  # position of smooth terms (if any) in rhs involve id
  
  smoothLocate <- function(id, rhs) {
    pos <- sapply(rhs, function(str) {
      # deal with linear terms and offset
      wk <- strsplit(str, "(", fixed = TRUE)[[1]]
      if (length(wk) == 1 | wk[1] == "offset") return(FALSE)
      # now just have smooth terms
      id %in% eval(parse(text = str))$term
    })
    seq_along(rhs)[pos]
  }

  isSmooth <- function(id, rhs) length(smoothLocate(id, rhs)) > 0
  
  # and likewise for linear terms
  
  linearLocate <- function(id, rhs) {
    if (! id %in% rhs) 
      return(numeric(0))
    rhsTerms <- terms(rhsToFormula(rhs))
    out <- attr(rhsTerms, "factors")
    colnames(out) <- gsub(" ", "", colnames(out))
    if ("offset" %in% names(attributes(rhsTerms))) {
      out <- cbind(out, 0)
      colnames(out)[ncol(out)] <- as.character(attr(rhsTerms, "variables")[attr(rhsTerms, "offset") + 1])
    }
    seq_along(rhs)[out[id, rhs] == 1]
  }

  isLinear <- function(id, rhs) length(linearLocate(id, rhs)) > 0
  
      
  # turns all the smooth terms involving id into linear terms
  # needed when considering whether linear terms are adequate or 
  # to make use of drop.scope and add.scope functionality
  # 
  # take out$term[1] as a work around for 2D smooths
    
  smoothToLinear <- function(id, rhs) {
    pos <- smoothLocate(id, rhs)
    rhs[pos] <- sapply(rhs[pos], function(x) {
      out <- eval(parse(text = x))
      if (out$by == "NA")
        return(out$term[1])
      out <- paste(out$term[1], out$by, sep = "*")
      gsub(":", "*", out, fixed = TRUE)
    })
    expandTerms(rhs)
  }
  
  # and the other way round
  
  linearToSmooth <- function(id, rhs, smooth) {
    pos <- linearLocate(id, rhs)
    if (length(pos) == 0) 
      return(rhs)
    drop <- rhs[pos]
    ordTerms <- sapply(strsplit(drop, ":"), length)
    if (max(ordTerms) == 1) 
      add <- smooth
    else {
      add.extra <- sapply(drop[ordTerms > 1], function(x) {
        x <- sub(paste0(id, ":"), "", x)
        x <- sub(paste0(":", id), "", x)
        x <- paste0(",by=", x)
        y <- substr(smooth, start = 0, stop = nchar(smooth) - 1)
        paste0(y, x, ")")
      })
      if (byExpand) {
        if (any(ordTerms) >= 3)
          stop("can only have one covariate interacting with a smoother when byExpand = TRUE")
        add <- c(smooth, add.extra)
      }
      else
        add <- add.extra
    }
    rhs <- rhs[-pos]
    rhs <- append(rhs, add, pos[1]-1)
    unname(rhs)
  }
 

  # sort out model selection criterion
  
  if (!missing(criterion))
    critID <- deparse(substitute(criterion))
  else {
    critID <- "AIC"
    if ("gamm4" %in% class(model)) 
      criterion <- function(model) extractAIC(model$mer)[2] 
    else if ("gamm" %in% class(model)) 
      criterion <- function(model) AIC(model$lme) 
    else 
      criterion <- AIC
  }
  

  # and likelihood calculation
  
  if (missing(calcLogLik))
    calcLogLik <- logLik


  # get model call and formula and set up output structure, with model fit information
  # from starting model

  call <- getCall(model)
  formula <- as.formula(call$formula)
  
  summary <- data.frame(direction = "none", inTerms = "", outTerms = "", 
                        logLik = calcLogLik(model), value = criterion(model))
  names(summary)[names(summary) == "value"] <- critID
  
  
  # expand model formula to get all terms, including offset if there is one, 
  # remove any spaces (to simplify specification) 

  rhs <- expandTerms(formula)


  # get lhs of model formula for later model calls
  
  lhs <- paste(as.character(formula)[c(2, 1)], collapse = " ")


  # look at smoothed terms: remove spaces, work out whether present as smooths in starting model
  
  if (!missing(smoothTerms)) {
    smoothTerms <- lapply(smoothTerms, gsub, pattern = " ", replacement = "", fixed = TRUE)
    smoothID <- names(smoothTerms)
  } else
    smoothID <- character(0)

  
  # establish which terms could be linearised
  
  if (missing(notLinear))
    notLinear <- character(0)
  stopifnot(notLinear %in% smoothID)


  
  if (length(smoothID) == 0)
    canLinearise <- character(0)
  else {
    canLinearise <- sapply(smoothID, isSmooth, rhs = rhs)
    canLinearise <- smoothID[canLinearise]
    canLinearise <- setdiff(canLinearise, notLinear)
  }  

  
  # establish which terms are essential and hence obtain  a set of candidate terms to add or drop
  # when random terms are in formula, need to deal with leading and trailing brackets that
  # get dropped in rhs
  
  if (!missing(essential)) { 
    essentialRT <- essential <- gsub(essential, pattern = " ", replacement = "", fixed = TRUE)
    isRT <- substring(essential, 1, 1) == "(" & substring(essential, nchar(essential)) == ")"
    essentialRT[isRT] <- substring(essentialRT[isRT], 2, nchar(essentialRT[isRT]) - 1)
    stopifnot(essentialRT %in% rhs)
    rhsLin <- setdiff(rhs, essentialRT)
  } else {
    essential <- character(0)
    rhsLin <- rhs
  }
  
  
  # work out which terms to add or drop - to use drop.scope and add.scope functionality
  # first make rhs linear

  for (id in smoothID) {
    if (isSmooth(id, rhs)) 
      rhsLin <- smoothToLinear(id, rhsLin)
  }


  addTerms <- add.scope(rhsToFormula(rhsLin), scope)

  dropTerms <- drop.scope(rhsToFormula(rhsLin))
  dropTerms <- dropTerms[! dropTerms %in% canLinearise]
    
  addRHS <- lapply(addTerms, function(x) {
    if (length(rhsLin) == 1 && rhsLin == "1") 
      x 
    else 
      expandTerms(c(rhsLin, x))
  })

  dropRHS <- lapply(dropTerms, function(x) {
    out <- setdiff(rhsLin, x)
    if (length(out)) out else "1"
  })
  
  # convert rhs for smooth terms back

  for (id in smoothID) {
    if (isSmooth(id, rhs) || id %in% notLinear) {
      addRHS <- lapply(addRHS, linearToSmooth, id = id, smooth = smoothTerms[[id]])
      dropRHS <- lapply(dropRHS, linearToSmooth, id = id, smooth = smoothTerms[[id]])
    }
  }


  # add in conversion from linear to smooth and demotion from smooth to linear (except if notLinear)
  
  for (id in setdiff(smoothID, notLinear)) {
    if (isSmooth(id, rhs))
      dropRHS <- c(dropRHS, list(smoothToLinear(id, rhs)))
    else if (isLinear(id, rhs))
      addRHS <- c(addRHS, list(linearToSmooth(id, rhs, smooth = smoothTerms[[id]])))
  }  
  

  # add back in essential terms
  
  addRHS <- lapply(addRHS, append, essential)
  dropRHS <- lapply(dropRHS, append, essential)
  
  
  # combine to single string
  
  addRHS <- sapply(addRHS, pasteRHS)
  dropRHS <- sapply(dropRHS, pasteRHS)

  newModels <- rbind(
    if (length(addRHS)) data.frame(direction = "add", rhs = addRHS, stringsAsFactors = FALSE), 
    if (length(dropRHS)) data.frame(direction = "drop", rhs = dropRHS, stringsAsFactors = FALSE)
  )


  newModels$inTerms <- sapply(newModels$rhs, function(x) {
    out <- setdiff(expandTerms(x), rhs)
    pasteRHS(out)
  })

  newModels$outTerms <- sapply(newModels$rhs, function(x) {
    out <- setdiff(rhs, expandTerms(x))
    if (length(out) == 1 && out == "1") "" else pasteRHS(out)
  })

  newModels <- newModels[c("direction", "inTerms", "outTerms", "rhs")]
 
  if (!fitModels)
    return(newModels) 

  
  # now refit model with new model formulae and get AIC
  
  refit <- function(x) {
    formula <- paste(c(lhs, x), collapse = " ")
    call$formula <- as.formula(formula)
    fit <- try(eval(call, parent.frame()))
    error <- class(fit) == "try-error"
    out <- list(formula = formula, error = error)
    out[["logLik"]] <- if (!error) calcLogLik(fit) else NA
    out[[critID]] <- if (!error) criterion(fit) else NA
    if (retainModels)
      out[["fit"]] <- fit
    out
  }
  
  fits <- pblapply(newModels$rhs, refit, cl = cluster)
  
  id <- c("formula", "error", "logLik", critID)
  newModels[id] <- lapply(id, function(x) sapply(fits, "[[", x))

  error <- any(newModels$error)

  
  # tidy up output
  
  pos <- match(c("rhs", "formula", "error"), names(newModels))
  summary <- rbind(summary, newModels[-pos])
  summary <- summary[order(summary[[critID]]), ]
  
  out <- list(summary = summary, error = error)
  
  if (out$error) {
    warning("there were errors in some model fits")
    return(out)
  }
  

  # select optimal model
  
  if (min(newModels[[critID]]) < criterion(model)) {
    update = TRUE
    newID <- which.min(newModels[[critID]])
    nextFormula <- newModels$formula[newID]
    if (retainModels)
      nextModel <- fits[[newID]][["fit"]]
  }
  else {
    update = FALSE 
    nextFormula <- model$call$formula
    nextFormula <- as.character(nextFormula)
    nextFormula <- paste(nextFormula[c(2, 3)], collapse = " ~ ")
    if (retainModels)
      nextModel <- model
  }

  out[["update"]] <- update
  out[["nextFormula"]] <- nextFormula
  if (retainModels)
    out[["nextModel"]] <- nextModel
  out
}

