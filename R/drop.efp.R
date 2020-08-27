#' @export
drop1.efp <- function(object, scope, scale = 1, test = AIC, ...) {
  #check_exact(object)
  x <- model.matrix(object)
  #offset <- model.offset(model.frame(object))
  iswt <- !is.null(wt <- object$weights)
  n <- nrow(x)
  #asgn <- attr(x, "assign")
  terms <- terms(formula(object))
  tl <- attr(terms, "term.labels")
  if (missing(scope))
    scope <- drop.scope(formula(object))
  else {
    if (!is.character(scope))
      scope <- attr(terms(update.formula(object, scope)),
                    "term.labels")
    if (!all(match(scope, tl, 0L) > 0L))
      stop("scope is not a subset of term labels")
  }
  ndrop <- match(scope, tl)
  ns <- length(scope)
  ## keep base model AICs etc
  #rdf <- object$df.residual
  #chisq <- deviance.lm(object)
  ## containers for model summaries
  dfs <- numeric(ns)
  RSS <- numeric(ns)

  # set up fit call
  fcall <- object$call

  for (i in seq_len(ns)) {
    # construct one model
    tts <- as.formula(paste("~ 1 + ", paste(scope[-i], collapse = " + "), sep = " "))
    form <- update.formula(object, tts)

        z <-
    dfs[i] <- z$rank
    oldClass(z) <- "lm"
    RSS[i] <- deviance(z)
  }
  scope <- c("<none>", scope)
  dfs <- c(object$rank, dfs)
}
