# Penalised smooth terms are fitted as random effects. A smooth with
# coefficients b, design matrix X and penalty matrix S is written as
#
#   b = N g + E u,    u ~ N(0, sigma^2 * S[keep, keep]^-1)
#
# where the columns of N span the null space of the penalty, and E selects
# the coefficients that remain (keep) when one coefficient is removed for
# each dimension of the null space. The parameters g are fixed effects and u
# are random effects with a proper distribution. For a random walk this is
# the usual formulation in which the starting values are fixed unknown
# parameters. Nothing is done to X or S other than removing rows and
# columns, so if they are sparse they stay sparse.
#
# The marginal likelihood depends (slightly) on which coefficients are
# removed, as the fixed effects are estimated by maximum likelihood rather
# than integrated out, so the choice has to be reproducible.

# orthonormal basis for the null space of a positive semi-definite matrix
null_space <- function(S, tol = 1e-8) {
  k <- ncol(S)
  if (k <= 1500) {
    es <- eigen(as.matrix(S), symmetric = TRUE)
    return(es$vectors[, es$values < tol * max(es$values), drop = FALSE])
  }

  # too big for a dense decomposition: use subspace iteration on the
  # inverse of a slightly shifted S, which only needs sparse solves
  scale <- max(abs(diag(S)))
  ch <- Cholesky(forceSymmetric(S + Diagonal(k, tol * scale)))
  q <- 8
  repeat {
    V <- qr.Q(qr(matrix(stats::rnorm(k * q), k, q)))
    for (i in 1:10) {
      V <- qr.Q(qr(as.matrix(solve(ch, V))))
    }
    es <- eigen(as.matrix(crossprod(V, S %*% V)), symmetric = TRUE)
    is_null <- es$values < tol * scale
    if (!all(is_null) || q >= k) break
    q <- min(2 * q, k)
  }
  V %*% es$vectors[, is_null, drop = FALSE]
}

# split a smooth into fixed and random effects
#   X: design matrix, S: list of penalty matrices
#   centred: if TRUE the constant function is removed from the smooth,
#     which requires that it is in the null space of the penalty. If it is
#     not NULL is returned.
smooth_block <- function(X, S, centred = FALSE) {
  X <- as(as(X, "generalMatrix"), "CsparseMatrix")
  S <- lapply(S, function(s) forceSymmetric(as(s, "CsparseMatrix")))
  k <- ncol(X)

  N <- null_space(Reduce(`+`, lapply(S, function(s) s / max(abs(s@x)))))
  Ng <- N
  const <- NULL

  if (centred) {
    # find the coefficients which give the constant function
    if (ncol(N) == 0) return(NULL)
    M <- as.matrix(X %*% N)
    a <- qr.coef(qr(M), rep(1, nrow(M)))
    a[is.na(a)] <- 0
    if (max(abs(M %*% a - 1)) > 1e-8) return(NULL)
    const <- drop(N %*% a)
    Ng <- N %*% qr.Q(qr(matrix(a)), complete = TRUE)[, -1, drop = FALSE]
  }

  # the coefficients to remove, one for each dimension of the null space:
  # these are the first coefficients which identify the null space, so for
  # a random walk they are the starting values
  drop_coef <- integer(0)
  if (ncol(N)) {
    drop_coef <- sort(qr(t(N))$pivot[seq_len(ncol(N))])
  }
  keep <- setdiff(seq_len(k), drop_coef)
  S <- lapply(S, function(s) s[keep, keep, drop = FALSE])

  out <-
    list(
      Xf = as.matrix(X %*% Ng),
      Zu = X[, keep, drop = FALSE],
      # coefficients of the smooth are B %*% c(g, u)
      B = cbind(Matrix(Ng, sparse = TRUE), Diagonal(k)[, keep, drop = FALSE]),
      const = const,
      keep = keep
    )

  if (length(S) == 1 && isDiagonal(S[[1]]) && all(diag(S[[1]]) == 1)) {
    # independent random effects
    out$S <- NULL
  } else {
    if (length(S) == 1) {
      out$logdet <- as.numeric(determinant(S[[1]], logarithm = TRUE)$modulus)
    }
    out$S <- lapply(S, function(s) as(s, "generalMatrix"))
  }
  out$npen <- length(S)

  out
}

# the matrix Z which mgcv uses to apply the identifiability constraint to a
# smooth, so that the constrained design matrix is X %*% Z, and the
# direction that is removed by the constraint
constraint_basis <- function(sm, k) {
  qrc <- attr(sm, "qrc")
  indi <- attr(sm, "indi")
  if (is.null(indi)) indi <- seq_len(k)
  nx <- length(indi)
  Q <- qr.Q(qrc, complete = TRUE)

  Z <- diag(k)
  Z[indi, indi] <- 0
  Z[indi, indi[-nx]] <- Q[, -1, drop = FALSE]
  q1 <- numeric(k)
  q1[indi] <- Q[, 1]

  list(Z = Z[, -indi[nx], drop = FALSE], q1 = q1)
}

# set up the penalised terms in a model
#   Gsetup: the output of gam(..., fit = FALSE)
#   sparse: if TRUE the smooths are constructed without the mgcv
#     identifiability constraints, which keeps sparse bases and penalties
#     sparse. If FALSE, or if this is not possible for a smooth, the
#     constrained (and dense) mgcv basis and penalty are used.
#' @importFrom mgcv smoothCon interpret.gam
smooth_setup <- function(Gsetup, formula, sparse = TRUE) {
  G <- Gsetup$X
  intercept <- which(Gsetup$term.names == "(Intercept)")

  penalised <- which(vapply(Gsetup$smooth, function(x) length(x$S) > 0, TRUE))
  if (!length(penalised)) return(list())

  # the smooths without identifiability constraints
  full <- NULL
  if (sparse) {
    full <-
      try(
        unlist(
          lapply(
            interpret.gam(formula)$smooth.spec, smoothCon,
            data = Gsetup$mf, knots = NULL, absorb.cons = FALSE,
            scale.penalty = FALSE
          ),
          recursive = FALSE
        ),
        silent = TRUE
      )
    labels <- vapply(Gsetup$smooth, function(x) x$label, "")
    if (inherits(full, "try-error") ||
        !identical(vapply(full, function(x) x$label, ""), labels)) {
      full <- NULL
    }
  }

  lapply(penalised, function(i) {
    sm <- Gsetup$smooth[[i]]
    cols <- sm$first.para:sm$last.para
    ncons <- attr(sm, "nCons")
    if (is.null(ncons)) ncons <- 0

    block <- NULL
    if (!is.null(full) && ncons == 0 && ncol(full[[i]]$X) == length(cols)) {
      block <- smooth_block(full[[i]]$X, full[[i]]$S)
      block$TG <- block$B
    } else if (!is.null(full) && ncons == 1 && length(intercept) == 1 &&
               inherits(attr(sm, "qrc"), "qr")) {
      X <- full[[i]]$X
      con <- constraint_basis(sm, ncol(X))
      if (ncol(con$Z) == length(cols) &&
          max(abs(X %*% con$Z - G[, cols])) < 1e-8 * max(abs(G[, cols]))) {
        block <- smooth_block(X, full[[i]]$S, centred = TRUE)
      }
      if (!is.null(block)) {
        # move the constant part of the smooth to the intercept, so that the
        # coefficients satisfy the mgcv constraint
        shift <- drop(con$q1 %*% block$B) / sum(con$q1 * block$const)
        block$TG <- crossprod(con$Z, as.matrix(block$B) - outer(block$const, shift))
        block$Tint <- shift
      }
    }

    if (is.null(block)) {
      # use the constrained mgcv basis and penalties, undoing the scaling of
      # the penalties
      S <- Map(function(s, scale) s * scale, sm$S, sm$S.scale)
      block <- smooth_block(G[, cols, drop = FALSE], S)
      block$TG <- block$B
    }

    block$cols <- cols
    block$label <- sm$label
    block$sparse <- !is.null(block$Tint) || (!is.null(full) && ncons == 0)

    if (ncol(block$Xf)) {
      colnames(block$Xf) <- paste0(sm$label, ".null", seq_len(ncol(block$Xf)))
    }
    if (length(block$keep) == length(cols) && is.null(block$Tint)) {
      colnames(block$Zu) <- Gsetup$term.names[cols]
    } else {
      colnames(block$Zu) <- paste0(sm$label, ".u", block$keep)
    }

    block
  })
}
