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
null_space <- function(S) {
  es <- eigen(as.matrix(S), symmetric = TRUE)
  tol <- max(es$values) * .Machine$double.eps^0.75
  es$vectors[, es$values < tol, drop = FALSE]
}

# split a smooth into fixed and random effects
#   X: design matrix, S: list of penalty matrices
#   centred: if TRUE the constant function is removed from the smooth,
#     which requires that it is in the null space of the penalty. If it is
#     not NULL is returned.
#   N: a basis for the null space of the penalty, if known
smooth_block <- function(X, S, centred = FALSE, N = NULL) {
  X <- as(as(X, "generalMatrix"), "CsparseMatrix")
  S <- lapply(S, function(s) forceSymmetric(as(s, "CsparseMatrix")))
  k <- ncol(X)

  if (is.null(N)) {
    Ssum <- Reduce(`+`, lapply(S, function(s) s / max(abs(s@x))))
    if (isDiagonal(Ssum) && all(diag(Ssum) > 0)) {
      # everything is penalised, as for independent random effects
      N <- matrix(0, k, 0)
    } else {
      N <- null_space(Ssum)
    }
  }
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

# the design matrix and penalty of a Markov random field smooth, in the
# order of the levels in the penalty supplied by the user. For a random walk
# this is the order of the walk, mgcv sorts the levels alphabetically.
mrf_parts <- function(sm, mf) {
  S <- sm$xt$penalty
  levels <- rownames(S)
  x <- as.character(mf[[sm$term]])
  list(
    X = sparseMatrix(
      i = seq_along(x), j = match(x, levels), x = 1,
      dims = c(length(x), length(levels))
    ),
    S = list(Matrix(unname(S), sparse = TRUE)),
    N = attr(S, "null.space"),
    # the position of each level in the mgcv basis
    order = match(levels, levels(sm$knots))
  )
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
  ncons <- vapply(Gsetup$smooth, function(x) max(0, attr(x, "nCons")), 0)
  is_mrf <-
    vapply(Gsetup$smooth, function(x) {
      inherits(x, "mrf.smooth") && x$by == "NA" && !is.null(rownames(x$xt$penalty)) &&
        x$bs.dim == nrow(x$xt$penalty)
    }, TRUE)

  # the smooths without identifiability constraints, only needed for
  # constrained smooths which are not Markov random fields
  full <- NULL
  if (sparse && any(ncons[penalised] == 1 & !is_mrf[penalised])) {
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

    block <- NULL
    if (sparse && ncons[i] == 1 && length(intercept) == 1 &&
        inherits(attr(sm, "qrc"), "qr") && (is_mrf[i] || !is.null(full))) {
      if (is_mrf[i]) {
        parts <- mrf_parts(sm, Gsetup$mf)
      } else {
        parts <- list(X = full[[i]]$X, S = full[[i]]$S, order = seq_len(ncol(full[[i]]$X)))
      }
      con <- constraint_basis(sm, ncol(parts$X))
      con$Z <- con$Z[parts$order, , drop = FALSE]
      con$q1 <- con$q1[parts$order]
      if (ncol(con$Z) == length(cols) &&
          max(abs(parts$X %*% con$Z - G[, cols])) < 1e-8 * max(abs(G[, cols]))) {
        block <- smooth_block(parts$X, parts$S, centred = TRUE, N = parts$N)
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
      # use the mgcv basis and penalties, undoing the scaling of the
      # penalties. These are dense if the smooth is constrained.
      S <- Map(function(s, scale) s * scale, sm$S, sm$S.scale)
      block <- smooth_block(G[, cols, drop = FALSE], S)
      block$TG <- block$B
    }

    block$cols <- cols
    block$label <- sm$label

    if (ncol(block$Xf)) {
      colnames(block$Xf) <- paste0(sm$label, ".null", seq_len(ncol(block$Xf)))
    }
    if (is.null(block$Tint) && length(block$keep) == length(cols)) {
      colnames(block$Zu) <- Gsetup$term.names[cols]
    } else if (is_mrf[i] && !is.null(block$Tint)) {
      colnames(block$Zu) <- paste0(sm$label, ".", rownames(sm$xt$penalty)[block$keep])
    } else {
      colnames(block$Zu) <- paste0(sm$label, ".u", block$keep)
    }

    block
  })
}
