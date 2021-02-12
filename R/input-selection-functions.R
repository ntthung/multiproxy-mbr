#' Weighted PCA based on correlation, a la Cook et al (2010). Each column of X will be weighted by an exponent r^p, i.e.
#' \eqn{W_i = X_i*r_i^p}
#' @param X A matrix of variables
#' @param rho Vector of correlations between each variable and the output, ignored if `rho == 1`
#' @param p Weight exponent
#' @param min.variance Minimum fraction of variance explained by the retained PCs
#' @param use.eigen IF TRUE, only PCs whose eigenvalue > 1 are kept.
#' @param return.matrix If `TRUE`, a matrix is returned, otherwise a `data.table` is returned.
#' @export
wPCA <- function(X, rho = 1, p = 0, min.variance = 0.95, use.eigen = TRUE, return.matrix = FALSE) {

  # Step 1: weighted PCA
  # Negative correlation will cause a problem when raised to a real power
  if (length(rho) > 1) X <- sapply(1:dim(X)[2], function(k) X[, k]*abs(rho[k])^p)
  # According to Cook et al (2010) SI, when p = 0 we do PCA with correlation matrix form,
  # and when p != 0 we do PCA with covariance matrix form.
  pca.model <- summary(stats::prcomp(X, scale. = (p == 0)))

  # Step 2: retain PCs
  n <- which(pca.model$importance['Cumulative Proportion', ] > min.variance)[1]
  if (use.eigen) {
    n2 <- max(which(pca.model$sdev >= 1))
    if (n2 < n) n <- n2
  }
  # If there is only one column, '[' returns a vector instead of a matrix if drop = TRUE
  out <- pca.model$x[ , 1:n, drop = FALSE]
  if (return.matrix) out else data.table(out)
}

#' Input variable selection
#'
#' Four methods are in use: Variable Selection Using Random Forests (VSURF) from the `VSURF` package, forward stepwise selection from the `FWDselect` package, and forward and backward selection using branch and bound from the `leaps` package. The latter three use BIC as the selection criterion.
#' @param X Input variables, can be a matrix or a data.frame. One variable per column. Must have column names.
#' @param Y Output variable
#' @param method Either "VSURF", "FWDselect", leaps forward", or "leaps backward".
#' @param nvmax Maximum number of variables; not applicable for the VSURF method.
#' @param parallel Only applicable if VSURF or FWDselect is used. Default = `FALSE`.
#' @return If X has column names and `use.name == TRUE`, the names of the selected variables are returned, otherwise the selected column numbers are returned.
#' @export
 input_selection <- function(X, Y, method, nvmax, parallel = FALSE) {

  colNames <- colnames(X)
  if (is.null(colNames)) {
    stop("Error in input_selection: X must have column names")
  }
  colNums <- if (ncol(X) == 1) 1 else {
    switch(method,
           'VSURF' = {
             vsurf.fit <- VSURF::VSURF(X, Y, na.action = stats::na.omit, parallel = parallel, verbose = FALSE)
             sv <- if (is.null(vsurf.fit$varselect.pred)) vsurf.fit$varselect.interp else vsurf.fit$varselect.pred
             sort(sv)
           },
           'FWDselect' = {
             nv <- min(ncol(X) - 1, nvmax) # ncol(X) > 1 here
             #FWDselect doesn't work with data.table
             X <- as.data.frame(X)
             fwd.fit <- FWDselect::qselection(X, Y, 1:nv, criterion = 'bic', cluster = parallel, )
             sv <- fwd.fit$selection[which.min(fwd.fit$bic)]
             # qselection returns the selected variables in a single string,
             #$ where the variable names are separated by commas, so we need to do strsplit
             sv <- strsplit(as.character(sv), ', ')[[1]]
             sort(match(sv, colNames))
           },
           'leaps forward' = {
             ivs <- summary(leaps::regsubsets(X, Y, method = 'forward', nvmax = nvmax))
             idx <- which.min(ivs$bic)
             which(ivs$which[idx, -1])
           },
           'leaps backward' = {
             ivs <- summary(leaps::regsubsets(X, Y, method = 'backward', nvmax = nvmax))
             idx <- which.min(ivs$bic)
             which(ivs$which[idx, -1])
           },
           stop("Error in input_selection: method not supported; method can only be 'VSURF', 'FWDselect', 'leaps forward', or 'leaps backward'"))
  }
  colNames[colNums]
}

#' Calculate the objective function value in cross-validation
#'
#' This function returns the negative of the objective function value, so the outer optimization should minimize it. his is consistent with the behaviour of the GA package. Please take note of this when you use another optimizer.
#'
#' @param choice Binary vector of choices
#' @param pool Data.table listing the sites with significant correlations for each season
#' @param n.site.min Minimum number of sites to be kept for each season
#' @inheritParams cv_mb
#' @param Xpool The input matrix of all the sites in the pool
#' @param use.robust.mean Whether the arithmetic mean (if set to FALSE) or the Tukey's biweight robust mean is used to etermine the final value from cross-validation. Default is FALSE.
#' @export
cv_site_selection <- function(choice, pool, n.site.min = 5,
                             Xpool, instQ, cv.folds, start.year,
                             lambda = 1, log.trans, force.standardize = FALSE,
                             use.robust.mean = FALSE) {

 season <- N <- site <- NULL
 seasons <- levels(instQ$season)
 poolSubset <- pool[choice == 1L]
 num.targets <- length(seasons)

 # Constraint: at least n.site.min for each season
 # If not met, return -1e12
 n.sites <- poolSubset[, .N, by = season][, N]

 if (length(n.sites) == num.targets && all(n.sites > n.site.min)) {

   years     <- start.year:max(instQ$year)
   instInd   <- which(years %in% instQ$year)

   pcList <- lapply(seasons, function(s) {
     Qa <- instQ[s]
     X <- Xpool[, poolSubset[s, site]]
     PC <- wPCA(X, use.eigen = FALSE, return.matrix = TRUE)
     sv <- input_selection(
       PC[instInd, , drop = FALSE],
       Qa$Qa, 'leaps backward', nvmax = 8)
     PC[, sv, drop = FALSE]
   })

   cvFval <- cv_mb(instQ, pcList, cv.folds, start.year, lambda, log.trans, force.standardize, return.type = 'mb')

   if (use.robust.mean) -tbrm(cvFval) else -mean(cvFval)
 } else -1e12
}
