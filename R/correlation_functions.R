#' Get correlation matrix between two groups of variables.
#'
#' @param x A matrix, with variable names encoded in column names.
#' @param colInd1 An integer or character vector containing the indices or column names of the first group.
#' @param colInd2 An integer or character vector containing the indices or column names of the second group.
#' @return A correlation matrix as returned by `cor()`.

cor_mat <- function(x, colInd1, colInd2) {
  cor(x[, colInd1], x[, colInd2], use = 'pairwise.complete.obs')
}

#' Bootstrapped correlation
#'
#' Calculate the observed correlation and boostrap correlation quantiles
#' @inheritParams cor_mat
#' @param B Number of bootstrap replicates
#' @param alpha Significance level, e.g. 0.05, 0.1, etc (must be less than 0.5).
#' @param groupNames A character vector of length 2: the names of the two groups of variables, in that order.
#' @return A data.table with seven columns:
#' * The first two columns are the two groups, each row is a combination of variables.
#' * The next three columns are correlation results: low (alpha-th percentile), median, and high ((1-alpha)-th percentile).
#' * signif, which is TRUE if 0 is not between (low, high).
#' * rho0: the sampling correlation, which is the `t0` returned by `boot::tsboot()`.

cor_boot <- function(X, colInd1, colInd2, groupNames, B = 1000, alpha = 0.05) {
  stopifnot(alpha < 0.5)
  N <- nrow(X)
  # To handle the case where character vectors are supplied in colInd1 and colInd2
  names(varNames) <- varNames <- colnames(X)
  corBoot <- boot::tsboot(X, statistic = cor_mat,
                          colInd1 = colInd1, colInd2 = colInd2,
                          R = B, l = ceiling(sqrt(N)), sim = 'geom')
  rho <- matrixStats::colQuantiles(corBoot$t, probs = c(alpha, 0.5, 1 - alpha), type = 8, drop = FALSE)
  corDT <- cbind(CJ(varNames[colInd2], varNames[colInd1], sorted = FALSE), rho)
  setnames(corDT, c(rev(groupNames), 'low', 'median', 'high'))
  corDT[, ':='(rho0 = c(corBoot$t0),
               signif = !between(0, low, high))][]
}

grid_cor <- function(proxy.dt, clim.dt, proxi.var, clim.var, proxy.meta) {
  # Merge the two data sets
  DT <- merge(proxy.dt, clim.dt, by = 'year', allow.cartesian = TRUE)
  # Perform correlation for each (tree ring site, grid point) pair
  rho <- DT[, {
    ans <- cor.test(get(proxi.var), get(clim.var))
    list(r = ans$estimate, p = ans$p.value)
  },
  by = .(site, lon, lat)]
  # Get the coordinates of the tree ring sites for plotting
  rho <- merge(
    rho,
    proxy.meta[, .(site, tree_lon = long, tree_lat = lat)],
    by = 'site')
  rho[, signif := p < 0.05]
  # Determine the boundaries of the significant areas
  setkey(rho, lon, lat)
  rhoSignif <- rho[, signif_area(.SD, 0.5, 0.5), by = site]
  list(rho = rho, rhoSignif = rhoSignif)
}

plot_cor <- function(cor.list, cor.limits = NULL) {
  rho <- cor.list$rho
  rhoSignif <- cor.list$rhoSignif
  if (is.null(cor.limits)) cor.limits <- abs_range(rho$r)
  ggplot(rho) +
    # Correlation grid
    geom_tile(aes(lon, lat, fill = r)) +
    # Country borders
    geom_polygon(
      aes(long, lat, group = group),
      # This function calls the build-in world data set of ggplot2
      # and extracts the countries we need
      map_data('world', c('Vietnam', 'Laos', 'Cambodia', 'Thailand', 'Myanmar')),
      size = 0.1,
      colour = 'gray80',
      fill = NA) +
    # Significant area boundaries
    geom_segment(
      aes(lon - 0.25, lat + 0.25, xend = lon + 0.25, yend = lat + 0.25),
      rhoSignif[{top}],
      size = 0.1) +
    geom_segment(
      aes(lon - 0.25, lat - 0.25, xend = lon + 0.25, yend = lat - 0.25),
      rhoSignif[{bottom}],
      size = 0.1) +
    geom_segment(
      aes(lon - 0.25, lat - 0.25, xend = lon - 0.25, yend = lat + 0.25),
      rhoSignif[{left}],
      size = 0.1) +
    geom_segment(
      aes(lon + 0.25, lat - 0.25, xend = lon + 0.25, yend = lat + 0.25),
      rhoSignif[{right}],
      size = 0.1) +
    # Sites
    geom_point(aes(tree_lon, tree_lat), colour = 'black') +
    # Colour scale
    scale_fill_distiller(
      name = 'Correlation',
      palette = 'RdBu',
      direction = 1,
      limits = cor.limits,
      guide = guide_colorbar(barheight = 12)) +
    # Facets
    facet_wrap(
      vars(site),
      ncol = 5,
      labeller = as_labeller(siteLab)) +
    # Coordinate system
    coord_quickmap(expand = FALSE) +
    # Theme
    theme_map(font_size = 11) +
    panel_border('black', 0.2)
}

#' Get the boundaries of areas of significant correlation
#'
#' Returns columsn of top, down, bottom, left, TRUE for that respective border
signif_area <- function(cor.dt, dx, dy) {

  cor.dt[{signif}, list(lon, lat,
                        top    = .SD[.(lon,      lat + dy), is.na(r)],
                        bottom = .SD[.(lon,      lat - dy), is.na(r)],
                        left   = .SD[.(lon - dx, lat     ), is.na(r)],
                        right  = .SD[.(lon + dx, lat     ), is.na(r)])]
}

#' Get correlatin matrix between two groups of variables.
#'
#' @param x A matrix, with variable names encoded in column names.
#' @param colInd1 An integer or character vector containing the indices or column names of the first group.
#' @param colInd2 An integer or character vector containing the indices or column names of the second group.
#' @return A correlation matrix as returned by `cor()`.
