# ---- Load packages and source files-------------

library(magrittr)
library(glue)
library(data.table)

library(qmap)
library(ldsr)
library(missMDA)
library(mbr)
library(GA)

library(cowplot)
library(ggplot2)
library(ggtext)
library(patchwork)

# Plot utilities ---------------------------------

## Coordinate labels
pasteLong <- function(x) {
  lest::case_when(x < 0   ~ paste0(-x, "\u00b0W"),
                  x > 180 ~ paste0(360 - x, "\u00b0W"),
                  x %in% c(0, 180)  ~ paste0(x, '\u00b0'),
                  TRUE ~ paste0(x, "\u00b0E"))
}

pasteLat <- function(x) {
  ifelse(x < 0, paste0(-x, "\u00b0S"), ifelse(x > 0, paste0(x, "\u00b0N"), paste0(x, '\u00b0')))
}

skip_label <- function(n) {
  function(x) {
    idx <- seq_along(x)
    x[idx %% n != 1] <- ' '
    x
  }
}

my_theme <- theme_cowplot(font_size = 11) +
  theme(
    line = element_line(size = 0.2),
    plot.tag = element_text(face = 'bold'),
    strip.background = element_rect('gray95'),
    strip.text = element_text(face = 'bold'),
    axis.line = element_line(size = 0.2),
    axis.ticks = element_line(size = 0.2))
theme_set(my_theme)

monthLab <- month.abb
names(monthLab) <- 1:12

lambdaLab <- function(x) sapply(x, function(v) paste0('\u03bb = ', substr(v, 1, 1)))

dark2 <- RColorBrewer::brewer.pal(3, 'Dark2')[c(2, 1, 3)]

pal <- wesanderson::wes_palette('Cavalcanti1')[c(2, 1)]

# Handling time -----------------------------------------------------------

doy_to_month <- function(d, type = 'char') {
  months <- month(as.Date(d - 1, origin = '1987-01-01'))
  if (type == 'numeric') months else month.abb[months]
}

doy_to_ddMMM <- function(d) format(as.Date(d - 1, origin = '1987-01-01'), '%b %d')

doy_month_label <- function(x) sapply(x, function(d) paste0(d, ' (', doy_to_ddMMM(d),')'))

get_seasons <- function(x, var, out.name, months, months.fwd = NULL, months.bwd = NULL, func = sum) {
  dt <- copy(x[month %in% months])
  dt[month %in% months.fwd, year := year + 1]
  dt[month %in% months.bwd, year := year - 1]
  ans <- dt[, .(out = if (.N == length(months)) func(get(var)) else numeric(0)), by = year]
  setnames(ans, 'out', out.name)
  ans[]
}

# Other utilities ---------------------------------

distance <- function(a, b) sqrt(sum((a - b)^2))

normalize <- function(x, ...) {(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

rm_null <- function(x) x[!sapply(x, is.null)]

abs_range <- function(x) {
  absMax <- max(abs(range(x)))
  c(-absMax, absMax)
}

hinkley <- function(x, type = 2) {
  scale <- if (type == 1) sd(x, na.rm = TRUE) else diff(quantile(x, c(0.25, 0.75), na.rm = TRUE, type = 8, names = FALSE))
  (mean(x, na.rm = TRUE) - median(x, na.rm = TRUE)) / scale
}

roundDT <- function(DT, digits = 2, type = 'numeric') {
  DT[, lapply(.SD,
              function(x) {
                if (is.numeric(x)) {
                  rd <- round(x, digits)
                  if (type == 'numeric') rd else sprintf(glue("%.{digits}f"), rd)
                } else x
              })]
}

power_set <- function(N, min.size = 1, max.size = N) {
  if (max.size > N) max.size <- N
  if (min.size > N) NULL else {
    subsets <- foreach(k = 1:N) %dopar% combn(1:N, k, simplify = FALSE)
    subsets <- unlist(subsets, recursive = FALSE)
    lengths <- sapply(subsets, length)
    subsets[lengths %between% c(min.size, max.size)]
  }
}

'%ni%' <- Negate('%in%')

lapplyrbind <- function(x, fun, ..., id = NULL) rbindlist(lapply(x, fun, ...), idcol = id)

# Megadroughts

mgd <- data.table(sta = c(1756, 1790, 1876),
                  fin = c(1768, 1796, 1878),
                  name = c('Strange Parallels Drought',
                           'East India Drought',
                           'Victorian Great Drought'))

standardize <- function(x, ...) (x - mean(x, ...)) / sd(x, ...)

boxcox_lambda <- function(x) {
  bc <- MASS::boxcox(x ~ 1, plotit = FALSE)
  bc$x[which.max(bc$y)]
}

boxcox <- function(x, lambda = NULL) {
  if (is.null(lambda)) lambda <- boxcox_lambda(x)
  if (abs(lambda) < 0.01) log(x) else (x^lambda - 1) / lambda
}
