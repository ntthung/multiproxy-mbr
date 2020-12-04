# ---- Load packages and source files-------------

library(magrittr)
library(ggplot2)
library(patchwork)
library(ggtext)
library(data.table)
library(ldsr)
library(foreach)
library(missMDA)
library(doParallel)
library(glue)
library(mbr)
library(GA)
library(cowplot)
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

skip_even <- function(x) {
  idx <- seq_along(x)
  x[idx %% 2 == 0] <- ' '
  x
}

skip_5 <- function(x) {
  idx <- seq_along(x)
  x[idx %% 5 != 1] <- ' '
  x
}

skip_label <- function(n) {
  function(x) {
    idx <- seq_along(x)
    x[idx %% n != 1] <- ' '
    x 
  }
}

my_theme <- theme_classic() +
  theme(
    line = element_line(size = 0.2),
    plot.tag = element_text(size = 10, face = 'bold'),
    strip.background = element_blank(),
    strip.text = element_text(face = 'bold'),
    axis.title = element_text(size = 10),
    axis.line = element_line(size = 0.2),
    axis.ticks = element_line(size = 0.2))
theme_set(my_theme)

# Read met data -----------------------------------------------------------

#' Convert MATLAB data format to tidy data format
#' MATLAB data for each variable come in as a matrix nyears x 366 (one column for each day)
#' and a separate vector for the years
#' @param years As provided by MATALB output
#' @param mat As provided by MATLAB output
#' @param var.name The name of the variable (Tmax, Tmin, T, P)
#' @return A data.table (year, month, doy, variable, value)
year_mat_to_dt <- function(years, mat, var.name) {
  dt <- data.table(
    year = rep(years, 366),
    doy = rep(1:366, each = length(years)),
    value = c(mat)
  )
  dt[, month := doy_to_month(doy)]
  dt[, variable := var.name]
  dt[]
}

#' Read MATLAB data for one met station
#' 
#' @param code The stations' 3-letter code as provided by Brendan
#' @return A data.table (station, variable, year, month, doy, value)
read_station <- function(code) {
  file <- paste0('data/met/', code, '.mat')
  dat  <- R.matlab::readMat(file)
  Tmax <- year_mat_to_dt(dat$yr.dmax, dat$dmax, 'Tmax')
  Tmin <- year_mat_to_dt(dat$yr.dmin, dat$dmin, 'Tmin')
  TT   <- year_mat_to_dt(dat$yr.dtemp, dat$dtemp, 'T')
  P    <- year_mat_to_dt(dat$yr.drain, dat$drain, 'P')
  dt <- rbind(Tmax, Tmin, TT, P)
  dt[, station := code]
  dt[]
}


# Read Thai streamflow data -----------------------------------------------

read_TH_Q <- function(stationName) {
  dt <- fread(paste0('data/streamflow/', stationName, '-monthly-Q.csv'))
  dt <- melt(dt, id.var = 'year', variable.name = 'month.name', value.name = 'Q', variable.factor = FALSE)
  dt[, station := paste(substr(stationName, 1, 1), substr(stationName, 2, 3), sep = '.')
   ][, month := match(month.name, month.abb)
   ][month %in% 1:3, year := year + 1
   ][, month.name := factor(month.name, levels = month.abb)]
  setcolorder(dt, 'station')
  dt[order(year, month)]
}

# PCA ---------------------------------------------------------------------

plotPCA <- function(pcaModel) {
  pcaRes <- pcaModel$importance %>% 
    t %>% 
    as.data.table(keep.rownames = TRUE) %>% 
    setnames(c('PC', 'eigen', 'varProp', 'cumuProp'))
  pcaRes[, PC := factor(PC, levels = paste0('PC', 1:.N))]
  p1 <- ggplot(pcaRes, aes(PC, eigen)) +
    geom_bar(stat = 'identity', fill = 'gray') +
    geom_hline(yintercept = 1) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = 'eigenvalue')
  p2 <- ggplot(pcaRes, aes(PC, varProp*100)) +
    geom_bar(stat = 'identity', fill = 'gray') +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = 'Variance explained [%]')
  p3 <- ggplot(pcaRes, aes(PC, cumuProp*100)) +
    geom_bar(stat = 'identity', fill = 'gray') +
    geom_hline(yintercept = 95) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = 'Cumulative variance [%]')
  p1 + p2 + p3 & theme_classic()  
}

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
