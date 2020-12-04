lambda <- 0
pop <- 500
gen <- 500

source('R/init.R')
source('R/correlation_functions.R')
library(GA)
library(mbr)

if (!require(memoise)) install.packages('memoise')

if (!require(mbr)) {
  if (!require(remotes)) install.packages('remotes')
  remotes::install_github('ntthung/mbr')
}

# Streamflow ------------------------------
p1merge <- fread('results/p1merge.csv', key = 'season')

## Hinkley D and transformation type
p1HD <- p1merge[, .(D   = round(hinkley(Qa), 2),
                    Dlog = round(hinkley(log(Qa)), 2)),
                keyby = season
              ][, trans := fifelse(abs(D) > abs(Dlog), 'log', 'none')]

# Proxies -------------------------------

oxi <- fread('data/oxi.csv', key = 'site')
crn <- fread('data/crn.csv', key = 'site')

crnMeta <- fread('data/crnMeta.csv', key = 'site')
oxiMeta <- fread('data/oxiMeta.csv', key = 'site')

# In-filling
crnWideFull   <- crn[year %in% 1748:2005, dcast(.SD, year ~ site, value.var = 'rwi')]
crnMatFull    <- as.matrix(crnWideFull[, -'year'])
crnMatFilled  <- imputePCA(crnMatFull, ncp = 19)$completeObs
crnWideFilled <- as.data.table(crnMatFilled)

oxiWideFull   <- oxi[year %in% 1748:2005, dcast(.SD, year ~ site, value.var = 'do18')]
oxiMatFull    <- as.matrix(oxiWideFull[, -'year'])
oxiMatFilled  <- imputePCA(oxiMatFull, ncp = 3)$completeObs
oxiWideFilled <- as.data.table(oxiMatFilled)

# Input matrices
allLags <- -2:2 # Core years: 1750:2003, and shift back or forth
Xrw <- do.call(cbind, lapply(allLags, function(l) {
  x <- crnMatFilled[3:256 - l, ]
  colnames(x) <- paste0(colnames(x), l)
  x
}))
Xdo <- do.call(cbind, lapply(allLags, function(l) {
  x <- oxiMatFilled[3:256 - l, ]
  colnames(x) <- paste0(colnames(x), l)
  x
}))
Xrwdo <- cbind(Xrw, Xdo)

# Cross-validation -------------------------
set.seed(24)
cvFolds <- make_Z(1922:2003, nRuns = 50, frac = 0.25, contiguous = TRUE)

# Labelling ---------------------------------
oxiLab <- function(s) oxiMeta[s, display]
crnLab <- function(s) crnMeta[s, display]

lambdaLab <- function(x) sapply(x, function(v) paste0('\u03bb = ', substr(v, 3, 5)))
ssnLabBrief <- function(x) seasonClass[x, type]
names(seasons) <- seasons <- c('NJ', 'JO', 'WY')

seasonClass <- fread('data/season_class.csv', key = 'season')

# Streamflow: 1922 to 2003
# Targets: 1922 to 2003
QDry <- p1merge['NJ'][1:82]
QWet <- p1merge['JO'][1:82]
QAnn <- p1merge['WY'][1:82]
p1tar <- rbind(QDry, QWet, QAnn)
p1tar[, season := factor(season, seasons)]
setkey(p1tar, season)

QMat    <- log(as.matrix(dcast(p1tar, year ~ season, value.var = 'Qa')[, -'year']))
cor_lag <- function(X, Y, l, Xyears, Yyears, alpha) {
  ind <- which(Xyears %in% Yyears) - l
  cor_boot(cbind(X[ind, ], Y),
           1:ncol(X),
           1:ncol(Y) + ncol(X),
           c('site', 'season'),
           alpha = alpha)
}
set.seed(42)
rhoCrn <- lapplyrbind(allLags, function(l) {
  DT <- cor_lag(crnMatFilled, QMat, l, 1748:2005, 1922:2003, 0.05)
  DT[, lag := l]
  DT[, alpha := 0.05][]
})
rhoOxi <- lapplyrbind(allLags, function(l) {
  DT <- cor_lag(oxiMatFilled, QMat, l, 1748:2005, 1922:2003, 0.05)
  DT[, lag := l]
  DT[, alpha := 0.05][]
})

dt <- CJ(s = factor(seasons, seasons), l = allLags)
dt[, sl := paste0(s, '(', l, ')')]

rhoCrn[, season_lag := paste0(season, '(', lag, ')')]
rhoCrn[, season_lag := factor(season_lag, levels = dt$sl)]
rhoCrn[, season := factor(season, levels = seasons)]
rhoCrn[, lag := factor(lag, levels = allLags)]
rhoOxi[, season_lag := paste0(season, '(', lag, ')')]
rhoOxi[, season_lag := factor(season_lag, levels = dt$sl)]
rhoOxi[, season := factor(season, levels = seasons)]
rhoOxi[, lag := factor(lag, levels = allLags)]

poolDT <- rbind(rhoCrn[{signif} & abs(rho0) >= 0.2, .(site, lag), by = season],
                rhoOxi[{signif} & abs(rho0) >= 0.2, .(site, lag), by = season]
  )[, site := paste0(site, lag)]
setkey(poolDT, season)

Xused <- Xrwdo[, poolDT[, unique(site)]]

lambdaChar <- fcase(
  lambda == 0, '00',
  lambda == 1, '10')

siteOptim <- ga(
  type = 'binary',
  fitness = memoise::memoise(mbr::cv_site_selection),
  pool = poolDT,
  seasons = seasons,
  Xpool = Xused,
  instQ = p1tar,
  cv.folds = cvFolds,
  start.year = 1750,
  lambda = lambda,
  log.trans = 1:3,
  num.targets = 3,
  popSize = pop,
  maxiter = gen,
  run = min(c(gen, 100)),
  parallel = TRUE,
  monitor = FALSE,
  nBits = nrow(poolDT))

saveRDS(siteOptim, glue('results/pop{pop}_gen{gen}_lambda{lambdaChar}_refactored.RDS'))
