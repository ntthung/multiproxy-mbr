args <- commandArgs(trailingOnly = TRUE)
lambdas <- seq(0, 2, 0.2)
lambda <- lambdas[as.integer(args)]

pop <- 500
gen <- 500

library(magrittr)
library(glue)
library(data.table)
library(mbr)
library(GA)
library(mbr)

source('R/input_selection_functions.R')
p1tar <- readRDS('results/p1tar.RDS')
Xrwdo <- readRDS('results/Xrwdo.RDS')
poolDT <- readRDS('results/poolDT.RDS')

set.seed(24)
cvFolds <- make_Z(1922:2003, nRuns = 50, frac = 0.25, contiguous = TRUE)

lambdaChar <- if (lambda < 1) paste0('0', lambda * 10) else as.character(lambda * 10)

lapply(seq(100, 400), function(s) {

  siteOptim <- ga(
    type = 'binary',
    fitness = memoise::memoise(cv_site_selection),
    pool = poolDT,
    Xpool = Xrwdo,
    instQ = p1tar,
    cv.folds = cvFolds,
    start.year = 1750,
    lambda = lambda,
    log.trans = 1:3,
    popSize = pop,
    maxiter = gen,
    run = min(c(gen, 100)),
    parallel = TRUE,
    monitor = FALSE,
    nBits = nrow(poolDT),
    seed = s)

  saveRDS(siteOptim, glue('results/pop{pop}_gen{gen}_lambda{lambdaChar}_s{s}.RDS'))

})
