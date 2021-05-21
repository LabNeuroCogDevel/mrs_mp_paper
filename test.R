#!/usr/bin/env Rscript
library(testthat)
source('funcs.R')

test_that("remove mjaor met outliers, keep good",{
  # data will only two good rows
  # each row has only 1 exclusion criteral met
  xd<- data.frame(
    id          = c(1,  2, 3, 4, 5, 6),
    # >= 10, will be removed
    GPC.Cho.SD  = c(11, 1, 1, 1,NA, 1),
    NAA.NAAG.SD = c( 1,11, 1, 1, 1, 1),
    Cr.SD       = c( 1, 1,11, 1, 1, 1),
    # >= 3  will be removed
    MM20.Cr     = c( 1, 1, 1, 4, 1, 1))

  # only the last to rows will remain:
  expect_equal(remove_major_met_outliers(xd)$id,  c(5,6))
})

test_that("revese 216 row/col", { 
    expect_equal(reverse_index(216),1)
    expect_equal(reverse_index(108),109)
    expect_equal(reverse_index(1)  ,216)
})

test_that("read in data",{

})

test_that("mk_age_resid regression",{
    xd <- data.frame(ex=1:3, invage=1:3) # prefect fit
    xd_reg <-  mk_age_resid(xd, 'ex') 
    expect_equal(xd_reg$ex, c(0, 0, 0))
})
