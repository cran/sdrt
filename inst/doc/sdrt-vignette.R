## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  install.packages("sdrt")
#  library(sdrt)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  library(devtools)
#  install_github("TharinduPDeAlwis/sdrt")
#  library(sdrt)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  data("lynx")
#  y <- log10(lynx)
#  p_list=seq(2,5,by=1)
#  fit.model=pd.boots(y,p_list,w1=0.1,B=10)
#  fit.model$dis_pd
#  fit.model$p_hat
#  fit.model$d_hat

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  set.seed(1)
#  data("lynx")
#  y <- log10(lynx)
#  p <- 3
#  d <- 1
#  w1_list=seq(0.1,0.5,by=0.1)
#  Tunning.model=sigma_u(y, p, d, w1_list=w1_list, std=FALSE, B=10)
#  Tunning.model$sigma_u_hat

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  library(sdrt)
#  data("lynx")
#  y <- log10(lynx)
#  p <- 3
#  d <- 1
#  fit.model <- sdrt(y, p, d=1,method="FM",density = "kernel")
#  fit.model$eta_hat

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  library(pracma)
#  data("lynx")
#  y <- log10(lynx)
#  p <- 4
#  d <- 1
#  fit.model <- sdrt(y, p, d,method="NW")
#  fit.model$eta_hat

