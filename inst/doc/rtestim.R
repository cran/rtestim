## ----setup, collapse=TRUE-----------------------------------------------------
library(rtestim)
library(ggplot2)
theme_set(theme_bw())

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  # dpi = 300,
  collapse = FALSE,
  comment = "#>",
  fig.asp = 0.618,
  fig.width = 6,
  fig.align = "center",
  out.width = "80%"
)

## ----fig.align='center'-------------------------------------------------------
set.seed(12345)
case_counts <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
ggplot(data.frame(x = 1:101, case_counts), aes(x, case_counts)) +
  geom_point(colour = "cornflowerblue") +
  labs(x = "Time", y = "Case Counts")

## ----fig.align='center'-------------------------------------------------------
mod <- estimate_rt(observed_counts = case_counts, nsol = 20)
plot(mod)

## -----------------------------------------------------------------------------
mod_cv <- cv_estimate_rt(observed_counts = case_counts)

## ----fig.align='center'-------------------------------------------------------
plot(mod_cv)

## ----fig.align='center'-------------------------------------------------------
plot(mod_cv, which_lambda = "lambda.1se")

## -----------------------------------------------------------------------------
observation_incr <- rpois(101, lambda = 2)
observation_incr[observation_incr == 0] <- 1
observation_time <- cumsum(observation_incr)
ggplot(data.frame(x = observation_time, case_counts), aes(x, case_counts)) +
  geom_point(colour = "cornflowerblue") +
  labs(x = "Time", y = "Case Counts")

## ----fig.align='center'-------------------------------------------------------
mod <- estimate_rt(observed_counts = case_counts, x = observation_time)
plot(mod)

## ----fig.align='center'-------------------------------------------------------
mod <- estimate_rt(observed_counts = case_counts, korder = 0, nsol = 20)
plot(mod)

## ----warning=FALSE, fig.align='center'----------------------------------------
can <- estimate_rt(
  observed_counts = cancovid$incident_cases,
  x = cancovid$date,
  korder = 2,
  nsol = 20,
  maxiter = 1e5
)

plot(can) + coord_cartesian(ylim = c(0.5, 2))

## ----fig.align='center'-------------------------------------------------------
can_cb <- confband(can, lambda = can$lambda[10], level = c(.5, .8, .95))
can_cb
plot(can_cb) + coord_cartesian(ylim = c(0.5, 2))

