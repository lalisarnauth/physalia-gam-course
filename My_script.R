
pkgs <- c("gamm4", "tidyverse", "readxl", "mgcViz", "DHARMa", "gratia",
          "marginaleffects", "ggforce")

install.packages(pkgs)

# Day 1 ----

library(readr)
library(here)
library(ggplot2)
library(tidyverse)

#### Darlington ----

wasp <- read_csv(here("data", "darlingtonia.csv"), comment = "#",
                 col_types = "dl")
m <- glm(visited ~ leafHeight, data = wasp, family = binomial)
m
summary(m)

# data to predict at
pdat <- with(wasp,
             tibble(leafHeight = seq(min(leafHeight),
                                     max(leafHeight),
                                     length = 100)))
# predict
pred <- predict(m, pdat, type = "link", se.fit = TRUE)
ilink <- family(m)$linkinv # g-1()
pdat <- pdat %>%
  bind_cols(data.frame(pred)) %>%
  mutate(fitted = ilink(fit),
         upper = ilink(fit + (2 * se.fit)),
         lower = ilink(fit - (2 * se.fit)))
# plot
ggplot(wasp, aes(x = leafHeight,
                 y = as.numeric(visited))) +
  geom_point() +
  geom_ribbon(aes(ymin = lower, ymax = upper,
                  x = leafHeight), data = pdat,
              inherit.aes = FALSE, alpha = 0.2) +
  geom_line(data = pdat, aes(y = fitted)) +
  labs(x = "Leaf Height [cm]",
       y = "Probability of visitation")

#### Peat Bog ----

maddy <- read_csv(here("data", "maddy-peat.csv"), col_types = "cdddddd")
maddy <- mutate(maddy, midDepth = upperDepth - (0.5 * abs(upperDepth - lowerDepth)),
                calMid = calUpper - (0.5 * abs(calUpper - calLower)))
maddy

ggplot(maddy, aes(x = midDepth, y = calMid)) +
  geom_point() +
  labs(y = "Calibrated Age", x = "Depth")

m_gamma <- glm(calMid ~ midDepth, data = maddy, family = Gamma(link = "identity"))
summary(m_gamma)

# data to predict at
pdat <- with(maddy,
             tibble(midDepth = seq(min(midDepth),
                                   max(midDepth),
                                   length = 100)))
# predict
p_gamma <- predict(m_gamma, pdat, type = "link",
                   se.fit = TRUE)
ilink <- family(m_gamma)$linkinv
# confidence interval
p_gamma <- pdat %>%
  bind_cols(data.frame(p_gamma)) %>%
  mutate(fitted = ilink(fit),
         upper = ilink(fit + (2 * se.fit)),
         lower = ilink(fit - (2 * se.fit)))
# plot
p1 <- ggplot(maddy, aes(x = midDepth, y = calMid)) +
  geom_ribbon(aes(ymin = lower, ymax = upper,
                  x = midDepth), data = p_gamma,
              inherit.aes = FALSE, alpha = 0.2) +
  geom_line(data = p_gamma, aes(y = fitted)) +
  geom_point() +
  labs(y = "Calibrated Age", x = "Depth",
       title = "Gamma GLM")
p1

# fit gaussian GLM
m_gaus <- glm(calMid ~ midDepth, data = maddy,
              family = gaussian)
# predict
p_gaus <- predict(m_gaus, pdat, type = "link",
                  se.fit = TRUE)
ilink <- family(m_gaus)$linkinv
# prep confidence interval
p_gaus <- pdat %>%
  bind_cols(data.frame(p_gaus)) %>%
  mutate(fitted = ilink(fit),
         upper = ilink(fit + (2 * se.fit)),
         lower = ilink(fit - (2 * se.fit)))
# plot
p2 <- ggplot(maddy, aes(x = midDepth, y = calMid)) +
  geom_ribbon(aes(ymin = lower, ymax = upper,
                  x = midDepth), data = p_gaus,
              inherit.aes = FALSE, alpha = 0.2) +
  geom_line(data = p_gaus, aes(y = fitted)) +
  geom_point() +
  labs(y = "Calibrated Age",
       x = "Depth",
       title = "Linear Model")
p2

library("patchwork")
p1 + p2

#### Diamondback moths ----

# DBM are pests
# wasps control them
# Treatment greenhouses have chemical products (HIPVs) that attracts wasps
# They used traps to count the number of wasps and DBMs (unit of analysis)
# We will model numbers of C. vestalis against numbers of DBM and treatment using each trap as the units of analysis.
# Response variable is the number of DBM

moth <- readr::read_csv("data/uefunex.csv")

library("dplyr")
moth |>
  group_by(treatment) |>
  summarise(n = n(), mean = mean(parasitoid), median = median(parasitoid),
            sd = sd(parasitoid))

moth |>
  ggplot(aes(x = treatment, y = parasitoid)) +
  geom_violin(aes(fill = treatment))

# Scaled
moth |>
  ggplot(aes(x = treatment, y = parasitoid)) +
  geom_violin(aes(fill = treatment)) +
  scale_y_sqrt()

moth |>
  ggplot(aes(x = treatment, y = parasitoid)) +
  geom_violin(aes(fill = treatment)) +
  scale_y_continuous(trans = ggforce::power_trans((1/4)))
# What distribution would you expect the response variable parasitoid to follow?

# Poisson

moth_glm1 <- glm(parasitoid ~ moth + treatment + moth:treatment,
                 family = poisson, data = moth)

summary(moth_glm1)

# analysis of deviance table

anova(moth_glm1, test = "LRT")

# plot

moth |>
  ggplot(aes(y = parasitoid, x = moth, color = treatment)) +
  geom_jitter(stat = "identity", width = 0.05, height = 0.05) +
  geom_smooth(method = "glm", method.args = list(family = "poisson")) +
  theme(legend.position = "bottom")

# Problems
moth |>
  group_by(treatment) |>
  summarise(n = n(), mean = mean(parasitoid), median = median(parasitoid),
            sd = sd(parasitoid))

# Overdispersion
presid <- resid(moth_glm1, type = "pearson")
n <- nrow(moth)
params <- length(coef(moth_glm1))
disp <- sum(presid^2) / (n - params)
disp

moth_glm2 <- glm(parasitoid ~ moth + treatment + moth:treatment,
                 family = quasipoisson, data = moth)
summary(moth_glm2)

# negative binomial

library("mgcv")
moth_glm3 <- gam(parasitoid ~ moth + treatment + moth:treatment,
                 family = nb(), method = "ML", data = moth)
summary(moth_glm3)

# negbin 2

library("glmmTMB")
moth_glm4 <- glmmTMB(parasitoid ~ moth + treatment + moth:treatment,
                     family = nbinom2("log"), REML = FALSE, data = moth)
summary(moth_glm4)

library("marginaleffects")
moth_glm4 |>
  plot_predictions(condition = c("moth", "treatment"),
                   vcov = TRUE, type = "response")

#### GAMs ----

library("readr")
library("ggplot2")
library("dplyr")
library("mgcv")
library("gratia")
library("marginaleffects")

URL <-  "https://bit.ly/hadcrutv4"
gtemp <- read_table(URL, col_types = 'nnnnnnnnnnnn', col_names = FALSE) %>%
  select(num_range('X', 1:2)) %>% setNames(nm = c('Year', 'Temperature'))

m <- gam(Temperature ~ s(Year), data = gtemp, method = 'REML')
summary(m)

# plot
plt_labs <- labs(x = "Year", y = expression(Temperature ~ degree * C))
gtemp |>
  ggplot(
    aes(x = Year, y = Temperature)
  ) +
  geom_line() +
  geom_point() +
  plt_labs

# fit the Gaussian GAM
m_gtemp <- gam(
  Temperature ~ s(Year),
  data = gtemp,
  method = "REML",
  family = gaussian()
)

# model summary
summary(m_gtemp)
# or overview
m_gtemp |> overview()

# plot the estimate smooth
draw(m_gtemp, residuals = TRUE, rug = FALSE)

# plot the fitted smooth on the original data
#
# generate new data to predict at
newd <- m_gtemp |>
  data_slice(
    Year = evenly(Year, n = 200)
  )

# use fitted_values to get predictions on the response scale
fv_gtemp <- m_gtemp |>
  fitted_values(
    data = newd, scale = "response"
  )

# plot
fv_gtemp |>
  ggplot(aes(x = Year, y = .fitted)) +
  geom_point(data = gtemp, aes(x = Year, y = Temperature)) +
  geom_ribbon(
    aes(
      ymin = .lower_ci, ymax = .upper_ci, x = Year
    ),
    alpha = 0.4,
    inherit.aes = FALSE,
    fill = "#fdb338"
  ) +
  geom_line(linewidth = 1, colour = "#025196") +
  plt_labs

# or with conditional_values
m_gtemp |>
  conditional_values(
    condition = "Year"
  ) |>
  draw() +
  geom_point(
    data = gtemp,
    aes(x = Year, y = Temperature)
  ) +
  plt_labs

# or with plot_predictions
m_gtemp |>
  plot_predictions(
    condition = "Year",
    points = 1
  ) +
  plt_labs

# Day 2 ----

# simply moving to higher-order polynomials doesn't always improve accuracy

model <- gam(
  y ~ s(x1) + s(x2) + te(x3, x4), # formula describing model
  data = my_data_frame,           # your data
  method = "REML",                # or "ML"
  family = gaussian               # or something more exotic
)

# s() terms are smooths of one or more variables
# te() terms are the smooth equivalent of main effects + interactions

library("mgcv")
gam(Temperature ~ s(Year, k = 10), data = gtemp, method = "REML")
# Estimated degrees of freedom: 7.84 - not linear

#### Splines ####
# functions composed of simpler functions
# Simpler functions are basis functions & the set of basis functions is a basis
# The sum of simpler functions = complex model
# K is the max # of basis functions

#### Maximise penalised log-likelihood → β ####
# We estimate the coefficients β by finding the best fit to the data while penalizing overly complex (wiggly) curves
# Fit - penalty
# Wiggly curves are more penalized
# Once smoothing is applied, curves have fewer effective degrees of freedom (EDF)

# high λ = linear (10000 λ)
# low λ = wiggliness (1e-05)

## Portuguese crestlarks ----

library(here)

larks <-  read_csv(here("data", "larks.csv"),
                   col_types = "ccdddd")
larks <- larks |>
  mutate(crestlark = factor(crestlark),
         linnet = factor(linnet),
         e = x / 1000,
         n = y / 1000)
head(larks)
# e and n are spacial coordinates

crest <- gam(crestlark ~ s(e, n, k = 100),
             data = larks,
             family = binomial,
             method = "REML")

residuals_linpred_plot(crest, type = "quantile")

# He's asking: Does the probability of finding Crestellark change spatially depending on the coordinates?

# Alternatively we can aggregate data at the QUADRICULA level & fit a binomial count model
# He binds many obs inside a bigger quadricula

larks2 <- larks |>
  mutate(crestlark = as.numeric(as.character(crestlark)),  # to numeric
         linnet = as.numeric(as.character(linnet)),
         tet_n = rep(1, nrow(larks)),                      # counter for how many grid cells (tetrads) we sum
         N = rep(1, nrow(larks)),                          # number of observations, 1 per row currently
         N = if_else(is.na(crestlark), NA_real_, N)) |>   # set N to NA if observation taken
  group_by(QUADRICULA) |>                                 # group by the larger grid square
  summarise(across(c(N, crestlark, linnet, tet_n, e, n),
                   ~ sum(., na.rm = TRUE))) |>            # sum all needed variables
  mutate(e = e / tet_n, n = n / tet_n)                     # rescale to get avg E,N coords for QUADRICULA

## fit binomial GAM
crest2 <- gam(cbind(crestlark, N - crestlark) ~ s(e, n, k = 100),
              data = larks2, family = binomial, method = "REML")
summary(crest2)


crest3 <- gam(cbind(crestlark, N - crestlark) ~
                s(e, n, k = 100),
              data = larks2, family = quasibinomial,
              method = "REML")
crest3$scale # should be == 1
# check if bionomial is ok

ggplot(tibble(Fitted = fitted(crest2),
              Resid = resid(crest2)),
       aes(Fitted, Resid)) + geom_point()

## South Pole CO2 ----

south <- read_csv(here("data", "south_pole.csv"), col_types = "ddd")
south

ggplot(south, aes(x = c.month, y = co2)) + geom_line()

# Fit naive GAM

m_co2 <- gam(co2 ~ s(c.month, k = 300, bs = "cr"), data = south, method = "REML")
summary(m_co2) # tendency + seasonality

# Predict the next 36 months

new_df <- with(south, tibble(c.month = 1:(nrow(south) + 36)))
fv <- fitted_values(m_co2, data = new_df, scale = "response")
fv

ggplot(fv, aes(x = c.month, y = .fitted)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) +
  geom_line(data = south, aes(c.month, co2), col = "red") +
  geom_line(alpha = 0.4) # this model is not good because it doesn't isolate tendency of seasonality

# Decompose into a 1. seasonal smooth 2. long term trend

m2_co2 <- gam(co2 ~ s(month, bs = "cc") + s(c.month, bs = "cr", k = 300),
              data = south, method = "REML",
              knots = list(month = c(0.5, 12.5)))

summary(m2_co2) # s(month, bs = "cc") captures seasonality
# s(c.month) captures tendency on the long term

# predict

nr <- nrow(south)
new_df <- with(south,
               tibble(c.month = 1:(nr + 36),
                      month = rep(seq_len(12), length = nr + 36)))
fv2 <- fitted_values(m2_co2, data = new_df, scale = "response")
fv2

ggplot(fv2, aes(x = c.month, y = .fitted)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) +
  geom_line(data = south, aes(c.month, co2), col = "red") +
  geom_line(alpha = 0.4)

# A GAM is just a fancy GLM

## !! The type of smoother is controlled by the bs argument (think basis) !!


