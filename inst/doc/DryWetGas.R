## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = FALSE-----
library(Rmbal)
library(Rrelperm)
library(pracma)
library(minpack.lm)
library(ggplot2)
library(dplyr)
library(magrittr)

p_pvt <- c(1798, 1680, 1540, 1428, 1335)                 # psia

Bg <- c(0.00152, 0.00163, 0.00179, 0.00196, 0.00210)     # RB/SCF

Bo <- rep(0, length(p_pvt))    # RB/STB

Rv <- rep(0, length(p_pvt))    # STB/SCF

Rs <- rep(1e8, length(p_pvt))  # SCF/STB

Bw <- rep(1, length(p_pvt))    # RB/STB

muo <- rep(0.5, length(p_pvt))              # cp

muw <- rep(0.25, length(p_pvt))             # cp

mug <- rep(0.015, length(p_pvt))            # cp

pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg,
                             Bw = Bw, muo = muo, mug = mug, muw = muw)

p <- c(1798, 1680, 1540, 1428, 1335)

Np <- rep(0, length.out = length(p))

Gp <- c(0, 0.96, 2.12, 3.21, 3.92) * 1e9    # SCF

Wp <- rep(0, length.out = length(p))

We <- rep(0, length.out = length(p))

Wi <- rep(0, length.out = length(p))

wf <- rep(1, length.out = length(p))

mbal_optim_gas_lst <- mbal_optim_param_gas(input_unit = "Field", output_unit = "Field",  
                                          unknown_param = "G", aquifer_model = NULL, 
                                          G = NULL, phi = 0.13, swi = 0.52, 
                                          Np = Np, Gp = Gp, Wp = Wp, Wi = Wi, 
                                          We = We, pd = 0, p = p, pvt = pvt_table, 
                                          M = 0, cf = 3e-6, wf = wf, sgrw = 0.0)

time_lst <- mbal_time(c(0, 0.5, 1.0, 1.5, 2.0), "year")

# a number of plots will be automatically generated for quality check

optim_results <- mbal_optim_gas(mbal_optim_gas_lst, time_lst)

glimpse(optim_results)


## ----echo=TRUE, fig.align="center", fig.height=4, fig.width=6, warning=FALSE----

mbal_results <- mbal_perform_gas(optim_results, time_lst)

mbal_results

p1 <- mbal_results %>% ggplot(aes(`P (psia)`, `RF_gas`, color = "Forecast")) +
  geom_point(size = 3) +
  scale_color_manual(name="Data", values=c("Forecast" = "red")) +
  ggtitle("Gas Recovery Plot") +
  theme_bw()

p1


p2 <- mbal_results  %>%
  tidyr::pivot_longer(cols = c("Igd", "Ifwd"), names_to = "Drive Mechanism",
                      values_to = "Fraction", values_drop_na = TRUE) %>%
  ggplot(aes(`P (psia)`, Fraction, fill = `Drive Mechanism`)) +
  geom_area() +
  ggtitle("Energy Plot") +
  theme_bw()

p2

## ----echo=TRUE, fig.align="center", fig.height=4, fig.width=6, warning=FALSE----

# gas saturation above the bubble point is zero, however the mbal_forecast_param_oil()
# requires a table of relative permeabilities as an input for the gas-oil system. 
# Therefore, an arbitrary table is generated using the 'Rrelperm' package.
# The generated table does not impact the predictions above the bubble point.


rel_perm <- as.data.frame(Rrelperm::kr2p_gl(SWCON = 0.52, SOIRG = 0.15, SORG = 0.15, 
                                            SGCON = 0.05, SGCRIT = 0.05, KRGCL = 1, 
                                            KROGCG = 1, NG = 1, NOG = 1.0, NP = 101))

colnames(rel_perm) <- c("Sg", "Sl", "Krg", "Krog")

forecast_lst <- mbal_forecast_param_gas(input_unit = "Field", output_unit = "Field",
                                        G = 1.41e10, phi = 0.13, swi = 0.52, pd = 0, 
                                        p = p, pvt = pvt_table, M = 0, cf = 3e-6, 
                                        wf = wf, rel_perm = rel_perm)

glimpse(forecast_lst)

forecast_results <- mbal_forecast_gas(forecast_lst, time_lst)
 
forecast_results

p1 <- forecast_results %>% ggplot(aes(`P (psia)`, `RF_gas`, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = mbal_results, aes(`P (psia)`, `RF_gas`, color = "Field"))+
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("Oil Recovery Plot") +
  theme_bw()

p1


p2 <- forecast_results  %>% 
  tidyr::pivot_longer(cols = c("Igd", "Ifwd"), names_to = "Drive Mechanism", 
                      values_to = "Fraction", values_drop_na = TRUE) %>%
  ggplot(aes(`P (psia)`, Fraction, fill = `Drive Mechanism`)) +
  geom_area() +
  ggtitle("Energy Plot") +
  theme_bw()

p2


## ----fig.align="center", fig.height=4, fig.width=6, warning=FALSE-------------
library(Rmbal)
library(Rrelperm)
library(pracma)
library(minpack.lm)
library(ggplot2)
library(dplyr)
library(magrittr)

p_pvt <- c(6411, 5947, 5509, 5093, 4697, 4319, 3957, 3610, 3276, 2953, 2638)

Bg <- c(6279, 6587, 6933, 7327, 7778, 8300, 8910, 9628, 10487, 11532, 
        12829) / 10000   # RB/SCF

Bo <- rep(0, length(p_pvt))    # RB/STB

Rv <- rep(0, length(p_pvt))    # STB/SCF

Rs <- rep(1e8, length(p_pvt))  # SCF/STB

Bw <- c(1.0452, 1.0467, 1.0480, 1.0493, 1.0506, 1.0517, 1.0529, 1.0540, 1.0551, 
        1.0560, 1.0571)        # RB/STB

muo <- rep(0.5, length(p_pvt))               # cp

muw <- rep(0.25, length(p_pvt))              # cp

mug <- rep(0.025, length(p_pvt))             # cp

pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg,
                             Bw = Bw, muo = muo, mug = mug, muw = muw)

p <- c(6411, 5947, 5509, 5093, 4697, 4319, 3957, 3610, 3276, 2953, 2638)   # psia

Np <- rep(0, length.out = length(p))

Gp <- c(0, 5.475, 10.950, 16.425, 21.900, 27.375, 32.850, 38.325, 43.800, 49.275, 
        54.750) * 1e9                                                      # SCF

Wp <- c(0, 378, 1434, 3056, 5284, 8183, 11864, 16425, 22019, 28860, 37256) # STB

Wi <- rep(0, length.out = length(p))

wf <- rep(1, length.out = length(p))

mbal_optim_gas_lst <- mbal_optim_param_gas(input_unit = "Field", output_unit = "Field",  
                                          unknown_param = "We", aquifer_model = "pot", 
                                          G = 101e9, phi = 0.15, swi = 0.15, 
                                          Np = Np, Gp = Gp, Wp = Wp, Wi = Wi, 
                                          We = NULL, pd = 0, p = p, pvt = pvt_table, 
                                          M = 0, cf = 6e-6, wf = wf, sgrw = 0.0, 
                                          phi_a = 0.15, h_a = 200, r_a = 500, 
                                          r_R = 210.6, tetha = 360, cw_a = 3e-6, 
                                          cf_a = 6e-6, mult_len = 2)

time_lst <- mbal_time(c(1:length(p)), "year")

# a number of plots will be automatically generated for quality check

optim_results <- mbal_optim_gas(mbal_optim_gas_lst, time_lst)

glimpse(optim_results)


## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = FALSE-----

mbal_results <- mbal_perform_gas(optim_results, time_lst)

mbal_results

p1 <- mbal_results %>% ggplot(aes(`P (psia)`, `RF_gas`, color = "Forecast")) +
  geom_point(size = 3) +
  scale_color_manual(name="Data", values=c("Forecast" = "red")) +
  ggtitle("Gas Recovery Plot") +
  theme_bw()

p1


p2 <- mbal_results  %>%
  tidyr::pivot_longer(cols = c("Inwd", "Igd", "Ifwd"), names_to = "Drive Mechanism",
                      values_to = "Fraction", values_drop_na = TRUE) %>%
  ggplot(aes(`P (psia)`, Fraction, fill = `Drive Mechanism`)) +
  geom_area() +
  ggtitle("Energy Plot") +
  theme_bw()

p2

## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = FALSE-----
library(Rmbal)
library(Rrelperm)
library(pracma)
library(minpack.lm)
library(ggplot2)
library(dplyr)
library(magrittr)

p_pvt <- c(10300, 9800, 9300, 8800, 8300, 7800, 7300, 6800, 6300, 5800, 5300, 4800,
         4300, 3800, 3300, 2800, 2300, 1800, 1300, 1050, 800, 738, 300) # psia

Bg <- c(0.5, 0.51, 0.52, 0.532, 0.545, 0.560, 0.577, 0.597, 0.621, 0.648, 0.668,
      0.72, 0.77, 0.83, 0.92, 1.06, 1.27, 1.63, 2.26, 2.925, 3.59,
      4.19, 8.39) / 1000                                                # RB/SCF

Bo <- c(18.62, 18.97, 19.36, 19.8, 20.20, 20.85, 21.49, 22.23, 23.11, 24.15, 25.31,
      26.80, 28.66, 30.89, 34.24, 39.45, 47.26, 60.66, 84.11, 108.86, 133.60, 
      155.93, 312.24)                                                   # RB/STB

Rv <- rep(26.9, length(p_pvt)) / 1e6                                    # STB/SCF

Rs <- rep(37216, length(p_pvt))                                         # SCF/STB

cw <- 2e-6

Bwi <- 1.0

Bw <- Bwi * exp(cw * (p_pvt[1] - p_pvt))       # RB/STB

muo <- rep(0.5, length(p_pvt))                 # cp

muw <- rep(0.25, length(p_pvt))                # cp

mug <- c(0.047, 0.0455, 0.044, 0.0425, 0.041, 0.00395, 0.0379, 0.0363, 0.0347, 
         0.033, 0.0311, 0.0289, 0.02267, 0.0243, 0.0220, 0.0197, 0.0177, 0.0160, 
         0.0147, 0.0142, 0.0138, 0.0136, 0.0128)      # cp

pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg,
                             Bw = Bw, muo = muo, mug = mug, muw = muw)

p <- p_pvt

wf <- rep(1,length(p))

# in-situ oil saturation is zero, however the mbal_forecast_param_gas()
# requires a table of relative permeabilities as an input for the gas-oil system. 
# Therefore, an arbitrary table is generated using the 'Rrelperm' package.
# The generated table does not impact the predictions.

rel_perm <- as.data.frame(Rrelperm::kr2p_gl(SWCON = 0.2, SOIRG = 0.15, SORG = 0.15, 
                                            SGCON = 0.05, SGCRIT = 0.05, KRGCL = 1, 
                                            KROGCG = 1, NG = 2.0, NOG = 2.0, 
                                            NP = 101))

colnames(rel_perm) <- c("Sg", "Sl", "Krg", "Krog")


forecast_lst <- mbal_forecast_param_gas(input_unit = "Field", output_unit = "Field",
                                          G = 69.48e9, phi = 0.1, swi = 0.2, pd = 0, 
                                          p = p, pvt = pvt_table, cf = 3e-6,
                                          M = 0, wf = wf, rel_perm = rel_perm)

time_lst <- mbal_time(1:length(p_pvt), "year")

glimpse(forecast_lst)
 
forecast_results <- mbal_forecast_gas(forecast_lst, time_lst)

forecast_results

p1 <- forecast_results %>% ggplot(aes(`P (psia)`, `RF_gas`, color = "Forecast")) +
  geom_point(size = 3) +
  scale_color_manual(name="Data", values=c("Forecast" = "red")) +
  ggtitle("Gas Recovery Plot") +
  theme_bw()

p1


p2 <- forecast_results  %>%
  tidyr::pivot_longer(cols = c("Igd", "Ifwd"), names_to = "Drive Mechanism",
                      values_to = "Fraction", values_drop_na = TRUE) %>%
  ggplot(aes(`P (psia)`, Fraction, fill = `Drive Mechanism`)) +
  geom_area() +
  ggtitle("Energy Plot") +
  theme_bw()

p2


