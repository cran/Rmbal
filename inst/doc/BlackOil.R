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

 pvt_table_oil <- as.data.frame(Rpvt::pvt_oil("Field", "Field", "black_oil", 
                                              "Standing", "Beggs_Robinson", 200, 
                                              7150, 43, 0.7, c(0,0,0), pb = 4500, 
                                              warning = "no"))

 colnames(pvt_table_oil) <- c("t", "p", "Rs", "Bo", "dens_o", "co", "muo", "Z", 
                              "Bg", "dens_g", "cg", "mug", "m_p")

 pvt_table_water <- as.data.frame(Rpvt::pvt_water("Field", "Field", "water", 
                                                  "Spivey", "Spivey", 200, 7150, 
                                                  0, "no", "no"))

 colnames(pvt_table_water) <- c("t", "p", "Rsw", "Bw", "dens_w", "cw", "muw")

 pvt_table <- dplyr::left_join(pvt_table_oil, pvt_table_water, by = c("p", "t"))

 pvt_table$Rv <- rep(0, length.out = nrow(pvt_table))           # zero for black oil

 pvt_table <- pvt_table %>% dplyr::select(p, Bo, Rs, Rv, Bg, Bw, muo, mug, muw)
 
 p <- c(7150,6600,5800,4950,4500)

 We <- rep(0, length.out = length(p))

 Np <- c(0, 8.072, 22.549, 36.369, 43.473) * 1e6

 Rp <- c(0, rep(pvt_table$Rs[nrow(pvt_table)], 4))

 Wp <- rep(0, length.out = length(p))

 Wi <- rep(0, length.out = length(p))

 Gi <- rep(0, length.out = length(p))

 wf <- rep(1,length(p))

 mbal_optim_oil_lst <- mbal_optim_param_oil(input_unit = "Field", output_unit = "Field",  
                                            unknown_param = "N", aquifer_model = NULL, 
                                            m = 0, phi = 0.2, swi = 0.43, Np = Np, 
                                            Rp = Rp, Wp = Wp, Gi = Gi, Wi = Wi, We = We, 
                                            pb = 4500, p = p, pvt = pvt_table, cf = 4e-6, 
                                            wf = wf, sorg = 0, sorw = 0)

 time_lst <- mbal_time(1:5, "year")

 # a number of plots will be automatically generated for quality check
 
 optim_results <- mbal_optim_oil(mbal_optim_oil_lst, time_lst)
 
 glimpse(optim_results)

## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = FALSE-----

mbal_results <- mbal_perform_oil(optim_results, time_lst)

mbal_results

## ----echo=TRUE, fig.align="center", fig.height=4, fig.width=6, warning=FALSE----

# gas saturation above the bubble point is zero, however the mbal_forecast_param_oil()
# requires a table of relative permeabilities as an input for the gas-oil system. 
# Therefore, an arbitrary table is generated using the 'Rrelperm' package.
# The generated table does not impact the predictions above the bubble point.


rel_perm <- as.data.frame(Rrelperm::kr2p_gl(SWCON = 0.43, SOIRG = 0.15, SORG = 0.15, 
                                            SGCON = 0.05, SGCRIT = 0.05, KRGCL = 1, 
                                            KROGCG = 1, NG = 1.5, NOG = 1.0, NP = 101))

colnames(rel_perm) <- c("Sg", "Sl", "Krg", "Krog")

forecast_lst <- mbal_forecast_param_oil(input_unit = "Field", output_unit = "Field",
                                        N = 6.35e8, m = 0, phi = 0.1, swi = 0.43,
                                        Gi = Gi, pb = 4500, p = p, pvt = pvt_table,
                                        cf = 4e-6, wf = wf, sorg = 0, rel_perm = rel_perm)

glimpse(forecast_lst)

forecast_results <- mbal_forecast_oil(forecast_lst, time_lst)
 
forecast_results

p1 <- forecast_results %>% ggplot(aes(`P (psia)`, `RF_oil`, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = mbal_results, aes(`P (psia)`, `RF_oil`, color = "Field"))+
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("Oil Recovery Plot") +
  theme_bw()

p1


p2 <- forecast_results  %>% 
  tidyr::pivot_longer(cols = c("Isd", "Ifwd"), names_to = "Drive Mechanism", 
                      values_to = "Fraction", values_drop_na = TRUE) %>%
  ggplot(aes(`P (psia)`, Fraction, fill = `Drive Mechanism`)) +
  geom_area() +
  ggtitle("Energy Plot") +
  theme_bw()

p2


## ----echo=TRUE, fig.align="center", fig.height=4, fig.width=6, warning=FALSE----
library(Rmbal)
library(Rrelperm)
library(pracma)
library(ggplot2)
library(dplyr)
library(magrittr)

 pvt_table_oil <- as.data.frame(Rpvt::pvt_oil("Field", "Field", "black_oil", 
                                              "Standing", "Beggs_Robinson", 200, 
                                              7150, 43, 0.7, c(0,0,0), pb = 4500, 
                                              warning = "no"))

 colnames(pvt_table_oil) <- c("t", "p", "Rs", "Bo", "dens_o", "co", "muo", "Z", 
                              "Bg", "dens_g", "cg", "mug", "m_p")

 pvt_table_water <- as.data.frame(Rpvt::pvt_water("Field", "Field", "water", 
                                                  "Spivey", "Spivey", 200, 7150, 
                                                  3e-6, "no", "no"))

 colnames(pvt_table_water) <- c("t", "p", "Rsw", "Bw", "dens_w", "cw", "muw")

 pvt_table <- dplyr::left_join(pvt_table_oil, pvt_table_water, by = c("p", "t"))
 
 pvt_table$Rv <- 0.0          # zero for black oil

 pvt_table <- pvt_table %>% dplyr::select(p, Bo, Rs, Rv, Bg, Bw, muo, mug, muw)
 
 p <- c(7150, 6600, 5800, 4950, 4500, 4350, 4060, 3840, 3600, 3480, 3260, 3100, 2940, 2800)

 We <- rep(0, length.out = length(p))

 Np <- c(0, 8.072, 22.549, 36.369, 43.473, 49.182, 58.383, 64.812, 69.562, 74.572, 
         78.4, 81.275, 83.879, 86.401) * 1e6

 Rp <- c(0, rep(pvt_table$Rs[nrow(pvt_table)], 4), 1576, 1788, 1992, 2158, 2383, 2596, 2785, 2953, 3103)

 Wp <- rep(0, length.out = length(p))

 Wi <- rep(0, length.out = length(p))

 Gi <- rep(0, length.out = length(p))

 wf <- rep(1,length(p))

 mbal_optim_oil_lst <- mbal_optim_param_oil(input_unit = "Field", output_unit = "Field",  
                                            unknown_param = "N", aquifer_model = NULL, 
                                            m = 0, phi = 0.2, swi = 0.43, Np = Np, 
                                            Rp = Rp, Wp = Wp, Gi = Gi, Wi = Wi, We = We, 
                                            pb = 4500, p = p, pvt = pvt_table, cf = 4e-6, 
                                            wf = wf, sorg = 0, sorw = 0)

 glimpse(mbal_optim_oil_lst)
 
 
 time_lst <- mbal_time(1:14, "year")

 # a number of plots will be automatically generated for quality check
 
 optim_results <- mbal_optim_oil(mbal_optim_oil_lst, time_lst)
 
 glimpse(optim_results)

## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = FALSE-----

mbal_results <- mbal_perform_oil(optim_results, time_lst)

mbal_results

## ----echo=TRUE, fig.align="center", fig.height=4, fig.width=6, warning=FALSE----

# Step I: generating a set of pseudo relative permeability curves by minimizing 
# the difference between field and theoretical 'fg' values
# 'swcrit' and 'sgcrit' values are estimated from the 'mbal_results' data frame

fun_fg <- function(x, swcrit, sgcrit, muo, mug, sg, fg, krg_kro) {
  Kr_table <- Rrelperm::kr2p_gl(SWCON = swcrit, SOIRG = x[1], SORG = x[1], 
                                            SGCON = sgcrit, SGCRIT = sgcrit, KRGCL = x[2], 
                                            KROGCG = x[3], NG = x[4], NOG = x[5], NP = 101)
  l <- length(fg)
  krg_est_sub <- vector(length = l)
  kro_est_sub <- vector(length = l)
  krg_est_sub <- approx(x = Kr_table[,1], y = Kr_table[,3], xout = sg, rule = 2)$y
  kro_est_sub <- approx(x = Kr_table[,1], y = Kr_table[,4], xout = sg, rule = 2)$y
  fg_est <- (krg_est_sub / mug) / (krg_est_sub / mug + kro_est_sub / muo)
  krg_kro_est <- krg_est_sub / kro_est_sub
  error <- (fg - fg_est) ^ 2 + (krg_kro - krg_kro_est) ^ 2
  return(error)
}

swcrit <- 0.44
sgcrit <- 0.015
p <- mbal_results$`P (psia)`
sg <- mbal_results$SGo  # gas saturation in the oil leg
fg <- mbal_results$fg   # in-situ gas fractional flow
krg_kro <- mbal_results$`krg/kro`
muo <- approx(x = pvt_table$p, y = pvt_table$muo, xout = p, rule = 2)$y
mug <- approx(x = pvt_table$p, y = pvt_table$mug, xout = p, rule = 2)$y

par <- c(0.1, 1, 1, 2, 2)
lower = c(0, 0.1, 0.1, 0.1, 0.1)
upper = c(1 - swcrit - sgcrit, 1.0, 1.0, 10.0, 10.0)
opt_results <- minpack.lm::nls.lm(par = par, fn = fun_fg, swcrit = swcrit, sgcrit = sgcrit, 
                                  muo = muo, mug = mug, sg = sg, fg = fg, krg_kro = krg_kro, 
                                  lower = lower, upper = upper)

opt_results

sol <- opt_results$par

sol

rel_perm <- as.data.frame(Rrelperm::kr2p_gl(SWCON = swcrit, SOIRG = sol[1], SORG = sol[1], 
                                            SGCON = sgcrit, SGCRIT = sgcrit, KRGCL = sol[2], 
                                            KROGCG = sol[3], NG = sol[4], NOG = sol[5], 
                                            NP = 101))

colnames(rel_perm) <- c("Sg", "Sl", "Krg", "Krog")

krg_est <- approx(x = rel_perm[,1], y = rel_perm[,3], xout = sg, rule = 2)$y

kro_est <- approx(x = rel_perm[,1], y = rel_perm[,4], xout = sg, rule = 2)$y

fg_est <- (krg_est/ mug) / (krg_est / mug + kro_est / muo)
  
fg_df <- data.frame(Sg = sg)

fg_df$fg <- fg_est 

p_fg <- fg_df %>% ggplot(aes(sg, fg, color = "Model")) +
  geom_point(size = 3) +
  geom_point(data = mbal_results, aes(SGo, fg, color = "Field"), size = 3)+
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  scale_color_manual(name="Data", values=c("Model" = "red", "Field" = "black")) +
  ggtitle("Gas Frational Flow Plot") +
  theme_bw()

p_fg 

p_forecast <- c(p, 2600, 2400, 2200, 2000, 1800, 1600, 1400, 1200)

Gi_forecast <- c(Gi, 0, 0, 0, 0, 0, 0, 0, 0)

wf_forecast <- c(wf, 1, 1, 1, 1, 1, 1, 1, 1)

time_lst_forecast <-  mbal_time(1:22, "year")

forecast_lst <- mbal_forecast_param_oil(input_unit = "Field", output_unit = "Field",
                                        N = 6.35e8, m = 0, phi = 0.1, swi = 0.43,
                                        Gi = Gi_forecast, pb = 4500, p = p_forecast, pvt = pvt_table,
                                        cf = 4e-6, wf = wf_forecast, sorg = 0, rel_perm = rel_perm)

glimpse(forecast_lst)



forecast_results <- mbal_forecast_oil(forecast_lst, time_lst_forecast)
 
forecast_results

p1 <- forecast_results %>% ggplot(aes(`P (psia)`, SOo, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = mbal_results, aes(`P (psia)`, SOo, color = "Field"))+
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("Oil Saturation Plot") +
  theme_bw()

p1


p2 <- forecast_results %>% ggplot(aes(`P (psia)`, `GOR (SCF/STB)`, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = mbal_results, aes(`P (psia)`, `GOR (SCF/STB)`, color = "Field"))+
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("GOR Plot") +
  theme_bw()

p2


p3 <- forecast_results %>% ggplot(aes(`P (psia)`, `RF_oil`, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = mbal_results, aes(`P (psia)`, `RF_oil`, color = "Field"))+
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("Oil Recovery Plot") +
  theme_bw()

p3


p4 <- forecast_results %>% ggplot(aes(`P (psia)`, `Liq_volume`, color = "Forecast")) +
  geom_point(size = 3) +
  scale_color_manual(name="Data", values=c("Forecast" = "red")) +
  ggtitle("CCE Liquid Volume Plot") +
  theme_bw()

p4


p5 <- forecast_results  %>% 
  tidyr::pivot_longer(cols = c("Isd", "Ifwd"), names_to = "Drive Mechanism", 
                      values_to = "Fraction", values_drop_na = TRUE) %>%
  ggplot(aes(`P (psia)`, Fraction, fill = `Drive Mechanism`)) +
  geom_area() +
  ggtitle("Energy Plot") +
  theme_bw()

p5

## ----echo=TRUE, fig.align="center", fig.height=4, fig.width=6, warning=FALSE----
library(Rmbal)
library(Rrelperm)
library(pracma)
library(ggplot2)
library(dplyr)
library(magrittr)

p_pvt <- c(3330, 3150, 3000, 2850, 2700, 2550, 2400)

Bo <- c(1.2511, 1.2353, 1.2222, 1.2122, 1.2022, 1.1922, 1.1822)

Rs <- c(510, 477, 450, 425, 401, 375, 352)

Bg <- c(0.00087, 0.00092, 0.00096, 0.00101, 0.00107, 0.00113, 0.00120)

cw <- 2e-6

Bwi <- 1.0

Bw <- Bwi * exp(cw * (p_pvt[1] - p_pvt))

Rv <- rep(0, length(p_pvt))

muo <- rep(0.5, length(p_pvt))

muw <- rep(0.25, length(p_pvt))

mug <- rep(0.02, length(p_pvt))

pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg, 
                               Bw = Bw, muo = muo, mug = mug, muw = muw)

p <- c(3330, 3150, 3000, 2850, 2700, 2550, 2400)

We <- rep(0, length.out = length(p))

Np <- c(0, 3.295, 5.903, 8.852, 11.503, 14.513, 17.730) * 1e6

Rp <- c(0, 1050, 1060, 1160, 1235, 1265, 1300)

Wp <- rep(0, length.out = length(p))

Wi <- rep(0, length.out = length(p))

Gi <- rep(0, length.out = length(p))

wf <- c(1, 1, 1, 0, 1, 0, 1)

mbal_optim_oil_lst <- mbal_optim_param_oil(input_unit = "Field", output_unit = "Field",  
                                          unknown_param = "N_m", aquifer_model = NULL, 
                                          phi = 0.2, swi = 0.2, Np = Np, 
                                          Rp = Rp, Wp = Wp, Gi = Gi, Wi = Wi, We = We, 
                                          pb = 3330, p = p, pvt = pvt_table, cf = 0, 
                                          wf = wf, sorg = 0.2, sorw = 0)

time_lst <- mbal_time(c(0, 365, 730, 1095, 1460, 1825, 2190), "day")

# a number of plots will be automatically generated for quality check

optim_results <- mbal_optim_oil(mbal_optim_oil_lst, time_lst)

glimpse(optim_results)

## ----echo=TRUE, fig.align="center", fig.height=4, fig.width=6, warning=FALSE----

mbal_results <- mbal_perform_oil(optim_results, time_lst)

mbal_results

## ----echo=TRUE, fig.align="center", fig.height=4, fig.width=6, warning=FALSE----

# Step I: generating a set of pseudo relative permeability curves by minimizing 
# the difference between field and theoretical 'fg' values
# 'swcrit' and 'sgcrit' values are estimated from the 'mbal_results' data frame

fun_fg <- function(x, swcrit, sgcrit, muo, mug, sg, fg, krg_kro) {
  Kr_table <- Rrelperm::kr2p_gl(SWCON = swcrit, SOIRG = x[1], SORG = x[1], 
                                            SGCON = sgcrit, SGCRIT = sgcrit, KRGCL = x[2], 
                                            KROGCG = x[3], NG = x[4], NOG = x[5], NP = 101)
  l <- length(fg)
  krg_est_sub <- vector(length = l)
  kro_est_sub <- vector(length = l)
  krg_est_sub <- approx(x = Kr_table[,1], y = Kr_table[,3], xout = sg, rule = 2)$y
  kro_est_sub <- approx(x = Kr_table[,1], y = Kr_table[,4], xout = sg, rule = 2)$y
  fg_est <- (krg_est_sub / mug) / (krg_est_sub / mug + kro_est_sub / muo)
  krg_kro_est <- krg_est_sub / kro_est_sub
  error <- (fg - fg_est) ^ 2 + (krg_kro - krg_kro_est) ^ 2
  return(error)
}

swcrit <- 0.2
sgcrit <- 0.00
p <- mbal_results$`P (psia)`
sg <- mbal_results$SGo  # gas saturation in the oil leg
fg <- mbal_results$fg   # in-situ gas fractional flow
krg_kro <- mbal_results$`krg/kro`
muo <- approx(x = pvt_table$p, y = pvt_table$muo, xout = p, rule = 2)$y
mug <- approx(x = pvt_table$p, y = pvt_table$mug, xout = p, rule = 2)$y

par <- c(0.1, 0.3, 1.0, 3, 3)
lower = c(0, 0.3, 1.0, 0.1, 0.1)
upper = c(1 - swcrit - sgcrit, 0.3, 1.0, 10.0, 10.0)
opt_results <- minpack.lm::nls.lm(par = par, fn = fun_fg, swcrit = swcrit, sgcrit = sgcrit, 
                                  muo = muo, mug = mug, sg = sg, fg = fg, krg_kro = krg_kro, 
                                  lower = lower, upper = upper)

opt_results

sol <- opt_results$par

sol

rel_perm <- as.data.frame(Rrelperm::kr2p_gl(SWCON = swcrit, SOIRG = sol[1], SORG = sol[1], 
                                            SGCON = sgcrit, SGCRIT = sgcrit, KRGCL = sol[2], 
                                            KROGCG = sol[3], NG = sol[4], NOG = sol[5], 
                                            NP = 101))

colnames(rel_perm) <- c("Sg", "Sl", "Krg", "Krog")

krg_est <- approx(x = rel_perm[,1], y = rel_perm[,3], xout = sg, rule = 2)$y

kro_est <- approx(x = rel_perm[,1], y = rel_perm[,4], xout = sg, rule = 2)$y

fg_est <- (krg_est/ mug) / (krg_est / mug + kro_est / muo)
  
fg_df <- data.frame(Sg = sg)

fg_df$fg <- fg_est 

p_fg <- fg_df %>% ggplot(aes(sg, fg, color = "Model")) +
  geom_point(size = 3) +
  geom_point(data = mbal_results, aes(SGo, fg, color = "Field"), size = 3)+
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  scale_color_manual(name="Data", values=c("Model" = "red", "Field" = "black")) +
  ggtitle("Gas Frational Flow Plot") +
  theme_bw()

p_fg 

p_forecast <- p_pvt

Gi_forecast <- c(Gi)

wf_forecast <- c(wf)

time_lst_forecast <- mbal_time(c(0, 365, 730, 1095, 1460, 1825, 2190), "day")

forecast_lst <- mbal_forecast_param_oil(input_unit = "Field", output_unit = "Field",
                                        N = 1.37e8, m = 0.377, phi = 0.2, swi = 0.2,
                                        Gi = Gi_forecast, pb = 3330, p = p_forecast, pvt = pvt_table,
                                        cf = 0, wf = wf_forecast, sorg = 0.2, rel_perm = rel_perm)

glimpse(forecast_lst)


forecast_results <- mbal_forecast_oil(forecast_lst, time_lst_forecast)

forecast_results

p1 <- forecast_results %>% ggplot(aes(`P (psia)`, SOo, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = mbal_results, aes(`P (psia)`, SOo, color = "Field"))+
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("Oil Saturation Plot") +
  theme_bw()

p1


p2 <- forecast_results %>% ggplot(aes(`P (psia)`, `GOR (SCF/STB)`, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = mbal_results, aes(`P (psia)`, `GOR (SCF/STB)`, color = "Field"))+
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("GOR Plot") +
  theme_bw()

p2


p3 <- forecast_results %>% ggplot(aes(`P (psia)`, `RF_oil`, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = mbal_results, aes(`P (psia)`, `RF_oil`, color = "Field"))+
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("Oil Recovery Plot") +
  theme_bw()

p3


p4 <- forecast_results %>% ggplot(aes(`P (psia)`, `Liq_volume`, color = "Forecast")) +
  geom_point(size = 3) +
  scale_color_manual(name="Data", values=c("Forecast" = "red")) +
  ggtitle("CCE Liquid Volume Plot") +
  theme_bw()

p4


p5 <- forecast_results  %>%
  tidyr::pivot_longer(cols = c("Igd", "Isd"), names_to = "Drive Mechanism",
                      values_to = "Fraction", values_drop_na = TRUE) %>%
  ggplot(aes(`P (psia)`, Fraction, fill = `Drive Mechanism`)) +
  geom_area() +
  ggtitle("Energy Plot") +
  theme_bw()

p5

## ----echo=TRUE, fig.align="center", fig.height=4, fig.width=6, warning=FALSE----
library(Rmbal)
library(Rrelperm)
library(pracma)
library(ggplot2)
library(dplyr)
library(magrittr)

p_pvt <- c(2740, 2500, 2290, 2109, 1949, 1818, 1702, 1608, 1535, 1480, 1440)

Bo <- c(1.404, 1.374, 1.349, 1.329, 1.316, 1.303, 1.294, 1.287, 1.280, 1.276, 1.273)

Rs <- c(650, 592, 545, 507, 471, 442, 418, 398, 383, 371, 364)

Bg <- c(0.00093, 0.00098, 0.00107, 0.00117, 0.00128, 0.00139, 0.00150, 0.00160, 
        0.00170, 0.00176, 0.00182)

cw <- 3e-6

Bwi <- 1.0

Bw <- Bwi * exp(cw * (p_pvt[1] - p_pvt))

Rv <- rep(0, length(p_pvt))

muo <- rep(0.5, length(p_pvt))

muw <- rep(0.55, length(p_pvt))

mug <- rep(0.02, length(p_pvt))

pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg, 
                               Bw = Bw, muo = muo, mug = mug, muw = muw)

p <- c(2740, 2500, 2290, 2109, 1949, 1818, 1702, 1608, 1535, 1480, 1440)

We <- rep(0, length.out = length(p))

Np <- c(0, 7.88, 18.42, 29.15, 40.69, 50.14, 58.42, 65.39, 70.74, 74.54, 77.43) * 1e6

Rp <- c(0, 760, 845, 920, 975, 1025, 1065, 1095, 1120, 1145, 1160)

Wp <- rep(0, length.out = length(p))

Wi <- rep(0, length.out = length(p))

Gi <- rep(0, length.out = length(p))

wf <- rep(1, length.out = length(p))

mbal_optim_oil_lst <- mbal_optim_param_oil(input_unit = "Field", output_unit = "Field", 
                                           unknown_param = "We", 
                                           aquifer_model = "pss_rad_edge", N = 312e6, 
                                           m = 0, phi = 0.25, swi = 0.25, Np = Np, Rp = Rp, 
                                           Wp = Wp, Gi = Gi, Wi = Wi, We = NULL, pb = 2740, 
                                           p = p, pvt = pvt_table, cf = 4e-6, 
                                           phi_a = 0.25, perm_h_a = 200, h_a = 100, 
                                           r_a = 5 * 9200, r_R = 9200, tetha = 140, 
                                           muw_a = 0.55, cw_a = 3e-6, cf_a = 4e-6, 
                                           wf = wf, sorw = 0.2, sorg = 0, 
                                           mult_len = c(1,1), lower = c(-Inf, 1), 
                                           upper = c(Inf, 1), 
                                           control = list(maxiter = 100), )

time_lst <- mbal_time(0:10, "year")

# a number of plots will be automatically generated for quality check

optim_results <- mbal_optim_oil(mbal_optim_oil_lst, time_lst)

glimpse(optim_results)

## ----echo=TRUE, fig.align="center", fig.height=4, fig.width=6, warning=FALSE----

mbal_results <- mbal_perform_oil(optim_results, time_lst)

mbal_results

p1 <- mbal_results %>% ggplot(aes(`P (psia)`, SOo, color = "Field")) +
  geom_point(size = 3) +
  scale_color_manual(name="Data", values=c("Field" = "black")) +
  ggtitle("Oil Saturation Plot") +
  theme_bw()

p1


p2 <- mbal_results %>% ggplot(aes(`P (psia)`, `GOR (SCF/STB)`, color = "Field")) +
  geom_point(size = 3) +
  scale_color_manual(name="Data", values = c("Field" = "black")) +
  ggtitle("GOR Plot") +
  theme_bw()

p2


p3 <- mbal_results %>% ggplot(aes(`P (psia)`, `RF_oil`, color = "Field")) +
  geom_point(size = 3) +
  scale_color_manual(name="Data", values=c("Field" = "black")) +
  ggtitle("Oil Recovery Plot") +
  theme_bw()

p3


p4 <- mbal_results  %>%
  tidyr::pivot_longer(cols = c("Isd", "Inwd", "Ifwd"), names_to = "Drive Mechanism",
                      values_to = "Fraction", values_drop_na = TRUE) %>%
  ggplot(aes(`P (psia)`, Fraction, fill = `Drive Mechanism`)) +
  geom_area() +
  ggtitle("Energy Plot") +
  theme_bw()

p4

## ----echo=TRUE, fig.align="center", fig.height=4, fig.width=6, warning=FALSE----
library(Rmbal)
library(Rrelperm)
library(pracma)
library(ggplot2)
library(dplyr)
library(magrittr)

p_pvt <- c(2855,	2779,	2627,	2457,	2402,	2223,	2080,	1833,	1665,	1460)  # psia

Bo <- c(1.2665, 1.2677, 1.2681, 1.2554, 1.2512, 1.2383, 1.2278, 1.2074, 1.1949, 
        1.1802)  # RB/STB

Rs <- c(0.501, 0.501,	0.4973,	0.4671,	0.4574,	0.4269,	0.4024,	0.3579,	0.3277,	
        0.2908) * 1000  #SCF/STB

Bg <- c(0.9201,	0.9637, 1.0502, 1.0977, 1.1146, 1.201, 1.2825, 1.4584, 1.6112, 
        1.8526) / 1000  # RB/SCF

Bw <- c(1.0222,	1.0224,	1.0228,	1.0232,	1.0233,	1.0237,	1.024,	1.0246,	1.025,	
        1.0254)  # RB/STB

Rv <- rep(0, length(p_pvt))

muo <- rep(0.5, length(p_pvt))

muw <- rep(0.55, length(p_pvt))

mug <- rep(0.02, length(p_pvt))

pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg, 
                               Bw = Bw, muo = muo, mug = mug, muw = muw)

p <- c(2855,	2779,	2627,	2457,	2402,	2223,	2080,	1833,	1665,	1460)  # psia

Np <- c(0,	192821,	633942,	1314880,	1524400,	2152960,	2572000,	3200560,	
        3584680,	4003720)  # STB

Rp <- c(0,	501,	492.2595442,	540.4827817,	558.2091315,	629.7005054,
        708.8841369,	853.8693229,	947.8502963,	1053.050663)  # SCF/STB

Wp <- rep(0, length.out = length(p))

Wi <- rep(0, length.out = length(p))

Gi <- rep(0, length.out = length(p))

wf <- rep(1, length.out = length(p))

mbal_optim_oil_lst <- mbal_optim_param_oil(input_unit = "Field", output_unit = "Field",
                                           unknown_param = "We", 
                                           aquifer_model = "pot", N = 20e6, 
                                           m = 0, phi = 0.28, swi = 0.208, Np = Np, 
                                           Rp = Rp, Wp = Wp, Gi = Gi, Wi = Wi, 
                                           We = NULL, pb = 2648, p = p, pvt = pvt_table, 
                                           cf = 26e-6, phi_a = 0.28, h_a = 200, r_a = 1000, 
                                           r_R = 300, tetha = 360, cw_a = 2.28e-6, 
                                           cf_a = 26e-6, mult_len = 3, wf = wf, 
                                           sorw = 0.0, sorg = 0)

time_lst <- mbal_time(c(0,305, 700, 1285, 1465, 2005, 2365, 2905, 3235, 3595), "day")

# a number of plots will be automatically generated for quality check

optim_results <- mbal_optim_oil(mbal_optim_oil_lst, time_lst)

glimpse(optim_results)

## ----echo=TRUE, fig.align="center", fig.height=4, fig.width=6, warning=FALSE----

mbal_results <- mbal_perform_oil(optim_results, time_lst)

mbal_results

p1 <- mbal_results %>% ggplot(aes(`P (psia)`, SOo, color = "Field")) +
  geom_point(size = 3) +
  scale_color_manual(name="Data", values=c("Field" = "black")) +
  ggtitle("Oil Saturation Plot") +
  theme_bw()

p1


p2 <- mbal_results %>% ggplot(aes(`P (psia)`, `GOR (SCF/STB)`, color = "Field")) +
  geom_point(size = 3) +
  scale_color_manual(name="Data", values = c("Field" = "black")) +
  ggtitle("GOR Plot") +
  theme_bw()

p2


p3 <- mbal_results %>% ggplot(aes(`P (psia)`, `RF_oil`, color = "Field")) +
  geom_point(size = 3) +
  scale_color_manual(name="Data", values=c("Field" = "black")) +
  ggtitle("Oil Recovery Plot") +
  theme_bw()

p3


p4 <- mbal_results  %>%
  tidyr::pivot_longer(cols = c("Isd", "Inwd", "Ifwd"), names_to = "Drive Mechanism",
                      values_to = "Fraction", values_drop_na = TRUE) %>%
  ggplot(aes(`P (psia)`, Fraction, fill = `Drive Mechanism`)) +
  geom_area() +
  ggtitle("Energy Plot") +
  theme_bw()

p4

## ----echo=TRUE, fig.align="center", fig.height=4, fig.width=6, warning=FALSE----
library(Rmbal)
library(Rrelperm)
library(pracma)
library(ggplot2)
library(dplyr)
library(magrittr)

p_pvt <- c(2627, 2457, 2402, 2223, 2080, 1833, 1665, 1460, 1400)

Bo <- c(1.2681, 1.2554, 1.2512, 1.2383, 1.2278, 1.2074, 1.1949, 1.1802, 1.17)

Rs <- c(0.501, 0.4973, 0.4671, 0.4574, 0.4269, 0.4024, 0.3579, 0.3277, 0.2908) * 1e3

Bg <- c(1.0502, 1.0977, 1.1146, 1.2010, 1.2825, 1.4584, 1.6112, 1.8526, 2.1) / 1e3

Bw <- c(1.0228, 1.0232, 1.0233, 1.0237, 1.0240, 1.0246, 1.0250, 1.0254, 1.0258)

Rv <- rep(0, length(p_pvt))

muo <- rep(0.5, length(p_pvt))

muw <- rep(0.55, length(p_pvt))

mug <- rep(0.02, length(p_pvt))

pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg, 
                               Bw = Bw, muo = muo, mug = mug, muw = muw)

p <- c(2627, 2457, 2402, 2223, 2080, 1833, 1665, 1460)

We <- rep(0, length.out = length(p))

Np <- c(0, 192821, 633942, 1314880, 1524400, 2152960, 2572000, 3200560)

Rp <- c(0, 490.1593, 492.2595, 540.4828, 558.2091, 629.7005, 708.8841, 853.8693)

Wp <- Wp <- c(0, 0, 0, 4, 7, 26, 60, 822)

Wi <- rep(0, length.out = length(p))

Gi <- rep(0, length.out = length(p))

wf <- c(1, 0, 1, 1, 1, 1, 1, 1)

mbal_optim_oil_lst <- mbal_optim_param_oil(input_unit = "Field", output_unit = "Field", 
                                           unknown_param = "We", 
                                           aquifer_model = "pot", N = 198e5, 
                                           m = 0.05, phi = 0.28, swi = 0.25, Np = Np, Rp = Rp, 
                                           Wp = Wp, Gi = Gi, Wi = Wi, We = NULL, pb = 2627, 
                                           p = p, pvt = pvt_table, cf = 26e-6, 
                                           phi_a = 0.28, perm_h_a = 100, h_a = 100, 
                                           r_a = 2592.013, r_R = 1269.619, tetha = 360, 
                                           muw_a = 0.55, cw_a = 2.28e-6, cf_a = 26e-6, 
                                           wf = wf, sorw = 0.2, sorg = 0, 
                                           mult_len = c(2))

time_lst <- mbal_time(c(0,305, 700, 1285, 1465, 2005, 2365, 2905), "day")

# a number of plots will be automatically generated for quality check

optim_results <- mbal_optim_oil(mbal_optim_oil_lst, time_lst)

glimpse(optim_results)

## ----echo=TRUE, fig.align="center", fig.height=4, fig.width=6, warning=FALSE----

mbal_results <- mbal_perform_oil(optim_results, time_lst)

mbal_results

p1 <- mbal_results %>% ggplot(aes(`P (psia)`, SOo, color = "Field")) +
  geom_point(size = 3) +
  scale_color_manual(name="Data", values=c("Field" = "black")) +
  ggtitle("Oil Saturation Plot") +
  theme_bw()

p1


p2 <- mbal_results %>% ggplot(aes(`P (psia)`, `GOR (SCF/STB)`, color = "Field")) +
  geom_point(size = 3) +
  scale_color_manual(name="Data", values = c("Field" = "black")) +
  ggtitle("GOR Plot") +
  theme_bw()

p2


p3 <- mbal_results %>% ggplot(aes(`P (psia)`, `RF_oil`, color = "Field")) +
  geom_point(size = 3) +
  scale_color_manual(name="Data", values=c("Field" = "black")) +
  ggtitle("Oil Recovery Plot") +
  theme_bw()

p3


p4 <- mbal_results  %>%
  tidyr::pivot_longer(cols = c("Igd","Isd", "Inwd", "Ifwd"), names_to = "Drive Mechanism",
                      values_to = "Fraction", values_drop_na = TRUE) %>%
  ggplot(aes(`P (psia)`, Fraction, fill = `Drive Mechanism`)) +
  geom_area() +
  ggtitle("Energy Plot") +
  theme_bw()

p4

## ----echo=TRUE, fig.align="center", fig.height=4, fig.width=6, warning=FALSE----
library(Rmbal)
library(Rrelperm)
library(pracma)
library(ggplot2)
library(dplyr)
library(magrittr)

p_pvt <- c(2000,	1800,	1700,	1640,	1600,	1400,	1200,	1000,	800, 600, 400, 200) # psia

Bo <- c(1.467, 1.472, 1.475, 1.463, 1.453, 1.408, 1.359, 1.322, 1.278, 1.237,
        1.194, 1.141) # RB/STB

Rs <- c(838.5, 838.5, 838.5, 816.1, 798.4, 713.4, 621, 548, 464, 383.9, 297.4, 
        190.9) #SCF/STB

Bg <- c(1.749, 1.755, 1.758, 1.921, 1.977, 2.308, 2.73, 3.328, 4.163, 5.471, 
        7.786, 13.331) / 1000 # RB/SCF

Rv <- c(1192.6, 1192.6,	1192.6,	0.2, 0.2,	0, 0, 0, 0, 0, 0, 0) / 1e6 # STB/SCF

cw <- 2e-6

Bwi <- 1.0

Bw <- Bwi * exp(cw * (p_pvt[1] - p_pvt))


muo <- c(0.3201,	0.3114,	0.3071,	0.3123,	0.316, 0.34,	0.371, 0.397, 0.432, 0.471,	
         0.518,	0.589)  # cp

muw <- rep(0.25, length(p_pvt))

mug <- c(0.3201, 0.3114, 0.3071, 0.0157, 0.0155, 0.014,	0.0138,	0.0132,	0.0126,	0.0121,
         0.0116,	0.0108) # cp

pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg, 
                               Bw = Bw, muo = muo, mug = mug, muw = muw)

p <- c(2000,	1800,	1700,	1640,	1600,	1400,	1200,	1000,	800, 600, 400, 200) # psia

Gi <- rep(0, length.out = length(p))

wf <- rep(1, length.out = length(p))

# generating a set of pseudo relative permeability curves using 
# laboratory 'Kr' values

sg_lab <- c(0.05, 0.152, 0.248, 0.352, 0.448, 0.552, 0.65)

krg_lab <- c(0, 0.05, 0.09, 0.18, 0.3, 0.5, 1)

kro_lab <- c(1, 0.6, 0.35, 0.13, 0.04, 0.01, 0)

swcrit_lab <- 0.2

sgcrit_lab <- 0.05

sorgr_lab <- 0.15

fun_kr <- function(x, swcrit, sgcrit, sorg, sg, krg, kro) {
   
  kr_est <- Rrelperm::kr2p_gl(SWCON = swcrit, SOIRG = sorg, SORG = sorg, 
                               SGCON = sgcrit, SGCRIT = sgcrit, KRGCL = 1, 
                               KROGCG = 1, NG = x[1], NOG = x[2], NP = 101)
   
   krg_est_sub <- approx(x = kr_est[,1], y = kr_est[,3], xout = sg, rule = 2)$y
   
   kro_est_sub <- approx(x = kr_est[,1], y = kr_est[,4], xout = sg, rule = 2)$y
   
   error <- (krg - krg_est_sub) ^ 2 + (kro - kro_est_sub) ^ 2
   
   return(error)
}

par <- c(2, 2)

opt_results <- minpack.lm::nls.lm(par = par, fn = fun_kr, swcrit = swcrit_lab, 
                                  sgcrit = sgcrit_lab, sorg = sorgr_lab, sg = sg_lab, 
                                  krg = krg_lab, kro = kro_lab, 
                                  lower = c(0.1,0.1), upper = c(10,10))

sol <- opt_results$par

sol

rel_perm <- as.data.frame(Rrelperm::kr2p_gl(SWCON = swcrit_lab, SOIRG = sorgr_lab, 
                                            SORG = sorgr_lab, SGCON = sgcrit_lab, 
                                            SGCRIT = sgcrit_lab, KRGCL = 1, 
                                            KROGCG = 1, NG = sol[1], NOG = sol[2], 
                                            NP = 101))

colnames(rel_perm) <- c("Sg", "Sl", "Krg", "Krog")

p_forecast <- p

Gi_forecast <- Gi

wf_forecast <- wf

time_lst_forecast <- mbal_time(c(1:length(p_forecast)), "year")


forecast_lst <- mbal_forecast_param_oil(input_unit = "Field", output_unit = "Field",
                                        N = 1e6, m = 0.0, phi = 0.1, swi = 0.2,
                                        Gi = Gi_forecast, pb = 1688, p = p_forecast, pvt = pvt_table,
                                        cf = 2e-6, wf = wf_forecast, sorg = 0.15, rel_perm = rel_perm)

forecast_results <- mbal_forecast_oil(forecast_lst, time_lst_forecast)

forecast_results

reservoir_performance_table <- data.frame(p = p)

reservoir_performance_table$`RF_oil` <- c(0, 0.4,	0.5, 2.7,	4.4, 11.3, 16.1,
                                          19.3,	22.2,	24.3,	26.2,	27.9) / 100

reservoir_performance_table$`Sg` <- c(0, 0, 0, 2.9,	5.3, 14.8, 22.3, 27.3, 32.2,
                                      36.2,	39.9,	43.9) * 0.8 / 100

reservoir_performance_table$`GOR` <- c(0.84,	0.84,	0.84,	0.82,	0.8,	1.41,	2.17,
                                       2.7,	3.52,	4.58,	5.56,	6.79) * 1000 # SCF/STB

p1 <- forecast_results %>% ggplot(aes(`P (psia)`, SGo, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = reservoir_performance_table, aes(`p`, Sg, color = "Field"))+
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("Gas Saturation Plot") +
  theme_bw()

p1

 
p2 <- forecast_results %>% ggplot(aes(`P (psia)`, `GOR (SCF/STB)`, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = reservoir_performance_table, aes(`p`, GOR, color = "Field"))+
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("GOR Plot") +
  theme_bw()

p2


p3 <- forecast_results %>% ggplot(aes(`P (psia)`, `RF_oil`, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = reservoir_performance_table, aes(`p`, RF_oil, color = "Field"))+
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("Oil Recovery Plot") +
  theme_bw()

p3


p4 <- forecast_results %>% ggplot(aes(`P (psia)`, `Liq_volume`, color = "Forecast")) +
  geom_point(size = 3) +
  scale_color_manual(name="Data", values=c("Forecast" = "red")) +
  ggtitle("CCE Liquid Volume Plot") +
  theme_bw()

p4


p5 <- forecast_results  %>%
  tidyr::pivot_longer(cols = c("Isd", "Ifwd"), names_to = "Drive Mechanism",
                      values_to = "Fraction", values_drop_na = TRUE) %>%
  ggplot(aes(`P (psia)`, Fraction, fill = `Drive Mechanism`)) +
  geom_area() +
  ggtitle("Energy Plot") +
  theme_bw()

p5


