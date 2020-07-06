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

p_pvt <- c(3700, 3650, 3400, 3100, 2800, 2500, 2200, 1900, 1600, 1300, 1000, 700, 
           600, 400)
 
Bo <- c(10.057, 2.417, 2.192, 1.916, 1.736, 1.617, 1.504, 1.416, 1.326, 1.268, 1.205, 
        1.149, 1.131, 1.093)

Rv <- c(84.11765, 84.11765, 70.5, 56.2, 46.5, 39.5, 33.8, 29.9, 27.3, 25.5, 25.9, 
        28.3, 29.8, 33.5) / 1e6

Rs <- c(11566, 2378, 2010, 1569, 1272, 1067, 873, 719, 565, 461, 349, 249, 218, 141)

Bg <- c(0.87, 0.88, 0.92, 0.99, 1.08, 1.20, 1.35, 1.56, 1.85, 2.28, 2.95, 4.09, 
        4.68, 6.53) / 1000

cw <- 3e-6

Bwi <- 10.05

Bw <- Bwi * exp(cw * (p_pvt[1] - p_pvt))

muo <- c(0.0612, 0.062, 0.1338, 0.1826, 0.2354, 0.3001, 0.3764, 0.4781, 0.6041, 
         0.7746, 1.0295, 1.358, 1.855, 2.500)

mug <- c(0.0612, 0.062, 0.0554, 0.0436, 0.0368, 0.0308, 0.0261, 0.0222, 0.0191, 
         0.0166, 0.0148, 0.0135, 0.0125, 0.0115)

muw <- rep(0.25, length(p_pvt))

pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg, 
                               Bw = Bw, muo = muo, mug = mug, muw = muw)

p <- c(3700, 3650, 3400, 3100, 2800, 2500, 2200, 1900, 1600, 1300, 1000, 700, 
           600)

We <- rep(0, length.out = length(p))

Np <- c(0, 28.6, 93, 231, 270, 379, 481, 517.2, 549, 580, 675, 755, 803) *1e3

Gp <- c(0, 0.34, 1.2, 3.3, 4.3, 6.6, 9.1, 10.5, 12, 12.8, 16.4, 19.1, 20.5) * 1e9

Wp <- rep(0, length.out = length(p))

Wi <- rep(0, length.out = length(p))

wf <- rep(1, length.out = length(p))

mbal_optim_gas_lst <- mbal_optim_param_gas(input_unit = "Field", output_unit = "Field",  
                                          unknown_param = "G", aquifer_model = NULL, 
                                          phi = 0.1, swi = 0.2, Np = Np, 
                                          Gp = Gp, Wp = Wp, Wi = Wi, We = We, 
                                          pd = 3650, p = p, pvt = pvt_table, M = 0, 
                                          cf = 2e-6, wf = wf, sgrw = 0.15)

time_lst <- mbal_time(c(1:length(p)), "year")

# a number of plots will be automatically generated for quality check

optim_results <- mbal_optim_gas(mbal_optim_gas_lst, time_lst)

glimpse(optim_results)

## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = FALSE-----

mbal_results <- mbal_perform_gas(optim_results, time_lst)

mbal_results

## ----fig.align="center", fig.height=4, fig.width=6, warning=FALSE-------------

# Step I: generating a set of pseudo relative permeability curves using 
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

wf_forecast <- wf

time_lst_forecast <- mbal_time(c(1:length(p_forecast)), "year")

forecast_lst <- mbal_forecast_param_gas(input_unit = "Field", output_unit = "Field",
                                        G = 2.41e10, phi = 0.1, swi = 0.2,
                                        pd = 3650, p = p_forecast, pvt = pvt_table,
                                        M = 0, cf = 2e-6, wf = wf_forecast,
                                        rel_perm = rel_perm)

glimpse(forecast_lst)

forecast_results <- mbal_forecast_gas(forecast_lst, time_lst_forecast)

forecast_results

p1 <- forecast_results %>% ggplot(aes(`P (psia)`, SOg, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = mbal_results, aes(`P (psia)`, SOg, color = "Field"))+
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("Condensate Saturation Plot") +
  theme_bw()

p1


p2 <- forecast_results %>% ggplot(aes(`P (psia)`, `GOR (SCF/STB)`, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = mbal_results, aes(`P (psia)`, `GOR (SCF/STB)`, color = "Field")) +
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("GOR Plot") +
  theme_bw()

p2


p3 <- forecast_results %>% ggplot(aes(`P (psia)`, `RF_oil`, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = mbal_results, aes(`P (psia)`, `RF_oil`, color = "Field")) +
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("Condensate Recovery Plot") +
  theme_bw()

p3


p4 <- forecast_results %>% ggplot(aes(`P (psia)`, `Liq_volume`, color = "Forecast")) +
  geom_point(size = 3) +
  scale_color_manual(name="Data", values=c("Forecast" = "red")) +
  ggtitle("CCE Liquid Volume Plot") +
  theme_bw()

p4


p5 <- forecast_results  %>%
  tidyr::pivot_longer(cols = c("Igd", "Ifwd"), names_to = "Drive Mechanism",
                      values_to = "Fraction", values_drop_na = TRUE) %>%
  ggplot(aes(`P (psia)`, Fraction, fill = `Drive Mechanism`)) +
  geom_area() +
  ggtitle("Energy Plot") +
  theme_bw()

p5

## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = FALSE-----
library(Rmbal)
library(Rrelperm)
library(pracma)
library(ggplot2)
library(dplyr)
library(magrittr)

p_pvt <- c(8000, 7500, 7280, 7250, 7000, 6500, 6000, 5500, 5000, 4500, 4000, 3500,	
           3000,	2500,	2000,	1500,	1000) # psia

Bo <- c(12.732,	13.044,	13.192,	1.054, 1.041,	1.018,	1.002,	0.983,	0.965, 0.947,
        0.93,	0.913,	0.896,	0.877,	0.858,	0.839,	0.819) # RB/STB

Rs <- c(22527,	22527,	22527,	860,	819,	754,	704,	648,	593,	541, 490, 440,
        389,	336,	280,	228,	169) #SCF/STB

Bg <- c(0.565,	0.579,	0.586,	0.587,	0.595,	0.613,	0.634,	0.661, 0.694, 0.737,
        0.795,	0.817,	0.997,	1.178,	1.466,	1.963, 2.912) / 1000 # RB/SCF

Rv <- c(44.4,	44.4,	44.4,	44.3,	43.9,	40.3,	36.5,	32.9,	29.2,	25.4,	21.4,	17.6, 
        13.7, 10.5,	7.9, 5.8,	4.4) / 1e6 # STB/SCF

cw <- 2e-6

Bwi <- 1.0

Bw <- Bwi * exp(cw * (p_pvt[1] - p_pvt))


muo <- c(0.049,	0.047,	0.046,	19.541,	20.965,	23.958,	26.338,	29.633,	33.319, 
         37.401, 42.161,	47.465,	53.765,	61.887,	72.143,	83.478,	99.049)  # cp

muw <- rep(0.25, length(p_pvt))

mug <- c(0.049,	0.047,	0.046,	0.046,	0.0449,	0.042,	0.0393,	0.0366,	0.0339,	
         0.0312,	0.0283,	0.0254,	0.0225,	0.0198,	0.0174,	0.0154,	0.014) # cp

pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg, 
                               Bw = Bw, muo = muo, mug = mug, muw = muw)

p <- c(8000,	7500,	7280,	7250,	7000,	6500,	6000,	5500,	5000,	4500,	4000,	3500,	
           3000,	2500,	2000,	1500,	1000) # psia

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

forecast_lst <- mbal_forecast_param_gas(input_unit = "Field", output_unit = "Field",
                                        G = 1e6, phi = 0.1, swi = 0.2,
                                        pd = 7255, p = p_forecast, pvt = pvt_table,
                                        M = 0, cf = 2e-6, wf = wf_forecast,
                                        rel_perm = rel_perm)

forecast_results <- mbal_forecast_gas(forecast_lst, time_lst_forecast)

forecast_results


reservoir_performance_table <- data.frame(p = p)

reservoir_performance_table$`RF_oil` <- c(0,	2.4,	3.5,	4.2,	5.8, 7.6,	10.4,
                                                13.8,	16.6,	20,	23.1,	25.9,	28.5,	30.8,
                                                32.7,	34.1,	35.2) / 100

reservoir_performance_table$`Sg` <- c(100,	100,	100,	99.99,	99.91,	99.29,	98.68,	
                                      98.15,	97.66,	97.21,	96.79,	96.45,	96.15,	
                                      95.99,	95.9,	95.9, 95.95) * 0.8 / 100


reservoir_performance_table$`GOR` <- c(22.53,	22.53, 22.53, 22.55, 22.78, 24.81,
                                       27.4, 30.4, 34.25,	39.37, 46.73, 56.82, 72.99, 
                                       95.24, 126.58, 172.41, 227.27) * 1000 # SCF/STB


p1 <- forecast_results %>% ggplot(aes(`P (psia)`, SGg, color = "Forecast")) +
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
  tidyr::pivot_longer(cols = c("Igd", "Ifwd"), names_to = "Drive Mechanism",
                      values_to = "Fraction", values_drop_na = TRUE) %>%
  ggplot(aes(`P (psia)`, Fraction, fill = `Drive Mechanism`)) +
  geom_area() +
  ggtitle("Energy Plot") +
  theme_bw()

p5



## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = FALSE-----
library(Rmbal)
library(Rrelperm)
library(pracma)
library(ggplot2)
library(dplyr)
library(magrittr)

p_pvt <- c(5800,	5550,	5450,	5420,	5300,	4800,	4300,	3800,	3300,	2800,	2300,	
           1800,	1300,	800) # psia

Bo <- c(4.382,	4.441,	4.468,	2.378,	2.366,	2.032,	1.828,	1.674,	1.554,	
        1.448,	1.36,	1.279,	1.2,	1.131) # RB/STB

Rs <- c(6042,	6042,	6042,	2795,	2750,	2128,	1730,	1422,	1177,	960,	776,	607,	
        443,	293) #SCF/STB

Bg <- c(0.725,	0.735,	0.739,	0.74,	0.743,	0.758,	0.794,	0.854,	0.947,	
        1.09,	1.313,	1.677,	2.316,	3.695) / 1000 # RB/SCF

Rv <- c(165.5,	165.5,	165.5,	164.2,	156.6,	114,	89,	65.2,	48.3,	35,	25,	
        19,	15,	13.5) / 1e6 # STB/SCF

cw <- 2e-6

Bwi <- 1.0

Bw <- Bwi * exp(cw * (p_pvt[1] - p_pvt))


muo <- c(0.0612,	0.062,	0.0587,	0.135,	0.1338,	0.1826,	0.2354,	0.3001,	0.3764,	
         0.4781,	0.6041,	0.7746,	1.0295,	1.358)  # cp

muw <- rep(0.25, length(p_pvt))

mug <- c(0.0612,	0.062,	0.0587,	0.0581,	0.0554,	0.0436,	0.0368,	0.0308,	0.0261,	
         0.0222,	0.0191,	0.0166,	0.0148,	0.0135) # cp

pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg, 
                               Bw = Bw, muo = muo, mug = mug, muw = muw)

p <- c(5800,	5550,	5450,	5420,	5300,	4800,	4300,	3800,	3300,	2800,	2300,	
           1800,	1300,	800) # psia

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

forecast_lst <- mbal_forecast_param_gas(input_unit = "Field", output_unit = "Field",
                                        G = 1e6, phi = 0.1, swi = 0.2,
                                        pd = 5430, p = p_forecast, pvt = pvt_table,
                                        M = 0, cf = 2e-6, wf = wf_forecast,
                                        rel_perm = rel_perm)

forecast_results <- mbal_forecast_gas(forecast_lst, time_lst_forecast)

forecast_results

reservoir_performance_table <- data.frame(p = p)

reservoir_performance_table$`RF_oil` <- c(0,	1.3,	1.9,	2.1,	2.6,	7,	10.1,	13.3,
                                          16.2,	18.4,	20.2,	21.6,	22.8,	23.7) / 100

reservoir_performance_table$`Sg` <- c(100,	100,	100,	99.2,	95,	81.8,	78.7,	76.7,	76.3,
                                      76.6,	77.2,	78.3,	79.5,	80.6) * 0.8 / 100


reservoir_performance_table$`GOR` <- c(6.04,	6.04,	6.04,	6.09,	6.39,	8.77,	11.22,	
                                       15.29,	20.62, 28.45, 39.82, 52.42, 66.48,	
                                       73.99) * 1000 # SCF/STB


p1 <- forecast_results %>% ggplot(aes(`P (psia)`, SGg, color = "Forecast")) +
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
  tidyr::pivot_longer(cols = c("Igd", "Ifwd"), names_to = "Drive Mechanism",
                      values_to = "Fraction", values_drop_na = TRUE) %>%
  ggplot(aes(`P (psia)`, Fraction, fill = `Drive Mechanism`)) +
  geom_area() +
  ggtitle("Energy Plot") +
  theme_bw()

p5



