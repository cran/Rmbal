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

p_pvt <- c(5070, 4998, 4798, 4698, 4658, 4598, 4398, 4198, 3998, 3798, 3598, 3398, 
           3198, 2998, 2798, 2598, 2398, 2198, 1998, 1798, 1598, 1398, 1198, 998, 798, 598)
 
Bo <- c(2.704, 2.713, 2.740, 2.754, 2.707, 2.631, 2.338, 2.204, 2.093, 1.991, 1.905, 
        1.828, 1.758, 1.686, 1.632, 1.580, 1.534, 1.49, 1.45, 1.413, 1.367, 1.333, 1.305,
        1.272, 1.239, 1.205)  # RB/STB

Rv <- c(343, 343, 343, 343, 116, 111, 106, 94, 84, 74, 66, 60, 54, 49, 44, 39, 36, 
        33, 30, 28, 26, 25, 24.1, 23.9, 24.4, 26.4) / 1e6  # STB/SCF

Rs <- c(2909, 2909, 2909, 2909, 2834, 2711, 2247, 2019, 1828, 1651, 1500, 1364, 1237, 
        1111, 1013, 918, 833, 752, 677, 608, 524, 461, 406, 344, 283, 212)  # SCF/STB

Bg <- c(9.27472e-04, 9.30559e-04, 9.39820e-04, 9.44622e-04, 0.83, 0.835, 0.853, 
        0.874, 0.901, 0.933, 0.97, 1.015, 1.066, 1.125, 1.196, 1.281, 1.38, 1.498, 
        1.642, 1.819, 2.035, 2.315, 2.689, 3.19, 3.911, 5.034) / 1000  # RB/SCF

cw <- 3e-6

Bwi <- 1.05

Bw <- Bwi * exp(cw * (p_pvt[1] - p_pvt))

muo <- c(742, 735, 716, 706, 718, 739, 847, 906, 968, 1028, 1104, 1177, 1242, 1325, 
         1409, 1501, 1598, 1697, 1817, 1940, 2064, 2223, 2438, 2629, 2882, 3193) / 10000

mug <- c(742, 735, 716, 706, 375, 367, 350, 327, 306, 288, 271, 255, 240, 227, 214, 
         203, 193, 184, 175, 168, 161, 155, 150, 146, 142, 138) / 10000

muw <- rep(0.25, length(p_pvt))

liq_vol <- c(1000, 1000, 1000, 1000, 967, 847, 747, 683, 630, 584, 544, 508, 471, 433, 
             402, 368, 336, 305, 271, 239, 209, 177, 146, 117, 89, 63) / 1000

pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg, 
                               Bw = Bw, muo = muo, mug = mug, muw = muw)

p <- c(5070, 4998, 4798, 4698, 4658, 4598, 4398, 4198, 3998, 3798, 3598, 3398, 
           3198, 2998, 2798, 2598, 2398, 2198, 1998, 1798, 1598, 1398, 1198, 998, 798)

We <- rep(0, length.out = length(p))

Np <- c(0, 36, 130, 184, 227, 302, 582, 808, 1022, 1227, 1388, 1528, 1646, 1764, 1861, 
        1947, 2022, 2097, 2151, 2215, 2269, 2323, 2366, 2420, 2463) * 1e6 / 1000

Rp <- c(0, 2909.000, 2909.000, 2909.000, 2894.273, 2887.417, 2847.079, 2928.218, 3072.407, 
        3286.064, 3578.530, 3918.848, 4290.401, 4661.565, 5053.735, 5470.981, 5889.713, 
        6294.230, 6760.576, 7172.009, 7640.811, 8059.836, 8499.155, 8900.413, 9329.679)

Wp <- rep(0, length.out = length(p))

Wi <- rep(0, length.out = length(p))

Gi <- rep(0, length.out = length(p))

wf <- rep(1, length.out = length(p))

mbal_optim_oil_lst <- mbal_optim_param_oil(input_unit = "Field", output_unit = "Field",  
                                          unknown_param = "N", aquifer_model = NULL, 
                                          m = 0, phi = 0.1, swi = 0.2, Np = Np, 
                                          Rp = Rp, Wp = Wp, Gi = Gi, Wi = Wi, We = We, 
                                          pb = 4698, p = p, pvt = pvt_table, cf = 2e-6, 
                                          wf = wf, sorg = 0.15, sorw = 0.0)

time_lst <- mbal_time(c(1:length(p)), "year")

# a number of plots will be automatically generated for quality check

optim_results <- mbal_optim_oil(mbal_optim_oil_lst, time_lst)

glimpse(optim_results)

## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = FALSE-----

mbal_results <- mbal_perform_oil(optim_results, time_lst)

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

Gi_forecast <- Gi

wf_forecast <- wf

time_lst_forecast <- mbal_time(c(1:length(p_forecast)), "year")

forecast_lst <- mbal_forecast_param_oil(input_unit = "Field", output_unit = "Field",
                                        N = 10179044, m = 0.0, phi = 0.1, swi = 0.2,
                                        Gi = Gi_forecast, pb = 4698, p = p_forecast, 
                                        pvt = pvt_table, cf = 2e-6, wf = wf_forecast, 
                                        sorg = 0.15, rel_perm = rel_perm)

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
  geom_point(data = mbal_results, aes(`P (psia)`, `GOR (SCF/STB)`, color = "Field")) +
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("GOR Plot") +
  theme_bw()

p2


p3 <- forecast_results %>% ggplot(aes(`P (psia)`, `RF_oil`, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = mbal_results, aes(`P (psia)`, `RF_oil`, color = "Field")) +
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("Oil Recovery Plot") +
  theme_bw()

p3


liq_vol_CCE <- data.frame(P = p_pvt, liq_vol = liq_vol)

p4 <- forecast_results %>% ggplot(aes(`P (psia)`, `Liq_volume`, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = liq_vol_CCE, aes(P, `liq_vol`, color = "Lab")) +
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Lab" = "black")) +
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

## ----fig.align="center", fig.height=4, fig.width=6, warning=FALSE-------------

library(Rmbal)
library(Rrelperm)
library(pracma)
library(ggplot2)
library(dplyr)
library(magrittr)

p_pvt <- c(4998,	4798,	4698, 4658,	4598,	4398,	4198,	3998,	3798,	3598,	3398,	3198,	
           2998,	2798,	2598, 2398,	2198,	1998,	1798,	1598,	1398,	1198,	998, 798, 
           598) # psia

Bo <- c(2.71261,	2.73953,	2.75371, 2.70727,	2.63143,	2.33771,	2.20391,	2.09309,	
        1.99116,	1.90524,	1.82832, 1.75726,	1.68592,	1.63232,	1.58028,	1.53414,	
        1.49008,	1.44996,	1.41304, 1.36658,	1.33283,	1.30465,	1.27171,	1.23937,	
        1.20516) # RB/STB

Rs <- c(2909,	2909,	2909, 2834,	2711,	2247,	2019,	1828,	1651,	1500,	1364,	1237,	
        1111, 1013,	918,	833, 752,	677,	608,	524,	461,	406,	344,	283, 
        212) #SCF/STB

Bg <- c(0.932,	0.942,	0.947, 0.83,	0.835,	0.853,	0.874,	0.901,	0.933,	0.97,
        1.015,	1.066,	1.125, 1.196,	1.281,	1.38,	1.498,	1.642,	1.819,	2.035, 
        2.315,	2.689,	3.19, 3.911,	5.034) / 1000 # RB/SCF

Rv <- c(343,	343,	343, 116,	111,	106,	94,	84,	74,	66,	60,	54,	49,	44,	39,	
        36,	33,	30,	28,	26, 25,	24.1,	23.9,	24.4,	26.4) / 1e6 # STB/SCF

cw <- 3e-6

Bwi <- 1.05

Bw <- Bwi * exp(cw * (p_pvt[1] - p_pvt))

muo <- c(735,	716, 706, 718, 739, 847, 906, 968, 1028, 1104, 1177, 1242, 1325, 
         1409, 1501, 1598, 1697, 1817, 1940, 2064, 2223, 2438, 2629, 2882, 
         3193) / 10000

mug <- c(735,	716, 706, 375, 367, 350, 327, 306, 288, 271, 255, 240, 227, 214, 
         203, 193, 184, 175, 168, 161, 155, 150, 146, 142, 138) / 10000

muw <- rep(0.25, length(p_pvt))

pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg, 
                               Bw = Bw, muo = muo, mug = mug, muw = muw)

p <- c(4998,	4798,	4698, 4658,	4598,	4398,	4198,	3998,	3798,	3598,	3398,	3198,	
      2998,	2798,	2598, 2398,	2198,	1998,	1798,	1598,	1398,	1198,	998, 798, 
      598) # psia

Np <- c(0,	1.05,	1.63,	2.08,	2.88,	5.87,	8.36,	10.64,	12.63,	14.25,	15.65,	16.88,
        17.98,	18.9,	19.71,	20.41,	21.03,	21.59,	22.09,	22.57,	22.98,	23.34,
        23.69,	24.03, 24.38) / 100

Rp <- c(0,	2.909,	2.909, 2.909,	2.87,	2.86,	2.97,	3.18,	3.5,	3.87,	4.27,	4.68,	5.12,
        5.56,	6.02,	6.47,	6.91,	7.36,	7.8,	8.26,	8.68,	9.08,	9.48,	9.86,	
        10.25) * 1000   # SCF/STB

Wp <- rep(0, length.out = length(p))

We <- c(0, 0.0016, 0.0037, 0.005, 0.0067, 0.0149, 0.0261,	0.03096,	0.0559,	0.0742,	
        0.0952, 0.1179,	0.1424,	0.1688,	0.1966,	0.2261,	0.257, 0.2892,	0.3224,	
        0.357, 0.3923, 0.4285,	0.4658,	0.5038,	0.5425)   # RB

Gi <- rep(0, length.out = length(p))

Wi <- rep(0, length.out = length(p))

wf <- rep(1, length.out = length(p))

mbal_optim_oil_lst <- mbal_optim_param_oil(input_unit = "Field", output_unit = "Field",
                                           unknown_param = "N", aquifer_model = NULL,
                                           m = 0, phi = 0.1, swi = 0.2, Np = Np, 
                                           Rp = Rp, Wp = Wp, Gi = Gi, Wi = Wi, 
                                           We = We, pb = 4698, p = p, pvt = pvt_table, 
                                           cf = 2e-6, wf = wf, sorg = 0.15, sorw = 0.15)

time_lst <- mbal_time(c(1:length(p)), "year")

 # a number of plots will be automatically generated for quality check
 
 optim_results <- mbal_optim_oil(mbal_optim_oil_lst, time_lst)
 
 glimpse(optim_results)

## ---- fig.width = 6, fig.height= 4, fig.align = "center", warning = FALSE-----

mbal_results <- mbal_perform_oil(optim_results, time_lst)

mbal_results

reservoir_performance_table <- data.frame(p = p)

reservoir_performance_table$`RF_oil` <- c(0, 1.05, 1.63, 2.08, 2.88, 5.87, 8.36, 
                                          10.64, 12.63, 14.25,	15.65, 16.88, 17.98, 
                                          18.9, 19.71, 20.41, 21.03,	21.59, 22.09,	
                                          22.57, 22.98, 23.34, 23.69, 24.03, 
                                          24.38) / 100

reservoir_performance_table$`Sg` <- c(0,	0,	0,	3.4,	8.6,	26.1,	32.8,	37.7,	41.7,	
                                      44.7,	47.4,	49.7,	52.1,	53.8,	55.4,	57.1,	58.6,
                                      60.1,	61.6,	63.6,	65.1,	66.6,	68.4,	70.3,	
                                      72.5) * 0.8 / 100

reservoir_performance_table$`GOR` <- c(2.91,	2.91,	2.91,	2.83,	2.75,	2.98,	3.49,	4.48,
                                       5.9,	7.62,	9.14,	10.87,	12.96,	15.23,	18.09,
                                       20.54,	22.88,	25.63,	28.03,	31.09,	33.28,	
                                       35.48,	37.12,	37.3,	35.28) * 1000 # SCF/STB

p1 <- mbal_results %>% ggplot(aes(`P (psia)`, SGo, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = reservoir_performance_table, aes(`p`, Sg, color = "Field"))+
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("Gas Saturation Plot") +
  theme_bw()

p1


p2 <- mbal_results %>% ggplot(aes(`P (psia)`, `GOR (SCF/STB)`, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = reservoir_performance_table, aes(`p`, GOR, color = "Field"))+
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("GOR Plot") +
  theme_bw()

p2


p3 <- mbal_results %>% ggplot(aes(`P (psia)`, `RF_oil`, color = "Forecast")) +
  geom_point(size = 3) +
  geom_point(data = reservoir_performance_table, aes(`p`, RF_oil, color = "Field"))+
  scale_color_manual(name="Data", values=c("Forecast" = "red", "Field" = "black")) +
  ggtitle("Oil Recovery Plot") +
  theme_bw()

p3


p4 <- mbal_results  %>%
  tidyr::pivot_longer(cols = c("Inwd", "Isd", "Ifwd"), names_to = "Drive Mechanism",
                      values_to = "Fraction", values_drop_na = TRUE) %>%
  ggplot(aes(`P (psia)`, Fraction, fill = `Drive Mechanism`)) +
  geom_area() +
  ggtitle("Energy Plot") +
  theme_bw()

p4

