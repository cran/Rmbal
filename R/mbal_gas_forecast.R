
#' A list object of class 'forecast_gas' for material balance analysis
#'
#' Create an object of class 'forecast_gas'
#'
#' @param input_unit a unit system for parameters, only the character string 'Field' is accepted
#' @param output_unit a unit system for properties, only the character string 'Field' is accepted
#' @param G original gas in place, SCF.
#' @param phi reservoir porosity, a numeric fraction
#' @param swi initial water saturation in the reservoir, a numeric fraction
#' @param pd dew point pressure, a numeric value, psi
#' @param p reservoir pressure, a numeric vector, psi
#' @param pvt a data frame of PVT properties including pressure 'p' in 'psi', oil formation volume factor 'Bo' in 'bbl/stb', solution gas-oil ratio 'Rs' in 'scf/stb', oil viscosity 'muo' in 'cp', volatilized oil-gas ratio 'Rv' in 'stb/scf', gas formation volume factor 'Bg' in 'bbl/scf', gas viscosity 'mug' in 'cp', water formation volume factor 'Bw' in 'bbl/stb', and water viscosity 'muw' in 'cp'
#' @param cf formation compressibility, a numeric value or vector, 1/psi
#' @param M ratio of non-net-pay pore volume to the reservoir (net-pay) volume, a numeric fraction.
#' @param wf weight factor, a numeric vector of zeros and ones. A zero value excludes the entire row of reservoir history data at a particular time from the material balance analysis
#' @param rel_perm a data frame with four columns: gas saturation 'Sg', liquid saturation 'Sl', gas relative permeability 'Krg', and oil relative permeability 'Krog'
#'
#' @return a list of class ’forecast_gas’ with all the required parameters for the mbal_forecast_gas() S3 methods
#' @export
#'
#' @examples
#' p_pvt <- c(3700, 3650, 3400, 3100, 2800, 2500, 2200, 1900, 1600, 1300, 1000,
#' 700,  600, 400)

#' Bo <- c(10.057, 2.417, 2.192, 1.916, 1.736, 1.617, 1.504, 1.416, 1.326, 1.268,
#' 1.205, 1.149, 1.131, 1.093)
#'
#' Rv <- c(84.11765, 84.11765, 70.5, 56.2, 46.5, 39.5, 33.8, 29.9, 27.3, 25.5, 25.9,
#' 28.3, 29.8, 33.5) / 1e6
#'
#' Rs <- c(11566, 2378, 2010, 1569, 1272, 1067, 873, 719, 565, 461, 349, 249, 218,
#' 141)
#'
#' Bg <- c(0.87, 0.88, 0.92, 0.99, 1.08, 1.20, 1.35, 1.56, 1.85, 2.28, 2.95, 4.09,
#' 4.68, 6.53) / 1000
#'
#' cw <- 3e-6
#'
#' Bwi <- 10.05
#'
#' Bw <- Bwi * exp(cw * (p_pvt[1] - p_pvt))
#'
#' muo <- c(0.0612, 0.062, 0.1338, 0.1826, 0.2354, 0.3001, 0.3764, 0.4781, 0.6041,
#' 0.7746, 1.0295, 1.358, 1.855, 2.500)
#'
#' mug <- c(0.0612, 0.062, 0.0554, 0.0436, 0.0368, 0.0308, 0.0261, 0.0222, 0.0191,
#' 0.0166, 0.0148, 0.0135, 0.0125, 0.0115)
#'
#' muw <- rep(0.25, length(p_pvt))
#'
#' pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg, Bw = Bw,
#'  muo = muo, mug = mug, muw = muw)
#'
#' rel_perm <- as.data.frame(Rrelperm::kr2p_gl(SWCON = 0.2, SOIRG = 0.15,
#' SORG = 0.15, SGCON = 0.05, SGCRIT = 0.05, KRGCL = 1, KROGCG = 1,
#' NG = 3.16, NOG = 2.74, NP = 101))
#'
#' colnames(rel_perm) <- c("Sg", "Sl", "Krg", "Krog")
#'
#' p <- c(3700, 3650, 3400, 3100, 2800, 2500, 2200, 1900, 1600, 1300, 1000, 700,
#' 600)
#'
#' wf <- rep(1, length.out = length(p))
#'
#' forecast_lst <- mbal_forecast_param_gas(input_unit = "Field",
#' output_unit = "Field", G = 2.41e10, phi = 0.1, swi = 0.2, pd = 3650,
#' p = p, pvt = pvt_table, M = 0, cf = 2e-6, wf = wf,
#' rel_perm = rel_perm)
#'
#' dplyr::glimpse(forecast_lst)

mbal_forecast_param_gas <- function(input_unit = "Field", output_unit = "Field", G = NULL, phi = NULL, swi = NULL, pd = NULL, p = NULL, pvt = NULL, cf = NULL, M = NULL, wf = NULL, rel_perm = NULL) {

   if (!is.character(input_unit)) stop("'input_unit' must be the character string 'Field'.")
   if (input_unit != "Field") stop("'input_unit' must be the character string 'Field'.")
   if (!is.character(output_unit)) stop("'output_unit' must be the character string 'Field'.")
   if (output_unit != "Field") stop("'output_unit' must be the character string 'Field'.")
   if (is.null(G)) stop("'G' must be a numeric value.")
   if (!is.numeric(G)) stop("'G' must be a numeric value.")
   if (length(G) != 1) stop("'G' must be a numeric value.")
   if (is.null(phi)) stop("'phi' must be a numeric value.")
   if (!is.numeric(phi)) stop("'phi' must be a numeric value.")
   if (length(phi) != 1) stop("'phi' must be a numeric value.")
   if (is.numeric(phi)) {
      if (phi >= 1 | phi <= 0) {
         stop("Reservoir porosity must be greater than zero and less than one.")
      }
   }
   if (is.null(swi)) stop("'swi' must be a numeric value.")
   if (!is.numeric(swi)) stop("'swi' must be a numeric value.")
   if (length(swi) != 1) stop("'swi' must be a numeric value.")
   if (is.numeric(swi)) {
      if (swi >= 1 | swi <= 0) {
         stop("Reservoir initial water saturation must be greater than zero and less than one.")
      }
   }
   if (is.null(pd)) stop("'pd' must be a numeric value.")
   if (!is.numeric(pd)) stop("'pd' must be a numeric value.")
   if (length(pd) != 1) stop("'pd' must be a numeric value.")
   if (is.null(p)) stop("'p' must be a numeric vector.")
   if (!is.numeric(p)) stop("'p' must be a numeric vector.")
   l <- length(p)
   if (is.null(pvt)) stop("'pvt' must be a data frame with columns 'p', 'Bg', 'Rv', 'mug', 'Bo', 'Rs', 'muo', 'Bw', and 'muw'.")
   if (!is.data.frame(pvt)) stop("'pvt' must be a data frame with columns 'p', 'Bg', 'Rv', 'mug', 'Bo', 'Rs', 'muo', 'Bw', and 'muw'.")
   if (is.data.frame(pvt)) {
      if (nrow(pvt) < l) stop("Number of rows in the 'pvt' data frame must be equal or greater than the length of 'p'.")
      if (!('p' %in% colnames(pvt))) {
         stop("Column 'p' is missing in the 'pvt' data frame.")
      }
      if (!('Bo' %in% colnames(pvt))) {
         stop("Column 'Bo' is missing in the 'pvt' data frame.")
      }
      if (!('Rs' %in% colnames(pvt))) {
         stop("Column 'Rs' is missing in the 'pvt' data frame.")
      }
      if (!('muo' %in% colnames(pvt))) {
         stop("Column 'muo' is missing in the 'pvt' data frame.")
      }
      if (!('Rv' %in% colnames(pvt))) {
         stop("Column 'Rv' is missing in the 'pvt' data frame.")
      }
      if (any(is.na(pvt$Rv))) {
         stop("Column 'Rv' in the 'pvt' data frame does not accept 'NA' values.")
      }
      if (!('Bg' %in% colnames(pvt))) {
         stop("Column 'Bg' is missing in the 'pvt' data frame.")
      }
      if (any(pvt$Bg == 0, na.rm = TRUE)) {
         stop("Column 'Bg' in the 'pvt' data frame does not accept zero values.")
      }
      if (any(is.na(pvt$Bg))) {
         stop("Column 'Bg' in the 'pvt' data frame does not accept 'NA' values.")
      }
      if (!('mug' %in% colnames(pvt))) {
         stop("Column 'mug' is missing in the 'pvt' data frame.")
      }
      if (any(pvt$mug == 0, na.rm = TRUE)) {
         stop("Column 'mug' in the 'pvt' data frame does not accept zero values.")
      }
      if (any(is.na(pvt$mug))) {
         stop("Column 'mug' in the 'pvt' data frame does not accept 'NA' values.")
      }
      if (!('Bw' %in% colnames(pvt))) {
         stop("Column 'Bw' is missing in the 'pvt' data frame.")
      }
      if (max(p) > max(pvt$p)) {
         stop("Pressure range in the 'pvt' data frame does not cover the entire range of 'p' vector.")
      }
      if (min(p) < min(pvt$p)) {
         stop("Pressure range in the 'pvt' data frame does not cover the entire range of 'p' vector.")
      }
   }
   # if (is.null(cw)) stop("'cw' must be a numeric value or a numeric vector.")
   # if (!is.numeric(cw)) stop("'cw' must be a numeric value or a numeric vector.")
   # if (is.numeric(cw)) {
   #    if (!(length(cw) %in% c(1, l))) stop("'cw' must be a constant numeric value or a numeric vector with the same length as 'p'.")
   #    if (length(cw) == 1) {
   #       cw <- rep(cw,l)
   #    }
   # }
   if (is.null(cf)) stop("'cf' must be a numeric value or a numeric vector.")
   if (!is.numeric(cf)) stop("'cf' must be a numeric value or a numeric vector.")
   if (is.numeric(cf)) {
      if (!(length(cf) %in% c(1, l))) stop("'cf' must be a constant numeric value or a numeric vector with the same length as 'p'.")
      if (length(cf) == 1) {
         cf <- rep(cf,l)
      }
   }
   if (is.null(M)) stop("'M' must be a numeric value.")
   if (!is.numeric(M)) stop("'M' must be a numeric value.")
   if (length(M) != 1) stop("'M' must be a numeric value.")
   if (is.numeric(M)) {
      if (M < 0) {
         stop("The ratio of non-net-pay porve volume to net-pay pore volume must be equal to or greater than zero.")
      }
   }
   if (is.null(wf)) {
      wf <- rep(1,l)
   } else {
      if (!is.numeric(wf)) stop("'Wf' must be a numeric vector.")
      if (length(wf) != l) stop("'p' and 'Wf' vectors must have the same length.")
      wf[wf != 0] <- 1
   }
   if (is.null(rel_perm)) stop("'rel_perm' data frame is missing.")
   if (!(is.data.frame(rel_perm))) stop("'rel_perm' must be a data frame.")
   if (!("Sg" %in% names(rel_perm))) stop("'rel_perm' data frame must have columns 'Sg', 'Sl', 'Krg', and 'Krog'.")
   if (!("Krg" %in% names(rel_perm))) stop("'rel_perm' data frame must have columns 'Sg', 'Sl', 'Krg', and 'Krog'.")
   if (!("Krog" %in% names(rel_perm))) stop("'rel_perm' data frame must have columns 'Sg', 'Sl', 'Krg', and 'Krog'.")
   if (min(rel_perm$Sg) != 0) stop("'Sg' must cover saturation range [0,1] in 'rel_perm' data frame.")
   if (max(rel_perm$Sg) != 1) stop("'Sg' must cover saturation range [0,1] in 'rel_perm' data frame.")
   if (min(rel_perm$Sl) != 0) stop("'Sl' must cover saturation range [0,1] in 'rel_perm' data frame.")
   if (max(rel_perm$Sl) != 1) stop("'Sl' must cover saturation range [0,1] in 'rel_perm' data frame.")
   if (min(rel_perm$Krg) < 0) stop("'Krg' must cover relative permeability range [0,1] in 'rel_perm' data frame.")
   if (max(rel_perm$Krg) > 1) stop("'Krg' must cover relative permeability range [0,1] in 'rel_perm' data frame.")
   if (min(rel_perm$Krog) < 0) stop("'Krog' must cover relative permeability range [0,1] in 'rel_perm' data frame.")
   if (max(rel_perm$Krog) > 1) stop("'Krog' must cover relative permeability range [0,1] in 'rel_perm' data frame.")
   pvt <- pvt
   forecast_lst <- list(input_unit = input_unit, output_unit = output_unit, G = G, phi = phi, swi = swi, pd = pd, p = p, cf = cf, M = M, pvt = pvt, wf = wf, rel_perm = rel_perm)
   class(forecast_lst) <- c("volumetric_forecast_gas", "forecast_gas")

   return(forecast_lst)
}



# ******************************************************************************

#' Generic function for performance forecasting of a gas reservoir
#'
#' Generate a data frame of reservoir production estimates, and fluids saturations and liquid dropout in the gas leg according to the class of 'forecast_lst' and 'time_lst' objects
#'
#' @param forecast_lst a list object of class 'forecast_gas'
#' @param time_lst a list object of class 'time/date'
#'
#' @return a data frame with estimates for saturation of fluids, liquid dropout, gas-oil ratio, recovery factor, and drive indices over a range of given pressures
#'
#' @export
#'
#' @references
#' \insertRef{Walsh2003}{Rmbal}
#'
#' \insertRef{Walsh1994a}{Rmbal}
#'
#' \insertRef{Walsh1994}{Rmbal}
#'
#' \insertRef{Walsh1995}{Rmbal}
#'
#' \insertRef{Fetkovich1998}{Rmbal}
#'
#' @examples
#' p_pvt <- c(3700, 3650, 3400, 3100, 2800, 2500, 2200, 1900, 1600, 1300, 1000,
#' 700,  600, 400)
#'
#' Bo <- c(10.057, 2.417, 2.192, 1.916, 1.736, 1.617, 1.504, 1.416, 1.326, 1.268,
#' 1.205, 1.149, 1.131, 1.093)
#'
#' Rv <- c(84.11765, 84.11765, 70.5, 56.2, 46.5, 39.5, 33.8, 29.9, 27.3, 25.5, 25.9,
#' 28.3, 29.8, 33.5) / 1e6
#'
#' Rs <- c(11566, 2378, 2010, 1569, 1272, 1067, 873, 719, 565, 461, 349, 249, 218,
#' 141)
#'
#' Bg <- c(0.87, 0.88, 0.92, 0.99, 1.08, 1.20, 1.35, 1.56, 1.85, 2.28, 2.95, 4.09,
#' 4.68, 6.53) / 1000
#'
#' cw <- 3e-6
#'
#' Bwi <- 10.05
#'
#' Bw <- Bwi * exp(cw * (p_pvt[1] - p_pvt))
#'
#' muo <- c(0.0612, 0.062, 0.1338, 0.1826, 0.2354, 0.3001, 0.3764, 0.4781, 0.6041,
#' 0.7746, 1.0295, 1.358, 1.855, 2.500)
#'
#' mug <- c(0.0612, 0.062, 0.0554, 0.0436, 0.0368, 0.0308, 0.0261, 0.0222, 0.0191,
#' 0.0166, 0.0148, 0.0135, 0.0125, 0.0115)
#'
#' muw <- rep(0.25, length(p_pvt))
#'
#' pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg, Bw = Bw,
#'  muo = muo, mug = mug, muw = muw)
#'
#' rel_perm <- as.data.frame(Rrelperm::kr2p_gl(SWCON = 0.2, SOIRG = 0.15,
#' SORG = 0.15, SGCON = 0.05, SGCRIT = 0.05, KRGCL = 1, KROGCG = 1,
#' NG = 3.16, NOG = 2.74, NP = 101))
#'
#' colnames(rel_perm) <- c("Sg", "Sl", "Krg", "Krog")
#'
#' p <- c(3700, 3650, 3400, 3100, 2800, 2500, 2200, 1900, 1600, 1300, 1000, 700,
#' 600)
#'
#' Gi <- rep(0, length.out = length(p))
#'
#' wf <- rep(1, length.out = length(p))
#'
#' forecast_lst <- mbal_forecast_param_gas(input_unit = "Field",
#' output_unit = "Field", G = 2.41e10, phi = 0.1, swi = 0.2, pd = 3650,
#' p = p, pvt = pvt_table, M = 0, cf = 2e-6, wf = wf,
#' rel_perm = rel_perm)
#'
#' time_lst <- mbal_time(c(1:length(p)), "year")
#'
#' mbal_forecast_results <- mbal_forecast_gas(forecast_lst, time_lst)
#'
#' dplyr::glimpse(mbal_forecast_results)

mbal_forecast_gas <- function(forecast_lst, time_lst) {

   if (inherits(forecast_lst, "forecast_gas") == TRUE & inherits(time_lst, "time")) {
      UseMethod("mbal_forecast_gas")
   } else {
      if (!inherits(forecast_lst, "forecast_gas")) {
         stop("A class of 'forecast_gas' must be assigned to the 'forecast_lst' parameter of the mbal_forecast_gas() function.")
      }
      if (!inherits(time_lst, "time")) {
         stop("A class of 'time' must be assigned to the 'time_lst' parameter of the mbal_forecast_gas() function.")
      }
   }
}


# ******************************************************************************

#' S3 method for class 'mbal_forecast_gas'
#'
#' Return a data frame with estimates for saturation of fluids, liquid dropout, gas-oil ratio, recovery factor, and drive indices over a range of given pressures for a volumetric reservoir
#'
#' @param forecast_lst a list object of class 'forecast_gas'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame with estimates for saturation of fluids, liquid dropout, gas-oil ratio, recovery factor, and drive indices over a range of given pressures for a volumetric reservoir
#'
#' @export
#'
mbal_forecast_gas.volumetric_forecast_gas <- function(forecast_lst, time_lst) {
   volumetric_forecast_gas(forecast_lst, time_lst)
}



# ******************************************************************************

volumetric_forecast_gas <- function(forecast_lst, time_lst) {

   month_to_day <- 30.41667
   year_to_day <- 365
   if (time_lst$unit == "day") {
      time <- time_lst$t
   }
   if (time_lst$unit == "month") {
      time <- time_lst$t * month_to_day
   }
   if (time_lst$unit == "year") {
      time <- time_lst$t * year_to_day
   }
   G <- forecast_lst$G
   phi <- forecast_lst$phi
   swi <- forecast_lst$swi
   pd <- forecast_lst$pd
   p <- forecast_lst$p
   l <- length(p)
   if (length(time_lst$t) != length(p)) stop("Lenght of 'time_lst$t' vector must be equal to the length of 'p' vector.")
   cf <- forecast_lst$cf
   M <- forecast_lst$M
   pvt <- forecast_lst$pvt
   wf <- forecast_lst$wf
   rel_perm <- forecast_lst$rel_perm
   keep <- which(wf == 1)
   p <- p[keep]
   time <- time[keep]
   if (length(cf) != 1) cf <- cf[keep]
   l <- sum(wf)
   Bo <- vector(length = l)
   Rs <- vector(length = l)
   muo <- vector(length = l)
   Rv <- vector(length = l)
   Bg <- vector(length = l)
   mug <- vector(length = l)
   Bw <- vector(length = l)
   muw <- vector(length = l)
   for (i in 1:l) {
      Bo[i] <- approx(pvt$p, pvt$Bo, xout = p[i])$y
      Rs[i] <- approx(pvt$p, pvt$Rs, xout = p[i])$y
      muo[i] <- approx(pvt$p, pvt$muo, xout = p[i])$y
      Rv[i] <- approx(pvt$p, pvt$Rv, xout = p[i])$y
      Bg[i] <- approx(pvt$p, pvt$Bg, xout = p[i])$y
      mug[i] <- approx(pvt$p, pvt$mug, xout = p[i])$y
      Bw[i] <- approx(pvt$p, pvt$Bw, xout = p[i])$y
      muw[i] <- approx(pvt$p, pvt$muw, xout = p[i])$y
   }
   dp <- vector(length = l)
   Bto <- vector(length = l)
   Btg <- vector(length = l)
   Eo <- vector(length = l)
   Eowf <- vector(length = l)
   Eg <- vector(length = l)
   Egwf <- vector(length = l)
   Ew <- vector(length = l)
   Ef <- vector(length = l)
   Et <- vector(length = l)
   F_ <- vector(length = l)
   RF_oil <- vector(length = l)
   RF_gas <- vector(length = l)
   Igd <- vector(length = l)
   Inwd <- vector(length = l)
   Ifwd <- vector(length = l)
   Iawd <- vector(length = l)
   Itot <- vector(length = l)
   sw <- vector(length = l)
   so <- vector(length = l)
   sg <- vector(length = l)
   sw_t <- vector(length = l)
   so_t <- vector(length = l)
   sg_t <- vector(length = l)
   Nfoi <- 0
   Gfgi <- G
   Nfgi <- Gfgi * Rv[1]
   N <- Nfgi
   PV <- (Gfgi * Bg[1]) / (1 - swi)
   BV <- PV / phi
   W <- PV * swi / Bw[1]
   dp <- p[1] - p
   Bto <- (Bo * (1 - Rs[1] * Rv) + Bg * (Rs[1] - Rs)) / (1 - Rs * Rv)
   Btg <- (Bg * (1 - Rs * Rv[1]) + Bo * (Rv[1] - Rv)) / (1 - Rs * Rv)
   Eo <- Bto - Bto[1]
   Eg <- Btg - Btg[1]
   Ew <- Bw - Bw[1]
   for (i in 1:l) {
      Ef[i] <- -1 * trapz(p[1:i], cf[1:i])
   }
   Egwf <- Eg + Bg[1] * ((swi + M) * Ew / Bw[1] + (1 + M) * Ef) / (1 - swi)
   Et <- Gfgi * Egwf
   sgi <- 1- swi
   soi <- 1- swi - sgi
   denom <- vector(length = l)
   phi_o <- vector(length = l)
   phi_g <- vector(length = l)
   denom <- Egwf
   phi_o <- ifelse(dplyr::row_number(denom) == 1, 0, Rv[1] * (Bo - Rs * Bg) / (1 - Rs * Rv) / denom)
   phi_g <- ifelse(dplyr::row_number(denom) == 1, 0, (Bg - Rv * Bo) / (1 - Rs * Rv) / denom)
   ogr <- vector(length = l)
   RF_oil <- vector(length = l)
   RF_gas <- vector(length = l)
   Np <- vector(length = l)
   Gp <- vector(length = l)
   volfrac_liq <- vector(length = l)
   ogr[1] <- Rv[1]
   RF_oil[1] <- 0
   RF_gas[1] <- 0
   Gp[1] <- RF_gas[1] * G
   Np[1] <- RF_oil[1] * N
   swi_t <- swi
   sgi_t <- 1 - swi_t
   soi_t <- 1 - swi_t - sgi_t
   sw_t[1] <- swi
   so_t[1] <- ifelse(p[1] >= pd, 0, (1 - RF_oil[1] - Bg[1] * Rv[1] / Bg[1] / Rv[1]) * (1 - sw_t[1])  / (Bg[1] / Bo[1] / Rv[1] - Bg[1] * Rv[1] / Bg[1] / Rv[1]))
   sg_t[1] <- 1 - sw_t[1] - so_t[1]
   sw[1] <- swi
   so[1] <- so_t[1]
   sg[1] <- 1 - sw[1] - so[1]
   volfrac_liq[1] <- 0.0
   if (all(Rv == 0) | all(Rv == Rv[1])) {
      for (i in 2:l) {
         dGp_G <- Bg[1] / Bg[i - 1] - Bg[1] / Bg[i]
         RF_gas[i] <- dGp_G + RF_gas[i - 1]
         Gp[i] <- RF_gas[i] * G
         dNp_N <- dGp_G
         RF_oil[i] <- dNp_N + RF_oil[i - 1]
         Np[i] <- RF_oil[i] * N
         sw_t[i] <- swi
         so_t[i] <- 0
         sg_t[i] <- 1 - sw_t[i] - so_t[i]
         sw[i] <- swi
         so[i] <- 0
         sg[i] <- 1 - sw[i] - so[i]
         ogr[i] <- ogr[1]
         volfrac_liq[i] <- 0.0
      }

   } else {
      for (i in 2:l) {
         if (p[i] >= pd) {
            ogr_g <- ogr[1]
         } else {
            ogr_g <- ogr[i - 1] / 1.2
         }
         error <- 1e3
         while (abs(error) > 1e-6) {
            ogr_m <- mean(c(ogr_g, ogr[i - 1]))
            dGp_G <- (1 - RF_oil[i - 1] * phi_o[i] - RF_gas[i - 1] * phi_g[i]) /
               (phi_g[i] + ogr_m  * phi_o[i] / Rv[1])
            RF_gas[i] <- dGp_G + RF_gas[i - 1]
            Gp[i] <- RF_gas[i] * G
            dNp_N <- ogr_m * dGp_G / Rv[1]
            RF_oil[i] <- dNp_N + RF_oil[i - 1]
            Np[i] <- RF_oil[i] * N
            sw_t[i] <- swi
            so_t[i] <- ifelse(p[i] >= pd, 0, (1 - RF_oil[i] - Bg[1] * Rv[i] / Bg[i] / Rv[1]) *
                                 (1 - sw_t[i]) / (Bg[1] / Bo[i] / Rv[1] - Bg[1] * Rv[i] / Bg[i] / Rv[1]))
            sg_t[i] <- 1 - sw_t[i] - so_t[i]
            sw[i] <- sw_t[i]
            so[i] <- so_t[i]
            sg[i] <- 1 - sw[i] - so[i]
            krg <- approx(x = rel_perm$Sg, y = rel_perm$Krg, xout = sg[i], rule = 2)$y
            kro <- approx(x = rel_perm$Sg, y = rel_perm$Krog, xout = sg[i], rule = 2)$y
            gor_g_n <- (kro * Rs[i] / muo[i] / Bo[i] + krg / mug[i] / Bg[i]) / (kro / muo[i] / Bo[i] + krg * Rv[i] / mug[i] / Bg[i])
            ogr_g_n <- 1 / gor_g_n
            error <- ogr_g_n - ogr_g
            ogr_g <- ogr_g_n
            ogr[i] <- ogr_g
         }
         volfrac_liq[i] <- 1 / (1 + Bg[i] * (1 / Rv[1] - Rs[i]) / (Bo[i] * (1 - 1 / Rv[1] * Rv[i])))
      }
   }
   Igd <- ifelse(dplyr::row_number(Igd) == 1, NA, Gfgi * Eg / Et)
   Inwd <- ifelse(dplyr::row_number(Inwd) == 1, NA, 0)
   Ifwd <- ifelse(dplyr::row_number(Ifwd) == 1, NA, (W * Ew + PV * Ef) / Et)
   Iawd <- ifelse(dplyr::row_number(Iawd) == 1, NA, 0)
   Itot <- Igd + Inwd + Ifwd + Iawd
   names <- c("P (psia)", "SOg", "SGg", "SWg", "SOT", "SGT", "SWT", "GOR (SCF/STB)", "RF_oil", "RF_gas", "Liq_volume", "Igd", "Inwd", "Ifwd", "Iawd", "Itot")
   results <- data.frame(p = p,  SO = so, SG = sg, SW = sw, SOT = so_t, SGT = sg_t, SWT = sw_t, gor = 1/ ogr, RF_oil = RF_oil, RF_gas = RF_gas, volfrac_liq = volfrac_liq, Igd = Igd, Inwd = Inwd, Ifwd = Ifwd, Iawd = Iawd, Itot = Itot)
   colnames(results) <- names
   return(results)
}

