
#' A list object of class 'forecast_oil' for material balance analysis
#'
#' Create an object of class 'forecast_oil'
#'
#' @param input_unit a unit system for parameters, only the character string 'Field' is accepted
#' @param output_unit a unit system for properties, only the character string 'Field' is accepted
#' @param N original oil in place, STB
#' @param m ratio of original gas cap volume to original oil leg volume, a numeric fraction
#' @param phi reservoir porosity, a numeric fraction
#' @param swi initial water saturation in the reservoir, a numeric fraction
#' @param Gi cumulative gas injection, SCF
#' @param pb bubble point pressure, a numeric value, psi
#' @param p reservoir pressure, a numeric vector, psi
#' @param pvt a data frame of PVT properties including pressure 'p' in 'psi', oil formation volume factor 'Bo' in 'bbl/stb', solution gas-oil ratio 'Rs' in 'scf/stb', oil viscosity 'muo' in 'cp', volatilized oil-gas ratio 'Rv' in 'stb/scf', gas formation volume factor 'Bg' in 'bbl/scf', gas viscosity 'mug' in 'cp', water formation volume factor 'Bw' in 'bbl/stb', and water viscosity 'muw' in 'cp'
#' @param cf formation compressibility, a numeric value or vector, 1/psi
#' @param wf weight factor, a numeric vector of zeros and ones. A zero value excludes the entire row of reservoir history data at a particular time from the material balance analysis
#' @param sorg residual oil saturation in gas invaded zone (gas cap expansion or gas injection), a numeric fraction
#' @param rel_perm a data frame with four columns: gas saturation 'Sg', liquid saturation 'Sl', gas relative permeability 'Krg', and oil relative permeability 'Krog'
#'
#' @return a list of class ’forecast_oil’ with all the required parameters for the mbal_forecast_oil() S3 methods
#' @export
#'
#' @examples
#' p_pvt <- c(3330, 3150, 3000, 2850, 2700, 2550, 2400)
#'
#' Bo <- c(1.2511, 1.2353, 1.2222, 1.2122, 1.2022, 1.1922, 1.1822)
#'
#' Rs <- c(510, 477, 450, 425, 401, 375, 352)
#'
#' Bg <- c(0.00087, 0.00092, 0.00096, 0.00101, 0.00107, 0.00113, 0.00120)
#'
#' cw <- 2e-6
#'
#' Bwi <- 1.0
#'
#' Bw <- Bwi * exp(cw * (p_pvt[1] - p_pvt))
#'
#' Rv <- rep(0, length(p_pvt))
#'
#' muo <- rep(0.5, length(p_pvt))
#'
#' muw <- rep(0.25, length(p_pvt))
#'
#' mug <- rep(0.02, length(p_pvt))
#'
#' pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg,
#'                        Bw = Bw, muo = muo, mug = mug, muw = muw)
#'
#' rel_perm <- as.data.frame(Rrelperm::kr2p_gl(SWCON = 0.2, SOIRG = 0.10,
#' SORG = 0.10, SGCON = 0.05, SGCRIT = 0.05, KRGCL = 0.3, KROGCG = 1,
#' NG = 0.93, NOG = 10, NP = 101))
#'
#' colnames(rel_perm) <- c("Sg", "Sl", "Krg", "Krog")
#'
#' p <- c(3330, 3150, 3000, 2850, 2700, 2550, 2400)
#'
#' Gi <- rep(0, length.out = length(p))
#'
#' wf <- c(1, 1, 1, 0, 1, 0, 1)
#'
#' forecast_lst <- mbal_forecast_param_oil(input_unit = "Field",
#' output_unit = "Field", N = 1.37e8, m = 0.377, phi = 0.2, swi = 0.2, Gi = Gi,
#' pb = 3330, p = p, pvt = pvt_table, cf = 0, wf = wf, sorg = 0.2,
#' rel_perm = rel_perm)
#'
#' dplyr::glimpse(forecast_lst)

mbal_forecast_param_oil <- function(input_unit = "Field", output_unit = "Field", N = NULL, m = NULL, phi = NULL, swi = NULL, Gi = NULL, pb = NULL, p = NULL, pvt = NULL, cf = NULL, wf = NULL, sorg = NULL, rel_perm = NULL) {

   if (!is.character(input_unit)) stop("'input_unit' must be the character string 'Field'.")
   if (input_unit != "Field") stop("'input_unit' must be the character string 'Field'.")
   if (!is.character(output_unit)) stop("'output_unit' must be the character string 'Field'.")
   if (output_unit != "Field") stop("'output_unit' must be the character string 'Field'.")
   if (is.null(N)) stop("'N' must be a numeric value.")
   if (!is.numeric(N)) stop("'N' must be a numeric value.")
   if (length(N) != 1) stop("'N' must be a numeric value.")
   if (is.null(m)) stop("'m' must be a numeric value.")
   if (!is.null(m)) {
      if (!is.numeric(m)) stop("'m' must be a numeric value.")
      if (m < 0) stop("'m' must be a numeric value equal or more than zero.")
   }
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
   if (is.null(pb)) stop("'pb' must be a numeric value.")
   if (!is.numeric(pb)) stop("'pb' must be a numeric value.")
   if (length(pb) != 1) stop("'pb' must be a numeric value.")
   if (is.null(p)) stop("'p' must be a numeric vector.")
   if (!is.numeric(p)) stop("'p' must be a numeric vector.")
   if (p[1] < pb) stop("Initial reservoir pressure must be equal to or greater than 'pb'.")
   if (m > 0) {
      if (max(p) != pb) {
         stop("Initial reservoir pressure must be equal to 'pb' in reservoirs with an associated gas cap.")
      }
   }
   l <- length(p)
   if (is.null(pvt)) stop("'pvt' must be a data frame with columns 'p', 'Bo', 'Rs', 'muo', 'Rv', 'Bg', 'mug', and 'Bw'.")
   if (!is.data.frame(pvt)) stop("'pvt' must be a data frame with columns 'p', 'Bo', 'Rs', 'muo', 'Rv', 'Bg', 'mug', and 'Bw'.")
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
      if (!('muw' %in% colnames(pvt))) {
         stop("Column 'muw' is missing in the 'pvt' data frame.")
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
   if (is.null(Gi)) {
      Gi <- rep(0, l)
   }
   if (!is.null(Gi)) {
      if (!is.numeric(Gi)) stop("'Gi' must be a numeric vector.")
      if (is.numeric(Gi)) {
         if (length(Gi) != l) stop("'p' and 'Gi' vectors must have the same length.")
      }
   }
   if (m > 0) {
      if (is.null(sorg)) stop("A numeric value must be assigned to 'sorg' for a 'gas_cap' drive reservoir (m > 0).")

   }
   if (m == 0) {
      if (is.null(sorg)) {
         if (any(Gi > 0)) stop("A numeric value must be assigned to 'sorg' for a gas injection case.")
      }
   }
   if (is.null(wf)) {
      wf <- rep(1,l)
   } else {
      if (!is.numeric(wf)) stop("'Wf' must be a numeric vector.")
      if (length(wf) != l) stop("'p' and 'Wf' vectors must have the same length.")
      wf[wf != 0] <- 1
   }
   if (m > 0) {
      if (wf[1] == 0) {
         stop("Information at 'pb' cannot be removed from the calculations. Change the corresponding 'wf' value to one.")
      }
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
   inj <- data.frame(Gi = Gi)
   forecast_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, inj = inj, wf = wf, sorg = sorg, rel_perm = rel_perm)
   if (m > 0) {
      class(forecast_lst) <- c("gas_cap_forecast_oil", "forecast_oil")
   }
   if (m == 0) {
      class(forecast_lst) <- c("volumetric_forecast_oil", "forecast_oil")
   }
   return(forecast_lst)
}



# ******************************************************************************

#' Generic function for performance forecasting of an oil reservoir
#'
#' Generate a data frame of reservoir production estimates, and fluids saturations and liquid dropout in the oil leg according to the class of 'forecast_lst' and 'time_lst' objects
#'
#' @param forecast_lst a list object of class 'forecast_oil'
#' @param time_lst a list object of class 'time/date'
#'
#' @return a data frame with estimates for saturation of fluids, liquid dropout, gas-oil ratio, recovery factor, and drive indices over a range of given pressures
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
#' @examples
#' p_pvt <- c(3330, 3150, 3000, 2850, 2700, 2550, 2400)
#'
#' Bo <- c(1.2511, 1.2353, 1.2222, 1.2122, 1.2022, 1.1922, 1.1822)
#'
#' Rs <- c(510, 477, 450, 425, 401, 375, 352)
#'
#' Bg <- c(0.00087, 0.00092, 0.00096, 0.00101, 0.00107, 0.00113, 0.00120)
#'
#' cw <- 2e-6
#'
#' Bwi <- 1.0
#'
#' Bw <- Bwi * exp(cw * (p_pvt[1] - p_pvt))
#'
#' Rv <- rep(0, length(p_pvt))
#'
#' muo <- rep(0.5, length(p_pvt))
#'
#' muw <- rep(0.25, length(p_pvt))
#'
#' mug <- rep(0.02, length(p_pvt))
#'
#' pvt_table <- data.frame(p = p_pvt, Bo = Bo, Rs = Rs, Rv = Rv, Bg = Bg,
#'                        Bw = Bw, muo = muo, mug = mug, muw = muw)
#'
#' rel_perm <- as.data.frame(Rrelperm::kr2p_gl(SWCON = 0.2, SOIRG = 0.10,
#' SORG = 0.10, SGCON = 0.05, SGCRIT = 0.05, KRGCL = 0.3, KROGCG = 1,
#' NG = 0.93, NOG = 10, NP = 101))
#'
#' colnames(rel_perm) <- c("Sg", "Sl", "Krg", "Krog")
#'
#' p <- c(3330, 3150, 3000, 2850, 2700, 2550, 2400)
#'
#' Gi <- rep(0, length.out = length(p))
#'
#' wf <- c(1, 1, 1, 0, 1, 0, 1)
#'
#' forecast_lst <- mbal_forecast_param_oil(input_unit = "Field",
#' output_unit = "Field", N = 1.37e8, m = 0.377, phi = 0.2, swi = 0.2, Gi = Gi,
#' pb = 3330, p = p, pvt = pvt_table, cf = 0, wf = wf, sorg = 0.2,
#' rel_perm = rel_perm)
#'
#' time_lst <- mbal_time(c(0, 365, 730, 1095, 1460, 1825, 2190), "day")
#'
#' mbal_forecast_results <- mbal_forecast_oil(forecast_lst, time_lst)
#'
#' dplyr::glimpse(mbal_forecast_results)

mbal_forecast_oil <- function(forecast_lst, time_lst) {

   if (inherits(forecast_lst, "forecast_oil") == TRUE & inherits(time_lst, "time")) {
      UseMethod("mbal_forecast_oil")
   } else {
      if (!inherits(forecast_lst, "forecast_oil")) {
         stop("A class of 'forecast_oil' must be assigned to the 'forecast_lst' parameter of the mbal_forecast_oil() function.")
      }
      if (!inherits(time_lst, "time")) {
         stop("A class of 'time' must be assigned to the 'time_lst' parameter of the mbal_forecast_oil() function.")
      }
   }
}


# ******************************************************************************

#' S3 method for class 'mbal_forecast_oil'
#'
#' Return a data frame with estimates for saturation of fluids, liquid dropout, gas-oil ratio, recovery factor, and drive indices over a range of given pressures for a volumetric oil reservoir
#'
#' @param forecast_lst a list object of class 'forecast_oil'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame with estimates for saturation of fluids, liquid dropout, gas-oil ratio, recovery factor, and drive indices over a range of given pressures for a volumetric oil reservoir
#' @export
mbal_forecast_oil.volumetric_forecast_oil <- function(forecast_lst, time_lst) {
   volumetric_forecast_oil(forecast_lst, time_lst)
}


# ******************************************************************************

#' S3 method for class 'mbal_forecast_oil'
#'
#' Return a data frame with estimates for saturation of fluids, liquid dropout, gas-oil ratio, recovery factor, and drive indices over a range of given pressures for a gas_cap_drive oil reservoir
#'
#' @param forecast_lst a list object of class 'forecast_oil'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame with estimates for saturation of fluids, liquid dropout, gas-oil ratio, recovery factor, and drive indices over a range of given pressures for a gas_cap_drive oil reservoir
#' @export
mbal_forecast_oil.gas_cap_forecast_oil <- function(forecast_lst, time_lst) {
   gas_cap_forecast_oil(forecast_lst, time_lst)
}


# ******************************************************************************


volumetric_forecast_oil <- function(forecast_lst, time_lst) {

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
   N <- forecast_lst$N
   m <- forecast_lst$m
   phi <- forecast_lst$phi
   swi <- forecast_lst$swi
   pb <- forecast_lst$pb
   p <- forecast_lst$p
   l <- length(p)
   if (length(time_lst$t) != length(p)) stop("Lenght of 'time_lst$t' vector must be equal to the length of 'p' vector.")
   cf <- forecast_lst$cf
   pvt <- forecast_lst$pvt
   inj <- forecast_lst$inj
   wf <- forecast_lst$wf
   sorg <- forecast_lst$sorg
   if (is.null(sorg)) stop("'sorg' must be a numeric value.")
   rel_perm <- forecast_lst$rel_perm
   keep <- which(wf == 1)
   p <- p[keep]
   time <- time[keep]
   if (length(cf) != 1) cf <- cf[keep]
   inj <- inj[keep,]
   l <- sum(wf)
   Bo <- vector(length = l)
   Rs <- vector(length = l)
   muo <- vector(length = l)
   Rv <- vector(length = l)
   Bg <- vector(length = l)
   mug <- vector(length = l)
   Bw <- vector(length = l)
   for (i in 1:l) {
      Bo[i] <- approx(pvt$p, pvt$Bo, xout = p[i])$y
      Rs[i] <- approx(pvt$p, pvt$Rs, xout = p[i])$y
      muo[i] <- approx(pvt$p, pvt$muo, xout = p[i])$y
      Rv[i] <- approx(pvt$p, pvt$Rv, xout = p[i])$y
      Bg[i] <- approx(pvt$p, pvt$Bg, xout = p[i])$y
      mug[i] <- approx(pvt$p, pvt$mug, xout = p[i])$y
      Bw[i] <- approx(pvt$p, pvt$Bw, xout = p[i])$y
   }
   Gi <- inj
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
   Isd <- vector(length = l)
   Inwd <- vector(length = l)
   Ifwd <- vector(length = l)
   Iawd <- vector(length = l)
   Itot <- vector(length = l)
   PVgas_r <- vector(length = l)
   sw <- vector(length = l)
   so <- vector(length = l)
   sg <- vector(length = l)
   sw_t <- vector(length = l)
   so_t <- vector(length = l)
   sg_t <- vector(length = l)
   Nfoi <- N / (1 + (m * Bo[1] * Rv[1]) / Bg[1])
   Gfgi <- m * Nfoi * Bo[1] / Bg[1]
   G <- Gfgi + Nfoi * Rs[1]
   PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
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
   Eowf <- Eo + Bo[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
   Egwf <- Eg + Bg[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
   Et <- Gfgi * Eg + Nfoi * Eo + W * Ew + PV * Ef
   soi <- 1- swi
   sgi <- 1- swi - soi
   PVgasi <- (Gfgi * Bg[1] * soi - Nfoi * Bo[1] * sgi) / ((1 - sorg - swi) * soi - sorg * sgi)
   PVoili <- PV - PVgasi
   xp <- 1
   PVgas_r <- PVgasi / PV + PVgasi / PV * xp * ((sorg * Rs[1] / Bo[1] + (1 - sorg - swi) / Bg[1]) / (sorg * Rs[1] / Bo + (1 - sorg - swi) / Bg)- 1) + Gi / ((sorg * Rs / Bo) + ((1 - sorg - swi) / Bg)) / PV
   denom <- vector(length = l)
   phi_o <- vector(length = l)
   phi_g <- vector(length = l)
   denom <- (Bto - Bo[1]) + (m * Bo[1] / Bg[1]) * (Btg - Bg[1]) + Gi / G  * (m * Bo[1] / Bg[1] + Rs[1]) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + (W * Ew + PV * Ef) * (1 + m * Bo[1] * Rv[1] / Bg[1]) / N
   phi_o <- ifelse(dplyr::row_number(denom) == 1, 0, ((1 + m * Bo[1] * Rv[1] / Bg[1]) * ((Bo - Bg * Rs) / (1 - Rv * Rs))) / denom)
   phi_g <- ifelse(dplyr::row_number(denom) == 1, 0, ((m * Bo[1] / Bg[1] + Rs[1]) * ((Bg - Bo * Rv) / (1 - Rv * Rs))) / denom)
   gor <- vector(length = l)
   RF_oil <- vector(length = l)
   RF_gas <- vector(length = l)
   Np <- vector(length = l)
   Gp <- vector(length = l)
   volfrac_liq <- vector(length = l)
   gor[1] <- Rs[1]
   RF_oil[1] <- 0
   RF_gas[1] <- 0
   Gp[1] <- RF_gas[1] * G
   swi_t <- swi
   soi_t <- soi
   sgi_t <- 1 - swi_t - soi_t
   sw[1] <- swi
   so[1] <- 1 - sw[1]
   sg[1] <- 1 - sw[1] - so[1]
   sw_t[1] <- sw[1]
   so_t[1] <- so[1] * (1 - PVgas_r[1]) + sorg * PVgas_r[1]
   sg_t[1] <- 1 - so_t[1] - sw_t[1]
   volfrac_liq[1] <- 1.0
   for (i in 2:l) {
      if (p[i] >= pb) {
         gor_g <- gor[1]
      } else {
         gor_g <- 1.1 * gor[i - 1]
      }
      error <- 1e3
      for (k in 1:100) {
         gor_m <- mean(c(gor_g, gor[i - 1]))
         dNp_N <- (1 - RF_oil[i - 1] * phi_o[i] - RF_gas[i - 1] * phi_g[i]) /
            (phi_o[i] + gor_m * (1 + m * Bo[1] * Rv[1] / Bg[1]) * phi_g[i] / (m * Bo[1] / Bg[1] + Rs[1]))
         RF_oil[i] <- dNp_N + RF_oil[i - 1]
         Np[i] <- RF_oil[i] * N
         dGp_G <- gor_m * dNp_N * (1 + m * Bo[1] * Rv[1] / Bg[1]) / (m * Bo[1] / Bg[1] + Rs[1])
         RF_gas[i] <- dGp_G + RF_gas[i - 1]
         Gp[i] <- RF_gas[i] * G
         if (p[i] >= pb) {
            sw[i] <- swi
            so[i] <- 1 - sw[i]
            sg[i] <- 1 - sw[i] - so[i]
            sw_t[i] <- sw[i]
            so_t[i] <- so[i] * (1 - PVgas_r[i]) + sorg * PVgas_r[i]
            sg_t[i] <- 1 - so_t[i] - sw_t[i]
         } else {
            sw_t[i] <- swi
            so_t[i] <- ((1 - Np[i] / N) * (soi_t / Bo[1] + sgi_t * Rv[1] / Bg[1]) - (1 - sw_t[i]) * Rv[i] / Bg[i]) / (1 / Bo[i] - Rv[i] / Bg[i])
            sg_t[i] <- 1 - so_t[i] - sw_t[i]
            sw[i] <- swi
            so[i] <- (so_t[i] - sorg * PVgas_r[i]) / (1 - PVgas_r[i])
            sg[i] <- 1 - sw[i] - so[i]
         }
         krg <- approx(x = rel_perm$Sg, y = rel_perm$Krg, xout = sg[i], rule = 2)$y
         kro <- approx(x = rel_perm$Sg, y = rel_perm$Krog, xout = sg[i], rule = 2)$y
         gor_g_n <- (kro * Rs[i] / muo[i] / Bo[i] + krg / mug[i] / Bg[i]) / (kro / muo[i] / Bo[i] + krg * Rv[i] / mug[i] / Bg[i])
         error <- gor_g_n - gor_g
         gor_g <- gor_g_n
         gor[i] <- gor_g
         if (abs(error) <= 1e-6) {
            break
         }
         if (k == 100) {
            RF_oil[i] <- NA
            RF_gas[i] <- NA
            Gp[i] <- NA
            sw_t[i] <- NA
            so_t[i] <- NA
            sg_t[i] <- NA
            sw[i] <- NA
            so[i] <- NA
            sg[i] <- NA
         }
      }
      # while (abs(error) > 1e-6) {
      #    gor_m <- mean(c(gor_g, gor[i - 1]))
      #    dNp_N <- (1 - RF_oil[i - 1] * phi_o[i] - RF_gas[i - 1] * phi_g[i]) /
      #       (phi_o[i] + gor_m * (1 + m * Bo[1] * Rv[1] / Bg[1]) * phi_g[i] / (m * Bo[1] / Bg[1] + Rs[1]))
      #    RF_oil[i] <- dNp_N + RF_oil[i - 1]
      #    Np[i] <- RF_oil[i] * N
      #    dGp_G <- gor_m * dNp_N * (1 + m * Bo[1] * Rv[1] / Bg[1]) / (m * Bo[1] / Bg[1] + Rs[1])
      #    RF_gas[i] <- dGp_G + RF_gas[i - 1]
      #    Gp[i] <- RF_gas[i] * G
      #    if (p[i] >= pb) {
      #       sw[i] <- swi
      #       so[i] <- 1 - sw[i]
      #       sg[i] <- 1 - sw[i] - so[i]
      #       sw_t[i] <- sw[i]
      #       so_t[i] <- so[i] * (1 - PVgas_r[i]) + sorg * PVgas_r[i]
      #       sg_t[i] <- 1 - so_t[i] - sw_t[i]
      #    } else {
      #       sw_t[i] <- swi
      #       so_t[i] <- ((1 - Np[i] / N) * (soi_t / Bo[1] + sgi_t * Rv[1] / Bg[1]) - (1 - sw_t[i]) * Rv[i] / Bg[i]) / (1 / Bo[i] - Rv[i] / Bg[i])
      #       sg_t[i] <- 1 - so_t[i] - sw_t[i]
      #       sw[i] <- swi
      #       so[i] <- (so_t[i] - sorg * PVgas_r[i]) / (1 - PVgas_r[i])
      #       sg[i] <- 1 - sw[i] - so[i]
      #    }
      #    krg <- approx(x = rel_perm$Sg, y = rel_perm$Krg, xout = sg[i], rule = 2)$y
      #    kro <- approx(x = rel_perm$Sg, y = rel_perm$Krog, xout = sg[i], rule = 2)$y
      #    gor_g_n <- (kro * Rs[i] / muo[i] / Bo[i] + krg / mug[i] / Bg[i]) / (kro / muo[i] / Bo[i] + krg * Rv[i] / mug[i] / Bg[i])
      #    error <- gor_g_n - gor_g
      #    gor_g <- gor_g_n
      #    gor[i] <- gor_g
      # }
      volfrac_liq[i] <- 1 / (1 + Bg[i] * (Rs[1] - Rs[i]) / (Bo[i] * (1 - Rs[1] * Rv[i])))
   }
   Igd <- ifelse(dplyr::row_number(Igd) == 1, NA, Gfgi * Eg / Et)
   Isd <- ifelse(dplyr::row_number(Isd) == 1, NA, Nfoi * Eo / Et)
   Inwd <- ifelse(dplyr::row_number(Inwd) == 1, NA, 0)
   Ifwd <- ifelse(dplyr::row_number(Ifwd) == 1, NA, (W * Ew + PV * Ef) / Et)
   Iawd <- ifelse(dplyr::row_number(Iawd) == 1, NA, 0)
   Itot <- Igd + Isd + Inwd + Ifwd + Iawd
   names <- c("P (psia)", "SOo", "SGo", "SWo", "SOT", "SGT", "SWT", "GOR (SCF/STB)", "RF_oil", "RF_gas", "Liq_volume", "Igd", "Isd", "Inwd", "Ifwd", "Iawd", "Itot")
   results <- data.frame(p = p, SO = so, SG = sg, SW = sw, SOT = so_t, SGT = sg_t, SWT = sw_t, gor = gor, RF_oil = RF_oil, RF_gas = RF_gas, volfrac_liq = volfrac_liq, Igd = Igd, Isd = Isd, Inwd = Inwd, Ifwd = Ifwd, Iawd = Iawd, Itot = Itot)
   colnames(results) <- names
   return(results)
}



gas_cap_forecast_oil <- function(forecast_lst, time_lst) {

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
   N <- forecast_lst$N
   m <- forecast_lst$m
   phi <- forecast_lst$phi
   swi <- forecast_lst$swi
   pb <- forecast_lst$pb
   p <- forecast_lst$p
   l <- length(p)
   if (length(time_lst$t) != length(p)) stop("Lenght of 'time_lst$t' vector must be equal to the length of 'p' vector.")
   cf <- forecast_lst$cf
   pvt <- forecast_lst$pvt
   inj <- forecast_lst$inj
   wf <- forecast_lst$wf
   sorg <- forecast_lst$sorg
   if (is.null(sorg)) stop("'sorg' must be a numeric value.")
   rel_perm <- forecast_lst$rel_perm
   keep <- which(wf == 1)
   p <- p[keep]
   time <- time[keep]
   if (length(cf) != 1) cf <- cf[keep]
   inj <- inj[keep,]
   l <- sum(wf)
   Bo <- vector(length = l)
   Rs <- vector(length = l)
   muo <- vector(length = l)
   Rv <- vector(length = l)
   Bg <- vector(length = l)
   mug <- vector(length = l)
   Bw <- vector(length = l)
   for (i in 1:l) {
      Bo[i] <- approx(pvt$p, pvt$Bo, xout = p[i])$y
      Rs[i] <- approx(pvt$p, pvt$Rs, xout = p[i])$y
      muo[i] <- approx(pvt$p, pvt$muo, xout = p[i])$y
      Rv[i] <- approx(pvt$p, pvt$Rv, xout = p[i])$y
      Bg[i] <- approx(pvt$p, pvt$Bg, xout = p[i])$y
      mug[i] <- approx(pvt$p, pvt$mug, xout = p[i])$y
      Bw[i] <- approx(pvt$p, pvt$Bw, xout = p[i])$y
   }
   Gi <- inj
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
   Isd <- vector(length = l)
   Inwd <- vector(length = l)
   Ifwd <- vector(length = l)
   Iawd <- vector(length = l)
   Itot <- vector(length = l)
   PVgas_r <- vector(length = l)
   sw <- vector(length = l)
   so <- vector(length = l)
   sg <- vector(length = l)
   sw_t <- vector(length = l)
   so_t <- vector(length = l)
   sg_t <- vector(length = l)
   Nfoi <- N / (1 + (m * Bo[1] * Rv[1]) / Bg[1])
   Gfgi <- m * Nfoi * Bo[1] / Bg[1]
   G <- Gfgi + Nfoi * Rs[1]
   PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
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
   Eowf <- Eo + Bo[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
   Egwf <- Eg + Bg[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
   Et <- Gfgi * Eg + Nfoi * Eo + W * Ew + PV * Ef
   soi <- 1- swi
   sgi <- 1- swi - soi
   PVgasi <- (Gfgi * Bg[1] * soi - Nfoi * Bo[1] * sgi) / ((1 - sorg - swi) * soi - sorg * sgi)
   PVoili <- PV - PVgasi
   xp <- 1
   PVgas_r <- PVgasi / PV + PVgasi / PV * xp * ((sorg * Rs[1] / Bo[1] + (1 - sorg - swi) / Bg[1]) / (sorg * Rs[1] / Bo + (1 - sorg - swi) / Bg)- 1) + Gi / ((sorg * Rs / Bo) + ((1 - sorg - swi) / Bg)) / PV
   denom <- vector(length = l)
   phi_o <- vector(length = l)
   phi_g <- vector(length = l)
   denom <- (Bto - Bo[1]) + (m * Bo[1] / Bg[1]) * (Btg - Bg[1]) + Gi / G  * (m * Bo[1] / Bg[1] + Rs[1]) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + (W * Ew + PV * Ef) * (1 + m * Bo[1] * Rv[1] / Bg[1]) / N
   phi_o <- ifelse(dplyr::row_number(denom) == 1, 0, ((1 + m * Bo[1] * Rv[1] / Bg[1]) * ((Bo - Bg * Rs) / (1 - Rv * Rs))) / denom)
   phi_g <- ifelse(dplyr::row_number(denom) == 1, 0, ((m * Bo[1] / Bg[1] + Rs[1]) * ((Bg - Bo * Rv) / (1 - Rv * Rs))) / denom)
   gor <- vector(length = l)
   RF_oil <- vector(length = l)
   RF_gas <- vector(length = l)
   Np <- vector(length = l)
   Gp <- vector(length = l)
   volfrac_liq <- vector(length = l)
   gor[1] <- Rs[1]
   RF_oil[1] <- 0
   RF_gas[1] <- 0
   Gp[1] <- RF_gas[1] * G
   swi_t <- swi
   soi_t <- soi * (1 - PVgas_r[1]) + sorg * PVgas_r[1]
   sgi_t <- 1 - swi_t - soi_t
   sw_t[1] <- swi_t
   so_t[1] <- soi_t
   sg_t[1] <- 1 - sw_t[1] - so_t[1]
   sw[1] <- swi
   so[1] <- 1 - sw[1]
   sg[1] <- 1 - sw[1] - so[1]
   volfrac_liq[1] <- 1.0
   for (i in 2:l) {
      if (p[i] >= pb) {
         gor_g <- gor[1]
      } else {
         gor_g <- 1.1 * gor[i - 1]
      }
      error <- 1e3
      while (abs(error) > 1e-6) {
         gor_m <- mean(c(gor_g, gor[i - 1]))
         dNp_N <- (1 - RF_oil[i - 1] * phi_o[i] - RF_gas[i - 1] * phi_g[i]) /
            (phi_o[i] + gor_m * (1 + m * Bo[1] * Rv[1] / Bg[1]) * phi_g[i] / (m * Bo[1] / Bg[1] + Rs[1]))
         RF_oil[i] <- dNp_N + RF_oil[i - 1]
         Np[i] <- RF_oil[i] * N
         dGp_G <- gor_m * dNp_N * (1 + m * Bo[1] * Rv[1] / Bg[1]) / (m * Bo[1] / Bg[1] + Rs[1])
         RF_gas[i] <- dGp_G + RF_gas[i - 1]
         Gp[i] <- RF_gas[i] * G
         sw_t[i] <- swi
         so_t[i] <- ((1 - RF_oil[i]) * (soi_t / Bo[1] + sgi_t * Rv[1] / Bg[1]) - (1 - sw_t[i]) * Rv[i] / Bg[i]) / (1 / Bo[i] - Rv[i] / Bg[i])
         sg_t[i] <- 1 - sw_t[i] - so_t[i]
         sw[i] <- swi
         so[i] <- ifelse(p[i] >= pb, 1 - sw[i], (so_t[i] - sorg * PVgas_r[i]) / (1 - PVgas_r[i]))
         sg[i] <- 1 - sw[i] - so[i]
         krg <- approx(x = rel_perm$Sg, y = rel_perm$Krg, xout = sg[i], rule = 2)$y
         kro <- approx(x = rel_perm$Sg, y = rel_perm$Krog, xout = sg[i], rule = 2)$y
         gor_g_n <- (kro * Rs[i] / muo[i] / Bo[i] + krg / mug[i] / Bg[i]) / (kro / muo[i] / Bo[i] + krg * Rv[i] / mug[i] / Bg[i])
         error <- gor_g_n - gor_g
         gor_g <- gor_g_n
         gor[i] <- gor_g
      }
      volfrac_liq[i] <- 1 / (1 + Bg[i] * (Rs[1] - Rs[i]) / (Bo[i] * (1 - Rs[1] * Rv[i])))
   }
   Igd <- ifelse(dplyr::row_number(Igd) == 1, NA, Gfgi * Eg / Et)
   Isd <- ifelse(dplyr::row_number(Isd) == 1, NA, Nfoi * Eo / Et)
   Inwd <- ifelse(dplyr::row_number(Inwd) == 1, NA, 0)
   Ifwd <- ifelse(dplyr::row_number(Ifwd) == 1, NA, (W * Ew + PV * Ef) / Et)
   Iawd <- ifelse(dplyr::row_number(Iawd) == 1, NA, 0)
   Itot <- Igd + Isd + Inwd + Ifwd + Iawd
   names <- c("P (psia)", "SOo", "SGo", "SWo", "SOT", "SGT", "SWT", "GOR (SCF/STB)", "RF_oil", "RF_gas", "Liq_volume", "Igd", "Isd", "Inwd", "Ifwd", "Iawd", "Itot")
   results <- data.frame(p = p,  SO = so, SG = sg, SW = sw, SOT = so_t, SGT = sg_t, SWT = sw_t, gor = gor, RF_oil = RF_oil, RF_gas = RF_gas, volfrac_liq = volfrac_liq, Igd = Igd, Isd = Isd, Inwd = Inwd, Ifwd = Ifwd, Iawd = Iawd, Itot = Itot)
   colnames(results) <- names
   return(results)
}


