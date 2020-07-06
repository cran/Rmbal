

#***********************mbal_param_gas******************************************

#' A list object of class 'mbal_gas' for material balance analysis
#'
#' Create an object of class 'mbal_gas'
#'
#' @param input_unit a unit system for parameters, only the character string 'Field' is accepted
#' @param output_unit a unit system for properties, only the character string 'Field' is accepted
#' @param aquifer_model defaulted to `NULL`, otherwise must be a character string, one of the following eight options: 'uss_rad_edge', 'uss_rad_bottom', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom', 'pot'. For further information about each model, please see 'Raquifer' package reference manual (https://cran.r-project.org/web/packages/Raquifer/index.html)
#' @param G original gas in place, SCF.
#' @param phi reservoir porosity, a numeric fraction
#' @param swi initial water saturation in the reservoir, a numeric fraction
#' @param Np cumulative oil production, STB
#' @param Gp cumulative gas production, SCF
#' @param Wp cumulative water production, STB
#' @param Wi cumulative water injection, STB
#' @param We cumulative aquifer water influx, BBL. If unknown, a `NULL` value must be assigned
#' @param pd dew point pressure, a numeric value, psi
#' @param p reservoir pressure, a numeric vector, psi
#' @param pvt a data frame of PVT properties including pressure 'p' in 'psi', oil formation volume factor 'Bo' in 'bbl/stb', solution gas-oil ratio 'Rs' in 'scf/stb', oil viscosity 'muo' in 'cp', volatilized oil-gas ratio 'Rv' in 'stb/scf', gas formation volume factor 'Bg' in 'bbl/scf', gas viscosity 'mug' in 'cp', water formation volume factor 'Bw' in 'bbl/stb', and water viscosity 'muw' in 'cp'
#' @param cf formation compressibility, a numeric value or vector, 1/psi
#' @param M ratio of non-net-pay pore volume to the reservoir (net-pay) volume, a numeric fraction.
#' @param phi_a aquifer porosity, a numeric fraction
#' @param perm_h_a aquifer horizontal permeability, md. Used in 'uss_rad_edge', 'uss_rad_bottom', 'uss_lin_edge', 'pss_rad_edge', 'pss_lin_edge' and 'pot' aquifer models
#' @param perm_v_a vertical permeability, md. Used in 'uss_rad_bottom', 'uss_lin_bottom', 'pss_rad_bottom', and 'pss_lin_bottom' aquifer models
#' @param h_a aquifer height, ft
#' @param r_a aquifer radius, ft. Used in 'uss_rad_edge', 'uss_rad_bottom', 'pss_rad_edge', and 'pot' aquifer models
#' @param r_R reservoir radius, ft. Used in 'uss_rad_edge', 'uss_rad_bottom', 'pss_rad_edge', and 'pot' aquifer models
#' @param w_a aquifer width, ft. Used in 'uss_lin_edge', 'uss_lin_bottom', 'pss_lin_edge', and 'pss_lin_bottom' aquifer models
#' @param l_a aquifer length, ft. Used in 'uss_lin_edge', 'uss_lin_bottom', 'pss_lin_edge', and 'pss_lin_bottom' aquifer models
#' @param tetha fraction of reservoir encircled by the aquifer, degrees. Used in 'uss_rad_edge', 'pss_rad_edge', and 'pot' aquifer models
#' @param muw_a aquifer water viscosity, cp
#' @param cw_a aquifer water compressibility, a numeric value, 1/psi
#' @param cf_a aquifer formation compressibility, a numeric value, 1/psi
#' @param wf weight factor, a numeric vector of zeros and ones. A zero value excludes the entire row of reservoir history data at a particular time from the material balance analysis
#' @param sgrw residual gas saturation in water invaded zone (aquifer encroachment or water injection), a numeric fraction
#'
#' @return a list of class ’mbal_gas’ with all the required parameters for the mbal_perform_gas() S3 methods
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
#' p <- c(3700, 3650, 3400, 3100, 2800, 2500, 2200, 1900, 1600, 1300, 1000, 700,
#' 600)
#'
#' We <- rep(0, length.out = length(p))
#'
#' Np <- c(0, 28.6, 93, 231, 270, 379, 481, 517.2, 549, 580, 675, 755, 803) *1e3
#'
#' Gp <- c(0, 0.34, 1.2, 3.3, 4.3, 6.6, 9.1, 10.5, 12, 12.8, 16.4, 19.1, 20.5) * 1e9
#'
#' Wp <- rep(0, length.out = length(p))
#'
#' Wi <- rep(0, length.out = length(p))
#'
#' wf <- rep(1, length.out = length(p))
#'
#' mbal_param_gas_lst <- mbal_perform_param_gas(input_unit = "Field",
#' output_unit = "Field", G = 2.41e10, aquifer_model = NULL,
#' phi = 0.1, swi = 0.2, Np = Np, Gp = Gp, Wp = Wp, Wi = Wi, We = We, pd = 3650,
#' p = p, pvt = pvt_table, M = 0, cf = 2e-6, wf = wf, sgrw = 0.15)
#'
#' dplyr::glimpse(mbal_param_gas_lst)

mbal_perform_param_gas <- function(input_unit = "Field", output_unit = "Field", aquifer_model = NULL, G = NULL, phi = NULL, swi = NULL, Gp = NULL, Np = NULL, Wp = NULL, Wi = NULL, We = NULL, pd = NULL, p = NULL, pvt = NULL, cf = NULL, M = NULL, phi_a = NULL, perm_h_a = NULL, perm_v_a = NULL, h_a = NULL, r_a = NULL, r_R = NULL, w_a = NULL, l_a = NULL, tetha = NULL, muw_a = NULL, cw_a = NULL, cf_a = NULL, wf = NULL, sgrw = NULL) {

   if (!is.character(input_unit)) stop("'input_unit' must be the character string 'Field'.")
   if (input_unit != "Field") stop("'input_unit' must be the character string 'Field'.")
   if (!is.character(output_unit)) stop("'output_unit' must be the character string 'Field'.")
   if (output_unit != "Field") stop("'output_unit' must be the character string 'Field'.")
   if (!is.null(aquifer_model)) {
      if (!is.character(aquifer_model)) stop("'aquifer_model' must be a character string, one of the following eight options: 'uss_rad_edge', 'uss_rad_bottom', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom', 'pot'.")
      if (!(aquifer_model %in% c('uss_rad_edge', 'uss_rad_bottom', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom', 'pot'))) stop("'aquifer_model' must be one of the following eight options: 'uss_rad_edge', 'uss_rad_bottom', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom', 'pot'.")
      if (length(aquifer_model) > 1) stop("Only one 'aquifer-model' is acceptable for the material balance analysis.")
   }
   if (is.null(G)) stop("'G' must be a numeric value for 'gas' reservoirs.")
   if (!is.numeric(G)) stop("'G' must be a numeric value for 'gas' reservoirs.")
   if (length(G) != 1) stop("'G' must be a numeric value for 'gas' reservoirs.")
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
   if (is.null(Np)) stop("'Np' must be a numeric vector.")
   if (!is.numeric(Np)) stop("'Np' must be a numeric vector.")
   if (Np[1] != 0) stop("First reported 'Np' value must be zero.")
   if (is.null(Gp)) stop("'Gp' must be a numeric vector.")
   if (!is.numeric(Gp)) stop("'Gp' must be a numeric vector.")
   if (Gp[1] != 0) stop("First reported 'Gp' value must be zero.")
   if (is.null(pd)) stop("'pd' must be a numeric value.")
   if (!is.numeric(pd)) stop("'pd' must be a numeric value.")
   if (length(pd) != 1) stop("'pd' must be a numeric value.")
   if (is.null(p)) stop("'p' must be a numeric vector.")
   if (!is.numeric(p)) stop("'p' must be a numeric vector.")
   l <- length(p)
   if (p[1] < pd) stop("Initial reservoir pressure must be equal to or greater than 'pd'.")
   if (length(Np) != l) stop("'p' and 'Np' vectors must have the same length.")
   if (length(Gp) != l) stop("'p' and 'Gp' vectors must have the same length.")
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
      if (!('muw' %in% colnames(pvt))) {
         stop("Column 'muw' is missing in the 'pvt' data frame.")
      }
      if (max(p) > max(pvt$p)) {
         stop("Pressure range in the 'pvt' data frame does not cover the entire range of 'p' vector.")
      }
      if (min(p) < min(pvt$p)) {
         stop("Pressure range in the 'pvt' data frame does not cover the entire range of 'p' vector.")
      }
      Rv_pvt_max <- (pvt[pvt$p == max(pvt$p),]$Rv)[1]
      for (i in 1:l) {
         if (i == 1) {
            diff <- 0
         } else {
            if (p[i] >= pd) {
               diff <- abs(Rv_pvt_max - Np[i] / Gp[i])
               if (diff > 1e-9) stop("'Np/Gp' values are not equal to the 'Rv' values in the 'pvt' data frame at pressures equal to or greater than 'pd'.")
            }
         }
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
   if (is.null(Wp)) {
      Wp <- rep(0, l)
   }
   if (!is.null(Wp)) {
      if (!is.numeric(Wp)) stop("'Wp' must be a numeric vector.")
      if (is.numeric(Wp)) {
         if (length(Wp) != l) stop("'p' and 'Wp' vectors must have the same length.")
      }
      if (Wp[1] != 0) stop("First reported 'Wp' value must be zero.")
   }
   if (is.null(Wi)) {
      Wi <- rep(0, l)
   }
   if (!is.null(Wi)) {
      if (!is.numeric(Wi)) stop("'Wi' must be a numeric vector.")
      if (is.numeric(Wi)) {
         if (length(Wi) != l) stop("'p' and 'Wi' vectors must have the same length.")
      }
      if (Wi[1] != 0) stop("First reported 'Wi' value must be zero.")
   }
   if (is.null(aquifer_model)) {
      if (is.null(We)) {
         We <- rep(0, l)
         aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, We = We)
         class(aquifer_lst) <- c("NoA", "aquifer")
      } else {
         if (!is.numeric(We)) stop("'We' must be a numeric vector.")
         if (is.numeric(We)) {
            if (length(We) != l) stop("'p' and 'We' vectors must have the same length.")
         }
         aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, We = We)
         class(aquifer_lst) <- c("We", "aquifer")
      }
   } # end of if (is.null(aquifer_model)) {


   if (!is.null(aquifer_model)) {
      if (aquifer_model == "pot") {
         mult_len <- c(1)
      }
      if (aquifer_model %in% c('uss_rad_edge', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom')) {
         mult_len <- c(1, 1)
      }
      if (aquifer_model == 'uss_rad_bottom') {
         mult_len <- c(1, 1, 1)
      }
      if (!is.numeric(phi_a)) stop("'phi_a' must be a numeric value between zero and one.")
      if (phi_a >= 1 | phi_a <= 0) stop("'phi_a' must be a numeric value between zero and one.")
      if (aquifer_model == "uss_rad_edge") {
         if (!is.numeric(perm_h_a)) stop("'perm_h_a' must be a numeric.")
         if (!is.numeric(h_a)) stop("'h_a' must be a numeric.")
         if (!is.numeric(r_a)) stop("'r_a' must be a numeric.")
         if (!is.numeric(r_R)) stop("'r_R' must be a numeric.")
         if (!is.numeric(tetha)) stop("'tetha' must be a numeric.")
         if (!is.numeric(muw_a)) stop("'muw_a' must be a numeric.")
         if (!is.numeric(cw_a)) stop("'cw_a' must be a numeric.")
         if (length(cw_a) != 1) stop("'cw_a' must be a numeric.")
         if (!is.numeric(cf_a)) stop("'cf_a' must be a numeric.")
         if (length(cf_a) != 1) stop("'cf_a' must be a numeric.")
         aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, model = "uss_rad_edge", phi = phi_a,
                             perm_h = perm_h_a, h_a = h_a, r_a = r_a, r_R = r_R, tetha = tetha, mu_water = muw_a,
                             c_water = cw_a, c_rock = cf_a, pressure = p, mult_len = mult_len)
         class(aquifer_lst) <- c("uss_rad_edge", "aquifer")
      }
      if (aquifer_model == "uss_rad_bottom") {
         if (!is.numeric(perm_h_a)) stop("'perm_h_a' must be a numeric.")
         if (!is.numeric(perm_v_a)) stop("'perm_v_a' must be a numeric.")
         if (!is.numeric(h_a)) stop("'h_a' must be a numeric.")
         if (!is.numeric(r_a)) stop("'r_a' must be a numeric.")
         if (!is.numeric(r_R)) stop("'r_R' must be a numeric.")
         if (!is.numeric(muw_a)) stop("'muw_a' must be a numeric.")
         if (!is.numeric(cw_a)) stop("'cw_a' must be a numeric.")
         if (length(cw_a) != 1) stop("'cw_a' must be a numeric.")
         if (!is.numeric(cf_a)) stop("'cf_a' must be a numeric.")
         if (length(cf_a) != 1) stop("'cf_a' must be a numeric.")
         aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, model = "uss_rad_bottom", phi = phi_a,
                             perm_h = perm_h_a, perm_v = perm_v_a, h_a = h_a, r_a = r_a, r_R = r_R, mu_water = muw_a,
                             c_water = cw_a, c_rock = cf_a, pressure = p, mult_len = mult_len)
         class(aquifer_lst) <- c("uss_rad_bottom", "aquifer")
      }
      if (aquifer_model == "pss_rad_edge") {
         if (!is.numeric(perm_h_a)) stop("'perm_h_a' must be a numeric.")
         if (!is.numeric(h_a)) stop("'h_a' must be a numeric.")
         if (!is.numeric(r_a)) stop("'r_a' must be a numeric.")
         if (!is.numeric(r_R)) stop("'r_R' must be a numeric.")
         if (!is.numeric(tetha)) stop("'tetha' must be a numeric.")
         if (!is.numeric(muw_a)) stop("'muw_a' must be a numeric.")
         if (!is.numeric(cw_a)) stop("'cw_a' must be a numeric.")
         if (length(cw_a) != 1) stop("'cw_a' must be a numeric.")
         if (!is.numeric(cf_a)) stop("'cf_a' must be a numeric.")
         if (length(cf_a) != 1) stop("'cf_a' must be a numeric.")
         aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, model = "pss_rad_edge", phi = phi_a,
                             perm_h = perm_h_a, h_a = h_a, r_a = r_a, r_R = r_R, tetha = tetha, mu_water = muw_a,
                             c_water = cw_a, c_rock = cf_a, pressure = p, mult_len = mult_len)
         class(aquifer_lst) <- c("pss_rad_edge", "aquifer")
      }
      if (aquifer_model == "uss_lin_edge") {
         if (!is.numeric(perm_h_a)) stop("'perm_h_a' must be a numeric.")
         if (!is.numeric(h_a)) stop("'h_a' must be a numeric.")
         if (!is.numeric(w_a)) stop("'w_a' must be a numeric.")
         if (!is.numeric(l_a)) stop("'l_a' must be a numeric.")
         if (!is.numeric(muw_a)) stop("'muw_a' must be a numeric.")
         if (!is.numeric(cw_a)) stop("'cw_a' must be a numeric.")
         if (length(cw_a) != 1) stop("'cw_a' must be a numeric.")
         if (!is.numeric(cf_a)) stop("'cf_a' must be a numeric.")
         if (length(cf_a) != 1) stop("'cf_a' must be a numeric.")
         aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, model = "uss_lin_edge", phi = phi_a,
                             perm_h = perm_h_a, h_a = h_a, w_a = w_a, l_a = l_a, mu_water = muw_a,
                             c_water = cw_a, c_rock = cf_a, pressure = p, mult_len = mult_len)
         class(aquifer_lst) <- c("uss_lin_edge", "aquifer")
      }
      if (aquifer_model == "uss_lin_bottom") {
         if (!is.numeric(perm_v_a)) stop("'perm_v_a' must be a numeric.")
         if (!is.numeric(h_a)) stop("'h_a' must be a numeric.")
         if (!is.numeric(w_a)) stop("'w_a' must be a numeric.")
         if (!is.numeric(l_a)) stop("'l_a' must be a numeric.")
         if (!is.numeric(muw_a)) stop("'muw_a' must be a numeric.")
         if (!is.numeric(cw_a)) stop("'cw_a' must be a numeric.")
         if (length(cw_a) != 1) stop("'cw_a' must be a numeric.")
         if (!is.numeric(cf_a)) stop("'cf_a' must be a numeric.")
         if (length(cf_a) != 1) stop("'cf_a' must be a numeric.")
         aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, model = "uss_lin_bottom", phi = phi_a,
                             perm_v = perm_v_a, h_a = h_a, w_a = w_a, l_a = l_a, mu_water = muw_a,
                             c_water = cw_a, c_rock = cf_a, pressure = p, mult_len = mult_len)
         class(aquifer_lst) <- c("uss_lin_bottom", "aquifer")
      }
      if (aquifer_model == "pss_lin_edge") {
         if (!is.numeric(perm_h_a)) stop("'perm_h_a' must be a numeric.")
         if (!is.numeric(h_a)) stop("'h_a' must be a numeric.")
         if (!is.numeric(w_a)) stop("'w_a' must be a numeric.")
         if (!is.numeric(l_a)) stop("'l_a' must be a numeric.")
         if (!is.numeric(muw_a)) stop("'muw_a' must be a numeric.")
         if (!is.numeric(cw_a)) stop("'cw_a' must be a numeric.")
         if (length(cw_a) != 1) stop("'cw_a' must be a numeric.")
         if (!is.numeric(cf_a)) stop("'cf_a' must be a numeric.")
         if (length(cf_a) != 1) stop("'cf_a' must be a numeric.")
         aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, model = "pss_lin_edge", phi = phi_a,
                             perm_h = perm_h_a, h_a = h_a, w_a = w_a, l_a = l_a, mu_water = muw_a,
                             c_water = cw_a, c_rock = cf_a, pressure = p, mult_len = mult_len)
         class(aquifer_lst) <- c("pss_lin_edge", "aquifer")
      }
      if (aquifer_model == "pss_lin_bottom") {
         if (!is.numeric(perm_v_a)) stop("'perm_v_a' must be a numeric.")
         if (!is.numeric(h_a)) stop("'h_a' must be a numeric.")
         if (!is.numeric(w_a)) stop("'w_a' must be a numeric.")
         if (!is.numeric(l_a)) stop("'l_a' must be a numeric.")
         if (!is.numeric(muw_a)) stop("'muw_a' must be a numeric.")
         if (!is.numeric(cw_a)) stop("'cw_a' must be a numeric.")
         if (length(cw_a) != 1) stop("'cw_a' must be a numeric.")
         if (!is.numeric(cf_a)) stop("'cf_a' must be a numeric.")
         if (length(cf_a) != 1) stop("'cf_a' must be a numeric.")
         aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, model = "pss_lin_bottom", phi = phi_a,
                             perm_v = perm_v_a, h_a = h_a, w_a = w_a, l_a = l_a, mu_water = muw_a,
                             c_water = cw_a, c_rock = cf_a, pressure = p, mult_len = mult_len)
         class(aquifer_lst) <- c("pss_lin_bottom", "aquifer")
      }
      if (aquifer_model == "pot") {
         if (!is.numeric(h_a)) stop("'h_a' must be a numeric.")
         if (!is.numeric(r_a)) stop("'r_a' must be a numeric.")
         if (!is.numeric(r_R)) stop("'r_R' must be a numeric.")
         if (!is.numeric(tetha)) stop("'tetha' must be a numeric.")
         if (!is.numeric(cw_a)) stop("'cw_a' must be a numeric.")
         if (length(cw_a) != 1) stop("'cw_a' must be a numeric.")
         if (!is.numeric(cf_a)) stop("'cf_a' must be a numeric.")
         if (length(cf_a) != 1) stop("'cf_a' must be a numeric.")
         aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, model = "pot", phi = phi_a,
                             h_a = h_a, r_a = r_a, r_R = r_R, tetha = tetha, c_water = cw_a,
                             c_rock = cf_a, pressure = p, mult_len = mult_len)
         class(aquifer_lst) <- c("pot", "aquifer")
      }
   } # end of if (!is.null(aquifer_model)) {
   if (class(aquifer_lst)[1] == "NoA") {
      if (is.null(sgrw)) {
         if (any(Wi > 0)) stop("A numeric value must be assigned to 'sgrw' for a water injection case.")
      }
   }
   if (class(aquifer_lst)[1] != "NoA") {
      if (is.null(sgrw)) stop("A numeric value must be assigned to 'sgrw' for a 'water_drive' reservoir.")
   }
   if (is.null(wf)) {
      wf <- rep(1,l)
   } else {
      if (!is.numeric(wf)) stop("'Wf' must be a numeric vector.")
      if (length(wf) != l) stop("'p' and 'Wf' vectors must have the same length.")
      wf[wf != 0] <- 1
   }
   pvt <- pvt
   prod <- data.frame(Np = Np, Gp = Gp, Wp = Wp)
   inj <- data.frame(Wi = Wi)
   mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, G = G, phi = phi, swi = swi, pd = pd, p = p, cf = cf, M = M, pvt = pvt, prod = prod, inj = inj, aquifer = aquifer_lst, wf = wf, sgrw = sgrw)

   if (class(aquifer_lst)[1] == "NoA") {
      class(mbal_lst) <- c("volumetric_gas", "mbal_gas")
   } else {
      class(mbal_lst) <- c("water_drive_gas", "mbal_gas")
   }
   return(mbal_lst)
}



# ******************************************************************************

#' Generic function for performance predictions for a gas reservoir
#'
#' Generate a data frame of reservoir performance data according to the class of 'mbal_lst' and 'time_lst' objects
#'
#' @param mbal_lst a list object of class 'mbal_gas'
#' @param time_lst a list object of class 'time/date'
#'
#' @return a data frame with estimates for fluids saturation, drive indices, production rates, and gas-oil ratios over the pressure history of the reservoir
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
#' p <- c(3700, 3650, 3400, 3100, 2800, 2500, 2200, 1900, 1600, 1300, 1000, 700,
#' 600)
#'
#' We <- rep(0, length.out = length(p))
#'
#' Np <- c(0, 28.6, 93, 231, 270, 379, 481, 517.2, 549, 580, 675, 755, 803) *1e3
#'
#' Gp <- c(0, 0.34, 1.2, 3.3, 4.3, 6.6, 9.1, 10.5, 12, 12.8, 16.4, 19.1, 20.5) * 1e9
#'
#' Wp <- rep(0, length.out = length(p))
#'
#' Wi <- rep(0, length.out = length(p))
#'
#' wf <- rep(1, length.out = length(p))
#'
#' mbal_param_gas_lst <- mbal_perform_param_gas(input_unit = "Field",
#' output_unit = "Field", G = 2.41e10, aquifer_model = NULL,
#' phi = 0.1, swi = 0.2, Np = Np, Gp = Gp, Wp = Wp, Wi = Wi, We = We, pd = 3650,
#' p = p, pvt = pvt_table, M = 0, cf = 2e-6, wf = wf, sgrw = 0.15)
#'
#' time_lst <- mbal_time(c(1:length(p)), "year")
#'
#' mbal_results <- mbal_perform_gas(mbal_param_gas_lst, time_lst)
#'
#' dplyr::glimpse(mbal_results)

mbal_perform_gas <- function(mbal_lst, time_lst) {

   if (inherits(mbal_lst, "mbal_gas") == TRUE & inherits(time_lst, "time")) {
      UseMethod("mbal_perform_gas")
   } else {
      if (!inherits(mbal_lst, "mbal_gas")) {
         stop("A class of 'mbal_gas' must be assigned to the 'mbal_lst' parameter of the mbal_perform_gas() function.")
      }
      if (!inherits(time_lst, "time")) {
         stop("A class of 'time' must be assigned to the 'time_lst' parameter of the mbal_perform_gas() function.")
      }
   }
}


# ******************************************************************************

#' S3 method for class 'mbal_perform_gas'
#'
#' Return a data frame with estimates for fluids saturation, drive indices, production rates, and gas-oil ratios over the pressure history of a volumetric gas reservoir
#'
#' @param mbal_lst a list object of class 'mbal_gas'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame with estimates for fluids saturation, drive indices, production rates, and gas-oil ratios over the pressure history of a volumetric gas reservoir
#' @export
#'
mbal_perform_gas.volumetric_gas <- function(mbal_lst, time_lst) {
   volumetric_gas(mbal_lst, time_lst)
}


# ******************************************************************************

#' S3 method for class 'mbal_perform_gas'
#'
#' Return a data frame with estimates for fluids saturation, drive indices, production rates, and gas-oil ratios over the pressure history of a water_drive gas reservoir
#'
#' @param mbal_lst a list object of class 'mbal_gas'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame with estimates for fluids saturation, drive indices, production rates, and gas-oil ratios over the pressure history of a water_drive gas reservoir
#' @export
#'
mbal_perform_gas.water_drive_gas <- function(mbal_lst, time_lst) {
   water_drive_gas(mbal_lst, time_lst)
}






# ******************************************************************************

volumetric_gas <- function(mbal_lst, time_lst) {

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
   G <- mbal_lst$G
   phi <- mbal_lst$phi
   swi <- mbal_lst$swi
   pd <- mbal_lst$pd
   p <- mbal_lst$p
   l <- length(p)
   if (length(time_lst$t) != length(p)) stop("Lenght of 'time_lst$t' vector must be equal to the length of 'p' vector.")
   cf <- mbal_lst$cf
   M <- mbal_lst$M
   pvt <- mbal_lst$pvt
   prod <- mbal_lst$prod
   inj <- mbal_lst$inj
   aquifer <- mbal_lst$aquifer
   aqu_cls <- class(aquifer)[1]
   if (aqu_cls == "NoA") We <- aquifer$We
   wf <- mbal_lst$wf
   sgrw <- mbal_lst$sgrw
   if (is.null(sgrw)) stop("'sgrw' must be a numeric value.")
   keep <- which(wf == 1)
   p <- p[keep]
   time <- time[keep]
   if (length(cf) != 1) cf <- cf[keep]
   prod <- prod[keep,]
   inj <- inj[keep,]
   if (aqu_cls == "NoA") We <- We[keep]
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
   Np <- prod$Np
   Gp <- prod$Gp
   Wp <- prod$Wp
   Rpv <- ifelse(p >= pd, Rv, Np / Gp)
   Wi <- inj
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
   PVwater <- vector(length = l)
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
   Et <- Gfgi * Egwf + We + (Wi - Wp) * Bw
   F_ <- Gp * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw
   denom <- (Bo * (Rpv - Rv) + Bg * (1 - Rpv * Rs)) / (1 - Rv * Rs)
   RF_gas <- (((Wi - Wp) * Bw + We) / Gfgi + Egwf) / denom
   if (all(Rv == 0)) {
      RF_oil <- RF_gas
   } else {
      RF_oil <- Rpv * RF_gas / Rv[1]
   }
   Igd <- ifelse(dplyr::row_number(Igd) == 1, NA, Gfgi * Eg / Et)
   Inwd <- ifelse(dplyr::row_number(Inwd) == 1, NA, We / Et)
   Ifwd <- ifelse(dplyr::row_number(Ifwd) == 1, NA, (W * Ew + PV * Ef) / Et)
   Iawd <- ifelse(dplyr::row_number(Iawd) == 1, NA, (Wi - Wp) * Bw / Et)
   Itot <- Igd + Inwd + Ifwd + Iawd
   sgi <- 1 - swi
   soi <- 1 - swi - sgi
   PVwater_r <- ((Wi - Wp) * Bw + We) / (1 - sgrw - swi) / PV
   swi_t <- swi * (1 - PVwater_r[1]) + (1 - sgrw) * PVwater_r[1]
   sgi_t <- sgi * (1 - PVwater_r[1]) + sgrw * PVwater_r[1]
   soi_t <- 1 - swi_t - sgi_t
   sw_t <- swi * (1 - PVwater_r) + (1 - sgrw) * PVwater_r
   so_t <- ifelse(p >= pd, 0, (1 - RF_oil - Bg[1] * Rv / Bg / Rv[1]) * (1 - sw_t)  / (Bg[1] / Bo / Rv[1] - Bg[1] * Rv / Bg / Rv[1]))
   sg_t <- 1 - sw_t - so_t
   sw <- swi
   so <- (so_t - sgrw * PVwater_r) / (1 - PVwater_r)
   sg <- 1 - sw - so
   qo <- vector(length = l)
   qw <- vector(length = l)
   qg <- vector(length = l)
   qo_insitu <- vector(length = l)
   fw <- vector(length = l)
   fg <- vector(length = l)
   gor <- vector(length = l)
   krg_kro <- vector(length = l)
   qg <- ifelse(dplyr::row_number(Gp) == 1, 0, (Gp - dplyr::lag(Gp)) / (time - dplyr::lag(time)))
   qo <- ifelse(dplyr::row_number(Np) == 1, 0, (Np - dplyr::lag(Np)) / (time - dplyr::lag(time)))
   qo_insitu <- ifelse(dplyr::row_number(Np) == 1, 0, (Gp * (Rpv - Rv) - dplyr::lag(Gp) * (dplyr::lag(Rpv) - dplyr::lag(Rv))) / (time - dplyr::lag(time)))
   qo_insitu <- ifelse(qo_insitu < 0, 0, qo_insitu)
   qw <- ifelse(dplyr::row_number(Wp) == 1, 0, (Wp - dplyr::lag(Wp)) / (time - dplyr::lag(time)))
   fo <- ifelse(dplyr::row_number(qo_insitu) == 1, 0, qo_insitu * Bo / (qo_insitu * Bo + qg * Bg + qw * Bw))
   fw <- ifelse(dplyr::row_number(qw) == 1, 0, qw * Bw / (qo_insitu * Bo + qg * Bg + qw * Bw))
   gor <- ifelse(p >= pd, 1 / Rv, qg / qo)
   kro_krg <- 0
   names <- c("P (psia)", "Eo (bbl/STB)", "Eg (bbl/SCF)", "Ew (bbl/STB)", "Ef (bbl/bbl)", "Egwf (bbl/SCF)", "Et (bbl)", "F (bbl)", "We", "Igd", "Inwd", "Ifwd", "Iawd", "Itot", "RF_oil", "RF_gas", "SOg", "SGg", "SWg", "SOT", "SGT", "SWT", "qo (STB/day)", "qg (SCF/day)", "qw (STB/day)", "fo", "fw", "GOR (SCF/STB)", "kro/krg")
   results <- data.frame(p = p, Eo = Eo, Eg = Eg, Ew = Ew, Ef = Ef, Egwf = Egwf, Et = Et, F_ = F_, We = We, Igd = Igd, Inwd = Inwd, Ifwd = Ifwd, Iawd = Iawd, Itot = Itot, RF_oil = RF_oil, RF_gas = RF_gas, SO = so, SG = sg, SW = sw, SOT = so_t, SGT = sg_t, SWT = sw_t, qo = qo, qg = qg, qw = qw, fo = fo, fw = fw, gor = gor, kro_krg = kro_krg)
   colnames(results) <- names
   return(results)
}




water_drive_gas <- function(mbal_lst, time_lst) {

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
   G <- mbal_lst$G
   phi <- mbal_lst$phi
   swi <- mbal_lst$swi
   pd <- mbal_lst$pd
   p <- mbal_lst$p
   l <- length(p)
   if (length(time_lst$t) != length(p)) stop("Lenght of 'time_lst$t' vector must be equal to the length of 'p' vector.")
   cf <- mbal_lst$cf
   M <- mbal_lst$M
   pvt <- mbal_lst$pvt
   prod <- mbal_lst$prod
   inj <- mbal_lst$inj
   aquifer <- mbal_lst$aquifer
   aqu_cls <- class(aquifer)[1]
   if (aqu_cls == "We") We <- aquifer$We
   wf <- mbal_lst$wf
   sgrw <- mbal_lst$sgrw
   if (is.null(sgrw)) stop("'sgrw' must be a numeric value.")
   keep <- which(wf == 1)
   p <- p[keep]
   time <- time[keep]
   if (length(cf) != 1) cf <- cf[keep]
   prod <- prod[keep,]
   inj <- inj[keep,]
   if (aqu_cls == "We") We <- We[keep]
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
   Np <- prod$Np
   Gp <- prod$Gp
   Wp <- prod$Wp
   Rpv <- ifelse(p >= pd, Rv, Np / Gp)
   Wi <- inj
   if (aqu_cls == "We") {
      # We <- aquifer$We
   } else {
      We <- vector(length = l)
      if (aquifer$model %in% c("uss_rad_edge", "pss_rad_edge", "uss_rad_bottom")) {
         phi <- aquifer$phi
         perm_h <- aquifer$perm_h
         if (aquifer$model == "uss_rad_bottom") {
            perm_v <- aquifer$perm_v
         }
         h_a <- aquifer$h_a
         r_a <- aquifer$r_a
         r_R <- aquifer$r_R
         if (aquifer$model %in% c("uss_rad_edge", "pss_rad_edge")) {
            tetha <- aquifer$tetha
         }
         mu_water <- aquifer$mu_water
         c_water <- aquifer$c_water
         c_rock <- aquifer$c_rock
         pressure <- aquifer$pressure[keep]
         mult_len <- aquifer$mult_len
      }
      if (aquifer$model %in% c("uss_lin_edge", "pss_lin_edge", "uss_lin_bottom", "pss_lin_bottom")) {
         phi <- aquifer$phi
         if (aquifer$model %in% c("uss_lin_edge", "pss_lin_edge")) {
            perm_h <- aquifer$perm_h
         }
         if (aquifer$model %in% c("uss_lin_bottom", "pss_lin_bottom")) {
            perm_v <- aquifer$perm_v
         }
         h_a <- aquifer$h_a
         w_a <- aquifer$w_a
         l_a <- aquifer$l_a
         mu_water <- aquifer$mu_water
         c_water <- aquifer$c_water
         c_rock <- aquifer$c_rock
         pressure <- aquifer$pressure[keep]
         mult_len <- aquifer$mult_len
      }
      if (aquifer$model == "pot") {
         phi <- aquifer$phi
         h_a <- aquifer$h_a
         r_a <- aquifer$r_a
         r_R <- aquifer$r_R
         tetha <- aquifer$tetha
         c_water <- aquifer$c_water
         c_rock <- aquifer$c_rock
         pressure <- aquifer$pressure[keep]
         mult_len <- aquifer$mult_len
      }
      if (aquifer$model == "uss_rad_edge") {
         We <- veh_uss_rad_edge(phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure, mult_len)
      }
      if (aquifer$model == "uss_rad_bottom") {
         We <- yk_uss_rad_bottom(phi, perm_h, perm_v, h_a, r_a, r_R, mu_water, c_water, c_rock, time, pressure, mult_len)
      }
      if (aquifer$model == "pss_rad_edge") {
         We <- fetkovich_pss_rad_edge(phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure, mult_len)
      }
      if (aquifer$model == "uss_lin_edge") {
         We <- nb_uss_lin_edge(phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure, mult_len)
      }
      if (aquifer$model == "uss_lin_bottom") {
         We <- nb_uss_lin_bottom(phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure, mult_len)
      }
      if (aquifer$model == "pss_lin_edge") {
         We <- fetk_pss_lin_edge(phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure, mult_len)
      }
      if (aquifer$model == "pss_lin_bottom") {
         We <- fetk_pss_lin_bottom(phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure, mult_len)
      }
      if (aquifer$model == "pot") {
         We <- pot(phi, h_a, r_a, r_R, tetha, c_water, c_rock, pressure, mult_len)
      }
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
   PVwater <- vector(length = l)
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
   Et <- Gfgi * Egwf + We + (Wi - Wp) * Bw
   F_ <- (Gp) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw
   denom <- (Bo * (Rpv - Rv) + Bg * (1 - Rpv * Rs)) / (1 - Rv * Rs)
   RF_gas <- (((Wi - Wp) * Bw + We) / Gfgi + Egwf) / denom
   if (all(Rv == 0)) {
      RF_oil <- RF_gas
   } else {
      RF_oil <- Rpv * RF_gas / Rv[1]
   }
   Igd <- ifelse(dplyr::row_number(Igd) == 1, NA, Gfgi * Eg / Et)
   Inwd <- ifelse(dplyr::row_number(Inwd) == 1, NA, We / Et)
   Ifwd <- ifelse(dplyr::row_number(Ifwd) == 1, NA, (W * Ew + PV * Ef) / Et)
   Iawd <- ifelse(dplyr::row_number(Iawd) == 1, NA, (Wi - Wp) * Bw / Et)
   Itot <- Igd + Inwd + Ifwd + Iawd
   sgi <- 1 - swi
   soi <- 1 - swi - sgi
   PVwater_r <- ((Wi - Wp) * Bw + We) / (1 - sgrw - swi) / PV
   swi_t <- swi * (1 - PVwater_r[1]) + (1 - sgrw) * PVwater_r[1]
   sgi_t <- sgi * (1 - PVwater_r[1]) + sgrw * PVwater_r[1]
   soi_t <- 1 - swi_t - sgi_t
   sw_t <- swi * (1 - PVwater_r) + (1 - sgrw) * PVwater_r
   so_t <- ifelse(p >= pd, 0, (1 - RF_oil - Bg[1] * Rv / Bg / Rv[1]) * (1 - sw_t)  / (Bg[1] / Bo / Rv[1] - Bg[1] * Rv / Bg / Rv[1]))
   sg_t <- 1 - sw_t - so_t
   sw <- swi
   so <- (so_t - sgrw * PVwater_r) / (1 - PVwater_r)
   sg <- 1 - sw - so
   qo <- vector(length = l)
   qw <- vector(length = l)
   qg <- vector(length = l)
   qo_insitu <- vector(length = l)
   fw <- vector(length = l)
   fg <- vector(length = l)
   gor <- vector(length = l)
   krg_kro <- vector(length = l)
   qg <- ifelse(dplyr::row_number(Gp) == 1, 0, (Gp - dplyr::lag(Gp)) / (time - dplyr::lag(time)))
   qo <- ifelse(dplyr::row_number(Np) == 1, 0, (Np - dplyr::lag(Np)) / (time - dplyr::lag(time)))
   qo_insitu <- ifelse(dplyr::row_number(Np) == 1, 0, (Gp * (Rpv - Rv) - dplyr::lag(Gp) * (dplyr::lag(Rpv) - dplyr::lag(Rv))) / (time - dplyr::lag(time)))
   qo_insitu <- ifelse(qo_insitu < 0, 0, qo_insitu)
   qw <- ifelse(dplyr::row_number(Wp) == 1, 0, (Wp - dplyr::lag(Wp)) / (time - dplyr::lag(time)))
   fo <- ifelse(dplyr::row_number(qo_insitu) == 1, 0, qo_insitu * Bo / (qo_insitu * Bo + qg * Bg + qw * Bw))
   fw <- ifelse(dplyr::row_number(qw) == 1, 0, qw * Bw / (qo_insitu * Bo + qg * Bg + qw * Bw))
   gor <- ifelse(p >= pd, 1 / Rv, qg / qo)
   kro_krg <- 0
   names <- c("P (psia)", "Eo (bbl/STB)", "Eg (bbl/SCF)", "Ew (bbl/STB)", "Ef (bbl/bbl)", "Egwf (bbl/SCF)", "Et (bbl)", "F (bbl)", "We", "Igd", "Inwd", "Ifwd", "Iawd", "Itot", "RF_oil", "RF_gas", "SOg", "SGg", "SWg", "SOT", "SGT", "SWT", "qo (STB/day)", "qg (SCF/day)", "qw (STB/day)", "fo", "fw", "GOR (SCF/STB)", "kro/krg")
   results <- data.frame(p = p, Eo = Eo, Eg = Eg, Ew = Ew, Ef = Ef, Egwf = Egwf, Et = Et, F_ = F_, We = We, Igd = Igd, Inwd = Inwd, Ifwd = Ifwd, Iawd = Iawd, Itot = Itot, RF_oil = RF_oil, RF_gas = RF_gas, SO = so, SG = sg, SW = sw, SOT = so_t, SGT = sg_t, SWT = sw_t, qo = qo, qg = qg, qw = qw, fo = fo, fw = fw, gor = gor, kro_krg = kro_krg)
   colnames(results) <- names
   return(results)
}



