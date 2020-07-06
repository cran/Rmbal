
#' A list object of class 'optimization_gas' for material balance analysis
#'
#' Create an object of class 'optimization_gas'
#'
#' @param input_unit a unit system for parameters, only the character string 'Field' is accepted
#' @param output_unit a unit system for properties, only the character string 'Field' is accepted
#' @param unknown_param a character string showing the unknown parameter(s). One of the following options: 'G', 'We', 'M', or 'G_M'
#' @param aquifer_model defaulted to `NULL`, otherwise must be a character string, one of the following eight options: 'uss_rad_edge', 'uss_rad_bottom', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom', 'pot'. For further information about each model, please see 'Raquifer' package reference manual (https://cran.r-project.org/web/packages/Raquifer/index.html)
#' @param G original gas in place, SCF. If unknown, a `NULL` value must be assigned
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
#' @param M ratio of non-net-pay pore volume to the reservoir (net-pay) volume, a numeric fraction. If unknown, a `NULL` value must be assigned.
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
#' @param mult_len a numeric vector of initial estimates for the 'aquifer_model' parameters
#' A vector of length one for the 'pot' aquifer model. It applies as a multiplier to the radius of the aquifer
#' A vector of length two for the 'uss_rad_edge', and 'pss_rad_edge' aquifer models. The first parameter is applied as a multiplier to the aquifer radius, and the second parameter is applied as a multiplier to the aquifer horizontal permeability
#' A vector of length two for the 'uss_lin_edge', and 'pss_lin_edge' aquifer models. The first parameter is applied as a multiplier to the aquifer length, and the second parameter is applied as a multiplier to the aquifer horizontal permeability
#' A vector of length two for the 'uss_lin_bottom', and 'pss_lin_bottom' aquifer models. The first parameter is applied as a multiplier to the aquifer height, and the second parameter is applied as a multiplier to the aquifer vertical permeability
#' A vector of length three for the 'uss_rad_bottom' aquifer model. The first parameter is applied as a multiplier to the aquifer radius, the second parameter is applied as a multiplier to the aquifer horizontal permeability, and the third parameter is applied as a multiplier to the aquifer vertical permeability
#' @param lower an optional numeric vector of lower bounds for the 'aquifer_model' parameters. See 'minpack.lm' package for details
#' @param upper an optional numeric vector of upper bounds for the 'aquifer_model' parameters. See 'minpack.lm' package for details
#' @param control an optional list of control settings. See 'minpack.lm' package for details
#'
#' @return a list of class 'mbal_gas' with all the required parameters for the mbal_perform_gas() S3 methods
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
#' mbal_optim_gas_lst <- mbal_optim_param_gas(input_unit = "Field",
#' output_unit = "Field", unknown_param = "G", aquifer_model = NULL,
#' phi = 0.1, swi = 0.2, Np = Np, Gp = Gp, Wp = Wp, Wi = Wi, We = We, pd = 3650,
#' p = p, pvt = pvt_table, M = 0, cf = 2e-6, wf = wf, sgrw = 0.15)
#'
#' dplyr::glimpse(mbal_optim_gas_lst)

mbal_optim_param_gas <- function(input_unit = "Field", output_unit = "Field", unknown_param = NULL, aquifer_model = NULL, G = NULL, phi = NULL, swi = NULL, Np = NULL, Gp = NULL, Wp = NULL, Wi = NULL, We = NULL, pd = NULL, p = NULL, pvt = NULL, cf = NULL, M = NULL, phi_a = NULL, perm_h_a = NULL, perm_v_a = NULL, h_a = NULL, r_a = NULL, r_R = NULL, w_a = NULL, l_a = NULL, tetha = NULL, muw_a = NULL, cw_a = NULL, cf_a = NULL, wf = NULL, sgrw = NULL, mult_len = NULL, lower = NULL, upper = NULL, control = NULL) {

   if (!is.character(input_unit)) stop("'input_unit' must be the character string 'Field'.")
   if (input_unit != "Field") stop("'input_unit' must be the character string 'Field'.")
   if (!is.character(output_unit)) stop("'output_unit' must be the character string 'Field'.")
   if (output_unit != "Field") stop("'output_unit' must be the character string 'Field'.")
   if (!is.character(unknown_param))  stop("'unknown_param' must be one of the character strings from the following options: 'G', 'We', 'M', or 'G_M'.")
   if (is.null(unknown_param)) stop("One of the following options must be selected: 'G', 'We', 'M', or 'G_M'.")
   if (!(unknown_param %in% c('G', 'We', 'M', 'G_M'))) stop("One of the following options must be selected: 'G', 'We', 'M', or 'G_M'.")

   if (unknown_param == "G") {
      if (!is.null(G)) stop("'G' must be NULL for the 'unknown_param = G' case.")
      if (is.null(M)) stop("'M' must be a numeric for the 'unknown_param = G' case.")
      if (!is.numeric(M)) stop("'M' must be a numeric for the 'unknown_param = G' case.")
      if (M < 0) stop("'M' must be a numeric equal to or greater than zero for the 'unknown_param = G' case.")
      if (is.null(aquifer_model)) {
         if (is.null(We)) stop("Either 'We' or 'aquifer_model' must be known in advance for the unknown_param = 'G' case.")
      }
      if ((!is.null(aquifer_model)) & (!is.null(We))) {
         stop("Either 'We' or 'aquifer_model' must be NULL for the unknown_param = 'G' case.")
      }
   }
   if (unknown_param == "We") {
      if (is.null(aquifer_model)) stop("'aquifer_model' must be known for the 'unknown_param = We' case.")
      if (!is.null(We)) stop("'We' must be NULL for the 'unknown_param = We' case.")
      if (is.null(G)) stop("'G' must be a numeric for the 'unknown_param = We' case.")
      if (!is.numeric(G)) stop("'G' must be a numeric for the 'unknown_param = We' case.")
      if (G < 0) stop("'G' must be a numeric equal to or greater than zero for the 'unknown_param = We' case.")
      if (is.null(M)) stop("'M' must be a numeric for the 'unknown_param = We' case.")
      if (!is.numeric(M)) stop("'M' must be a numeric for the 'unknown_param = We' case.")
      if (M < 0) stop("'M' must be a numeric equal to or greater than zero for the 'unknown_param = We' case.")
   }
   if (unknown_param == "M") {
      if (!is.null(M)) stop("'M' must be NULL for the 'unknown_param = M' case.")
      if (is.null(G)) stop("'G' must be a numeric for the 'unknown_param = M' case.")
      if (!is.numeric(G)) stop("'G' must be a numeric for the 'unknown_param = M' case.")
      if (G < 0) stop("'G' must be a numeric equal to or greater than zero for the 'unknown_param = M' case.")
      if (is.null(aquifer_model)) {
         if (is.null(We)) stop("Either 'We' or 'aquifer_model' must be known in advance for the unknown_param = 'M' case.")
      }
      if ((!is.null(aquifer_model)) & (!is.null(We))) {
         stop("Either 'We' or 'aquifer_model' must be NULL for the unknown_param = 'M' case.")
      }
   }
   if (unknown_param == "G_M") {
      if (!is.null(G)) stop("'G' must be NULL for the 'unknown_param = G_M' case.")
      if (!is.null(M)) stop("'M' must be NULL for the 'unknown_param = G_M' case.")
      if (is.null(aquifer_model)) {
         if (is.null(We)) stop("Either 'We' or 'aquifer_model' must be known in advance for the unknown_param = 'G_M' case.")
      }
      if ((!is.null(aquifer_model)) & (!is.null(We))) {
         stop("Either 'We' or 'aquifer_model' must be NULL for the unknown_param = 'G_M' case.")
      }
   }
   if (!is.null(aquifer_model)) {
      if (!is.character(aquifer_model)) stop("'aquifer_model' must be a character string, one of the following eight options: 'uss_rad_edge', 'uss_rad_bottom', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom', 'pot'.")
      if (!(aquifer_model %in% c('uss_rad_edge', 'uss_rad_bottom', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom', 'pot'))) stop("'aquifer_model' must be one of the following eight options: 'uss_rad_edge', 'uss_rad_bottom', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom', 'pot'.")
      if (length(aquifer_model) > 1) stop("Only one 'aquifer-model' is acceptable for the material balance analysis.")
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
      if (!is.numeric(We)) stop("'We' must be a numeric vector.")
      if (is.numeric(We)) {
         if (length(We) != l) stop("'p' and 'We' vectors must have the same length.")
      }
      if (any(We != 0)) {
         aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, We = We)
         class(aquifer_lst) <- c("We", "aquifer")
      } else {
         aquifer_lst <- list(input_unit = input_unit, output_unit = output_unit, We = We)
         class(aquifer_lst) <- c("NoA", "aquifer")
      }
   } # end of if (is.null(aquifer_model)) {

   if (!is.null(aquifer_model)) {
      if (is.null(mult_len)) {
         if (aquifer_model == "pot") {
            mult_len <- c(1)
         }
         if (aquifer_model %in% c('uss_rad_edge', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom')) {
            mult_len <- c(1, 1)
         }
         if (aquifer_model == 'uss_rad_bottom') {
            mult_len <- c(1, 1, 1)
         }
      }
      if (!is.null(mult_len)) {
         if (aquifer_model == "pot") {
            if (length(mult_len) != 1) {
               stop("'mult_len' must be a numeric vector of length 1.")
            }
         }
         if (aquifer_model %in% c('uss_rad_edge', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom')) {
            if (length(mult_len) != 2) {
               stop("'mult_len' must be a numeric vector of length 2.")
            }
         }
         if (aquifer_model == 'uss_rad_bottom') {
            if (length(mult_len) != 3) {
               stop("'mult_len' must be a numeric vector of length 3.")
            }
         }
      }
      if (!is.null(lower)) {
         if (aquifer_model == "pot") {
            if (length(lower) != 1) {
               stop("'lower' must be a numeric vector of length 1.")
            }
         }
         if (aquifer_model %in% c('uss_rad_edge', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom')) {
            if (length(lower) != 2) {
               stop("'lower' must be a numeric vector of length 2.")
            }
         }
         if (aquifer_model == 'uss_rad_bottom') {
            if (length(lower) != 3) {
               stop("'lower' must be a numeric vector of length 3.")
            }
         }
      }
      if (!is.null(upper)) {
         if (aquifer_model == "pot") {
            if (length(upper) != 1) {
               stop("'upper' must be a numeric vector of length 1.")
            }
         }
         if (aquifer_model %in% c('uss_rad_edge', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom')) {
            if (length(upper) != 2) {
               stop("'upper' must be a numeric vector of length 2.")
            }
         }
         if (aquifer_model == 'uss_rad_bottom') {
            if (length(upper) != 3) {
               stop("'upper' must be a numeric vector of length 3.")
            }
         }
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
                             c_water = cw_a, c_rock = cf_a, pressure = p, mult_len = mult_len, lower = lower, upper = upper, control = control)
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
                             c_water = cw_a, c_rock = cf_a, pressure = p, mult_len = mult_len, lower = lower, upper = upper, control = control)
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
                             c_water = cw_a, c_rock = cf_a, pressure = p, mult_len = mult_len, lower = lower, upper = upper, control = control)
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
                             c_water = cw_a, c_rock = cf_a, pressure = p, mult_len = mult_len, lower = lower, upper = upper, control = control)
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
                             c_water = cw_a, c_rock = cf_a, pressure = p, mult_len = mult_len, lower = lower, upper = upper, control = control)
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
                             c_water = cw_a, c_rock = cf_a, pressure = p, mult_len = mult_len, lower = lower, upper = upper, control = control)
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
                             c_water = cw_a, c_rock = cf_a, pressure = p, mult_len = mult_len, lower = lower, upper = upper, control = control)
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
                             c_rock = cf_a, pressure = p, mult_len = mult_len, lower = lower, upper = upper, control = control)
         class(aquifer_lst) <- c("pot", "aquifer")
      }
   } # end of if (!is.null(aquifer_model)) {
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
   if (unknown_param == "G") {
      if (is.null(aquifer_model)) {
         if (all(We == 0)) res_type <- "volumetric"
         if (any(We != 0)) res_type <- "water_drive"
      } else {
         res_type <- "water_drive"
      }
   }
   if (unknown_param == "We") {
      res_type <- "water_drive"
   }
   if (unknown_param == "M") {
      if (is.null(aquifer_model)) {
         if (all(We == 0)) res_type <- "volumetric"
         if (any(We != 0)) res_type <- "water_drive"
      } else {
         res_type <- "water_drive"
      }
   }
   if (unknown_param == "G_M") {
      if (is.null(aquifer_model)) {
         if (all(We == 0)) res_type <- "volumetric"
         if (any(We != 0)) res_type <- "water_drive"
      } else {
         res_type <- "water_drive"
      }
   }
   if (res_type == "volumetric") {
      if (is.null(G) & is.null(M)) {
         optim_lst <- list(input_unit = input_unit, output_unit = output_unit, G = G, phi = phi, swi = swi, pd = pd, p = p, cf = cf, M = M, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sgrw = sgrw)
         class(optim_lst) <- c("volumetric_optim_gas", "optimization_gas")
      } else if (is.null(G)) {
         optim_lst <- list(input_unit = input_unit, output_unit = output_unit, G = G, phi = phi, swi = swi, pd = pd, p = p, cf = cf, M = M, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sgrw = sgrw)
         class(optim_lst) <- c("volumetric_M_optim_gas", "volumetric_optim_gas", "optimization_gas")
      } else {
         optim_lst <- list(input_unit = input_unit, output_unit = output_unit, G = G, phi = phi, swi = swi, pd = pd, p = p, cf = cf, M = M, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sgrw = sgrw)
         class(optim_lst) <- c("volumetric_G_optim_gas", "volumetric_optim_gas", "optimization_gas")
      }
   }
   if (res_type == "water_drive") {
      if (is.null(We)) {
         optim_lst <- list(input_unit = input_unit, output_unit = output_unit, G = G, phi = phi, swi = swi, pd = pd, p = p, cf = cf, M = M, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sgrw = sgrw)
         class(optim_lst) <- c("water_drive_G_M_optim_gas", "water_drive_optim_gas", "optimization_gas")
      } else {
         if (is.null(G) & is.null(M)) {
            optim_lst <- list(input_unit = input_unit, output_unit = output_unit, G = G, phi = phi, swi = swi, pd = pd, p = p, cf = cf, M = M, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sgrw = sgrw)
            class(optim_lst) <- c("water_drive_We_optim_gas", "water_drive_optim_gas", "optimization_gas")
         } else if(is.null(G)) {
            optim_lst <- list(input_unit = input_unit, output_unit = output_unit, G = G, phi = phi, swi = swi, pd = pd, p = p, cf = cf, M = M, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sgrw = sgrw)
            class(optim_lst) <- c("water_drive_M_We_optim_gas", "water_drive_optim_gas", "optimization_gas")
         } else {
            optim_lst <- list(input_unit = input_unit, output_unit = output_unit, G = G, phi = phi, swi = swi, pd = pd, p = p, cf = cf, M = M, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sgrw = sgrw)
            class(optim_lst) <- c("water_drive_G_We_optim_gas", "water_drive_optim_gas", "optimization_gas")
         }
      }
   }
   return(optim_lst)
}



# ******************************************************************************

#' Generic function for predicting unknown parameters of a material balance model
#'
#' Generate a list of class 'mbal_gas' with estimates for the unknown parameters of the material balance model according to the class of 'optim_lst' and 'time_lst' objects
#'
#' @param optim_lst a list object of class 'optimization_gas'
#' @param time_lst a list object of class 'time/date'
#'
#' @return a list of class 'mbal_gas' with estimates for the unknown parameters of the material balance model according to the class of 'optim_lst' and 'time_lst' objects
#'
#' @importFrom graphics abline mtext plot text
#' @importFrom stats approx lm
#' @importFrom pracma fzero trapz coth
#' @import gsl
#' @import minpack.lm
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
#' mbal_optim_gas_lst <- mbal_optim_param_gas(input_unit = "Field",
#' output_unit = "Field", unknown_param = "G", aquifer_model = NULL,
#' phi = 0.1, swi = 0.2, Np = Np, Gp = Gp, Wp = Wp, Wi = Wi, We = We, pd = 3650,
#' p = p, pvt = pvt_table, M = 0, cf = 2e-6, wf = wf, sgrw = 0.15)
#'
#' time_lst <- mbal_time(c(1:length(p)), "year")
#'
#' optim_results <- mbal_optim_gas(mbal_optim_gas_lst, time_lst)
#'
#' dplyr::glimpse(optim_results)

mbal_optim_gas <- function(optim_lst, time_lst) {

   if (inherits(optim_lst, "optimization_gas") == TRUE & inherits(time_lst, "time")) {
      UseMethod("mbal_optim_gas")
   } else {
      if (!inherits(optim_lst, "optimization_gas")) {
         stop("A class of 'optimization_gas' must be assigned to the 'optim_lst' parameter of the mbal_optim_gas() function.")
      }
      if (!inherits(time_lst, "time")) {
         stop("A class of 'time' must be assigned to the 'time_lst' parameter of the mbal_optim_gas() function.")
      }
   }
}


# ******************************************************************************

#' S3 method for class 'mbal_optim_gas'
#'
#' Generate a list of class 'mbal_gas' with estimates for the unknown parameters of a volumetric gas reservoir
#'
#' @param optim_lst a list object of class 'optimization_gas'
#' @param time_lst a list object of class 'time'
#'
#' @return a list of class 'mbal_gas' with estimates for the unknown parameters of a volumetric gas reservoir
#' @export
#'
mbal_optim_gas.volumetric_optim_gas <- function(optim_lst, time_lst) {
   volumetric_optim_gas(optim_lst, time_lst)
}


# ******************************************************************************

#' S3 method for class 'mbal_optim_gas'
#'
#' Generate a list of class 'mbal_gas' with estimates for the unknown parameters of a water_drive gas reservoir
#'
#' @param optim_lst a list object of class 'optimization_gas'
#' @param time_lst a list object of class 'time'
#'
#' @return a list of class 'mbal_gas' with estimates for the unknown parameters of a water_drive gas reservoir
#' @export
#'
mbal_optim_gas.water_drive_optim_gas <- function(optim_lst, time_lst) {
   water_drive_optim_gas(optim_lst, time_lst)
}


# ******************************************************************************

volumetric_optim_gas <- function(optim_lst, time_lst) {

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
   input_unit <- optim_lst$input_unit
   output_unit <- optim_lst$output_unit
   G <- optim_lst$G
   phi <- optim_lst$phi
   swi <- optim_lst$swi
   pd <- optim_lst$pd
   p <- optim_lst$p
   l <- length(p)
   if (length(time_lst$t) != length(p)) stop("Lenght of 'time_lst$t' vector must be equal to the length of 'p' vector.")
   cf <- optim_lst$cf
   M <- optim_lst$M
   pvt <- optim_lst$pvt
   prod <- optim_lst$prod
   inj <- optim_lst$inj
   We <- optim_lst$We
   aquifer <- optim_lst$aquifer
   cls <- class(optim_lst)[1]
   aqu_cls <- class(aquifer)[1]
   if (aqu_cls == "NoA") We <- aquifer$We
   wf <- optim_lst$wf
   sgrw <- optim_lst$sgrw
   if (is.null(sgrw)) stop("'sorw' must be a numeric value.")
   keep <- which(wf == 1)
   p_ <- p[keep]
   time_ <- time[keep]
   if (length(cf) != 1) {
      cf_ <- cf[keep]
   } else {
      cf_ <- cf
   }
   prod_ <- prod[keep,]
   inj_ <- inj[keep,]
   if (aqu_cls == "NoA") We_ <- We[keep]
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
      Bo[i] <- approx(pvt$p, pvt$Bo, xout = p_[i])$y
      Rs[i] <- approx(pvt$p, pvt$Rs, xout = p_[i])$y
      muo[i] <- approx(pvt$p, pvt$muo, xout = p_[i])$y
      Rv[i] <- approx(pvt$p, pvt$Rv, xout = p_[i])$y
      Bg[i] <- approx(pvt$p, pvt$Bg, xout = p_[i])$y
      mug[i] <- approx(pvt$p, pvt$mug, xout = p_[i])$y
      Bw[i] <- approx(pvt$p, pvt$Bw, xout = p_[i])$y
      muw[i] <- approx(pvt$p, pvt$muw, xout = p_[i])$y
   }
   Np <- prod_$Np
   Gp <- prod_$Gp
   Wp <- prod_$Wp
   Wi <- inj_
   dp <- vector(length = l)
   Bto <- vector(length = l)
   Btg <- vector(length = l)
   Eo <- vector(length = l)
   Eowf <- vector(length = l)
   Eg <- vector(length = l)
   Egwf <- vector(length = l)
   Ew <- vector(length = l)
   Ef <- vector(length = l)
   F_ <- vector(length = l)


   if (cls == "volumetric_optim_gas") {
      dp <- p_[1] - p_
      Bto <- (Bo * (1 - Rs[1] * Rv) + Bg * (Rs[1] - Rs)) / (1 - Rs * Rv)
      Btg <- (Bg * (1 - Rs * Rv[1]) + Bo * (Rv[1] - Rv)) / (1 - Rs * Rv)
      Eo <- Bto - Bto[1]
      Eg <- Btg - Btg[1]
      Ew <- Bw - Bw[1]
      for (i in 1:l) {
         Ef[i] <- -1 * trapz(p_[1:i], cf_[1:i])
      }
      A <- Eg + Bg[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      B <- Bg[1] * (Ew / Bw[1] + Ef) / (1 - swi)
      F_ <- Gp * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_
      data <- data.frame(LHS = F_, A = A, B = B)
      mbal_lm <- lm(LHS ~ 0 + A + B, data = data)
      Gfgi <- mbal_lm$coefficients[1]
      names(Gfgi) <- NULL
      M <- mbal_lm$coefficients[2] / Gfgi
      names(M) <- NULL
      Egwf <- Eg + Bg[1] * ((swi + M) * Ew / Bw[1] + (1 + M) * Ef) / (1 - swi)
      Nfoi <- 0
      N <- Nfoi + Rv[1] * Gfgi
      G <- Gfgi + Rs[1] * Nfoi
      PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
      BV <- PV / phi
      W <- PV * swi / Bw[1]
      Et <- Gfgi * Egwf + We_ + (Wi - Wp) * Bw
      STDERR <- sqrt(sum((F_ - Gfgi * Egwf) * ((F_ - Gfgi * Egwf))) / (l - 2))
      plot(Egwf, F_, pch = 21, bg = "blue" , xlab = "Egwf", ylab = "F")
      abline(a = 0, b = Gfgi)
      mtext(paste("slope = Gfgi = ", round(Gfgi,0), ", M = ", round(M,3), " , STDERR = ", round(STDERR,2)), side = 3)
      xmin <- min(Egwf, na.rm = TRUE) + 0.25 * diff(range(Egwf, na.rm = TRUE))
      ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
      text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F - G[fgi] * E[gwf])^2), n - 2))))
      p_calc <- vector(length = l)
      for (i in 1:l) {
         if (i == 1) {
            p_calc[i] <- p_[i]
         } else {
            fun <- function(x, pd, p_ini, Bg_ini, Btg_ini, Rs_ini, Rv_ini, Bw_ini, Gfgi, swi, p, Ef, Np, Gp, Wp, Wi, We, M, pvt) {
               Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear", ties = max)$y
               Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear", ties = max)$y
               Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear", ties = max)$y
               Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear", ties = max)$y
               Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear", ties = max)$y
               Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear", ties = max)$y
               dp <- p_ini - x
               Btg <- (Bg * (1 - Rs * Rv_ini) + Bo * (Rv_ini - Rv)) / (1 - Rs * Rv)
               Eg <- Btg - Btg_ini
               Ew <- Bw - Bw_ini
               Egwf <- Eg + Bg_ini * ((swi + M) * Ew / Bw_ini + (1 + M) * Ef_internal) / (1 - swi)
               F_lhs_calc <- Gp * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
               F_rhs_calc <- Gfgi * Egwf
               diff <- F_lhs_calc - F_rhs_calc
               return(diff)
            }
            p_calc[i] <- fzero(fun = fun, x = p_[i], pd = pd, p_ini = p_[1], Bg_ini = Bg[1], Btg_ini = Btg[1], Rs_ini = Rs[1], Rv_ini = Rv[1], Bw_ini = Bw[1], Gfgi = Gfgi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], We = We_[i], M = M, pvt = pvt, maxiter = 1000)$x
         }
      }
      STDERR <- sqrt(sum((p_ - p_calc) * ((p_ - p_calc)), na.rm = TRUE) / (l - 2))
      xmin <- min(c(p_,p_calc), na.rm = TRUE)
      xmax <- max(c(p_,p_calc), na.rm = TRUE)
      plot(p_, p_calc, pch = 21, bg = "blue" , xlab = "Pressure", ylab = "Pressure_estimated", xlim = c(xmin, xmax), ylim = c(xmin, xmax))
      abline(a = 0, b = 1)
      mtext(paste("slope = 1.0 " , "STDERR = ", round(STDERR,2)), side = 3)
      xrange <- xmax - xmin
      text(xmin + 0.25 * xrange, xmax - 0.25 * xrange, expression(STDERR == sqrt(frac(sum((p - p[est])^2), n - 2))))

      ymin <- max(0, min((F_[2:l]) / Egwf[2:l]) / 1.2)
      ymax <- 1.2 * max((F_[2:l]) / Egwf[2:l])
      plot(Gp, (F_)/ Egwf, pch = 21, bg = "blue" , xlab = "Gp", ylab = "F/Egwf", ylim = c(ymin, ymax))

      mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, G = G, phi = phi, swi = swi, pd = pd, p = p, cf = cf, M = M, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sgrw = sgrw, p_est = p_calc)
      class(mbal_lst) <- c("volumetric_gas", "mbal_gas")
      return(mbal_lst)
   }

   if (cls == "volumetric_M_optim_gas") {
      dp <- p_[1] - p_
      Bto <- (Bo * (1 - Rs[1] * Rv) + Bg * (Rs[1] - Rs)) / (1 - Rs * Rv)
      Btg <- (Bg * (1 - Rs * Rv[1]) + Bo * (Rv[1] - Rv)) / (1 - Rs * Rv)
      Eo <- Bto - Bto[1]
      Eg <- Btg - Btg[1]
      Ew <- Bw - Bw[1]
      for (i in 1:l) {
         Ef[i] <- -1 * trapz(p_[1:i], cf_[1:i])
      }
      Egwf <- Eg + Bg[1] * ((swi + M) * Ew / Bw[1] + (1 + M) * Ef) / (1 - swi)
      F_ <- Gp * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_
      data <- data.frame(Egwf = Egwf, `F_` = F_)
      mbal_lm <- lm(`F_` ~ 0 + Egwf, data = data)
      Gfgi <- mbal_lm$coefficients
      names(Gfgi) <- NULL
      Nfoi <- 0
      N <- Nfoi + Rv[1] * Gfgi
      G <- Gfgi + Rs[1] * Nfoi
      PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
      BV <- PV / phi
      W <- PV * swi / Bw[1]
      Et <- Gfgi * Egwf + We_ + (Wi - Wp) * Bw
      STDERR <- sqrt(sum((F_ - Gfgi * Egwf) * ((F_ - Gfgi * Egwf))) / (l - 2))
      plot(Egwf, F_, pch = 21, bg = "blue" , xlab = "Egwf", ylab = "F")
      abline(a = 0, b = Gfgi)
      mtext(paste("slope = Gfgi = ", round(Gfgi,0), " , STDERR = ", round(STDERR,2)), side = 3)
      xmin <- min(Egwf, na.rm = TRUE) + 0.25 * diff(range(Egwf, na.rm = TRUE))
      ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
      text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F - G[fgi] * E[gwf])^2), n - 2))))
      p_calc <- vector(length = l)
      for (i in 1:l) {
         if (i == 1) {
            p_calc[i] <- p_[i]
         } else {
            fun <- function(x, pd, p_ini, Bg_ini, Btg_ini, Rs_ini, Rv_ini, Bw_ini, Gfgi, swi, p, Ef, Np, Gp, Wp, Wi, We, M, pvt) {
               Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear", ties = max)$y
               Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear", ties = max)$y
               Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear", ties = max)$y
               Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear", ties = max)$y
               Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear", ties = max)$y
               Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear", ties = max)$y
               dp <- p_ini - x
               Btg <- (Bg * (1 - Rs * Rv_ini) + Bo * (Rv_ini - Rv)) / (1 - Rs * Rv)
               Eg <- Btg - Btg_ini
               Ew <- Bw - Bw_ini
               Egwf <- Eg + Bg_ini * ((swi + M) * Ew / Bw_ini + (1 + M) * Ef_internal) / (1 - swi)
               F_lhs_calc <- Gp * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
               F_rhs_calc <- Gfgi * Egwf
               diff <- F_lhs_calc - F_rhs_calc
               return(diff)
            }
            p_calc[i] <- fzero(fun = fun, x = p_[i], pd = pd, p_ini = p_[1], Bg_ini = Bg[1], Btg_ini = Btg[1], Rs_ini = Rs[1], Rv_ini = Rv[1], Bw_ini = Bw[1], Gfgi = Gfgi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], We = We_[i], M = M, pvt = pvt, maxiter = 1000)$x
         }
      }
      STDERR <- sqrt(sum((p_ - p_calc) * ((p_ - p_calc)), na.rm = TRUE) / (l - 2))
      xmin <- min(c(p_,p_calc), na.rm = TRUE)
      xmax <- max(c(p_,p_calc), na.rm = TRUE)
      plot(p_, p_calc, pch = 21, bg = "blue" , xlab = "Pressure", ylab = "Pressure_estimated", xlim = c(xmin, xmax), ylim = c(xmin, xmax))
      abline(a = 0, b = 1)
      mtext(paste("slope = 1.0 " , "STDERR = ", round(STDERR,2)), side = 3)
      xrange <- xmax - xmin
      text(xmin + 0.25 * xrange, xmax - 0.25 * xrange, expression(STDERR == sqrt(frac(sum((p - p[est])^2), n - 2))))

      ymin <- max(0, min((F_[2:l]) / Egwf[2:l]) / 1.2)
      ymax <- 1.2 * max((F_[2:l]) / Egwf[2:l])
      plot(Gp, F_ / Egwf, pch = 21, bg = "blue" , xlab = "Gp", ylab = "F/Egwf", ylim = c(ymin, ymax))

      mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, G = G, phi = phi, swi = swi, pd = pd, p = p, cf = cf, M = M, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sgrw = sgrw, p_est = p_calc)
      class(mbal_lst) <- c("volumetric_gas", "mbal_gas")
      return(mbal_lst)
   }
   if (cls == "volumetric_G_optim_gas") {
      dp <- p_[1] - p_
      Bto <- (Bo * (1 - Rs[1] * Rv) + Bg * (Rs[1] - Rs)) / (1 - Rs * Rv)
      Btg <- (Bg * (1 - Rs * Rv[1]) + Bo * (Rv[1] - Rv)) / (1 - Rs * Rv)
      Eo <- Bto - Bto[1]
      Eg <- Btg - Btg[1]
      Ew <- Bw - Bw[1]
      for (i in 1:l) {
         Ef[i] <- -1 * trapz(p_[1:i], cf_[1:i])
      }
      Gfgi <- G
      A <- Eg + Bg[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      B <- Bg[1] * (Ew / Bw[1] + Ef) / (1 - swi)
      F_ <- Gp * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_
      data <- data.frame(LHS = F_ / Gfgi - A, RHS = B)
      mbal_lm <- lm(LHS ~ 0 + RHS, data = data)
      M <- mbal_lm$coefficients
      names(M) <- NULL
      Egwf <- Eg + Bg[1] * ((swi + M) * Ew / Bw[1] + (1 + M) * Ef) / (1 - swi)
      Nfoi <- 0
      N <- Nfoi + Rv[1] * Gfgi
      PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
      BV <- PV / phi
      W <- PV * swi / Bw[1]
      Et <- Gfgi * Egwf + We_ + (Wi - Wp) * Bw
      STDERR <- sqrt(sum((F_ - Gfgi * Egwf) * ((F_ - Gfgi * Egwf))) / (l - 2))
      plot(Egwf, F_, pch = 21, bg = "blue" , xlab = "Egwf", ylab = "F")
      abline(a = 0, b = Gfgi)
      mtext(paste("slope = Gfgi = ", round(Gfgi,0), ", M = ", round(M,3), " , STDERR = ", round(STDERR,2)), side = 3)
      xmin <- min(Egwf, na.rm = TRUE) + 0.25 * diff(range(Egwf, na.rm = TRUE))
      ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
      text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F - G[fgi] * E[gwf])^2), n - 2))))
      p_calc <- vector(length = l)
      for (i in 1:l) {
         if (i == 1) {
            p_calc[i] <- p_[i]
         } else {
            fun <- function(x, pd, p_ini, Bg_ini, Btg_ini, Rs_ini, Rv_ini, Bw_ini, Gfgi, swi, p, Ef, Np, Gp, Wp, Wi, We, M, pvt) {
               Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear", ties = max)$y
               Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear", ties = max)$y
               Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear", ties = max)$y
               Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear", ties = max)$y
               Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear", ties = max)$y
               Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear", ties = max)$y
               dp <- p_ini - x
               Btg <- (Bg * (1 - Rs * Rv_ini) + Bo * (Rv_ini - Rv)) / (1 - Rs * Rv)
               Eg <- Btg - Btg_ini
               Ew <- Bw - Bw_ini
               Egwf <- Eg + Bg_ini * ((swi + M) * Ew / Bw_ini + (1 + M) * Ef_internal) / (1 - swi)
               F_lhs_calc <- Gp * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
               F_rhs_calc <- Gfgi * Egwf
               diff <- F_lhs_calc - F_rhs_calc
               return(diff)
            }
            p_calc[i] <- fzero(fun = fun, x = p_[i], pd = pd, p_ini = p_[1], Bg_ini = Bg[1], Btg_ini = Btg[1], Rs_ini = Rs[1], Rv_ini = Rv[1], Bw_ini = Bw[1], Gfgi = Gfgi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], We = We_[i], M = M, pvt = pvt, maxiter = 1000)$x
         }
      }
      STDERR <- sqrt(sum((p_ - p_calc) * ((p_ - p_calc)), na.rm = TRUE) / (l - 2))
      xmin <- min(c(p_,p_calc), na.rm = TRUE)
      xmax <- max(c(p_,p_calc), na.rm = TRUE)
      plot(p_, p_calc, pch = 21, bg = "blue" , xlab = "Pressure", ylab = "Pressure_estimated", xlim = c(xmin, xmax), ylim = c(xmin, xmax))
      abline(a = 0, b = 1)
      mtext(paste("slope = 1.0 " , "STDERR = ", round(STDERR,2)), side = 3)
      xrange <- xmax - xmin
      text(xmin + 0.25 * xrange, xmax - 0.25 * xrange, expression(STDERR == sqrt(frac(sum((p - p[est])^2), n - 2))))

      ymin <- max(0, min((F_[2:l]) / Egwf[2:l]) / 1.2)
      ymax <- 1.2 * max((F_[2:l]) / Egwf[2:l])
      plot(Gp, F_ / Egwf, pch = 21, bg = "blue" , xlab = "Gp", ylab = "F/Egwf", ylim = c(ymin, ymax))

      mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, G = G, phi = phi, swi = swi, pd = pd, p = p, cf = cf, M = M, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sgrw = sgrw, p_est = p_calc)
      class(mbal_lst) <- c("volumetric_gas", "mbal_gas")
      return(mbal_lst)
   }
}




# *************** Water Drive **************

water_drive_optim_gas <- function(optim_lst, time_lst) {

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
   input_unit <- optim_lst$input_unit
   output_unit <- optim_lst$output_unit
   G <- optim_lst$G
   phi <- optim_lst$phi
   swi <- optim_lst$swi
   pd <- optim_lst$pd
   p <- optim_lst$p
   l <- length(p)
   if (length(time_lst$t) != length(p)) stop("Lenght of 'time_lst$t' vector must be equal to the length of 'p' vector.")
   cf <- optim_lst$cf
   M <- optim_lst$M
   pvt <- optim_lst$pvt
   prod <- optim_lst$prod
   inj <- optim_lst$inj
   We <- optim_lst$We
   aquifer <- optim_lst$aquifer
   cls <- class(optim_lst)[1]
   aqu_cls <- class(aquifer)[1]
   if (aqu_cls == "We") We <- aquifer$We
   wf <- optim_lst$wf
   sgrw <- optim_lst$sgrw
   if (is.null(sgrw)) stop("'sorw' must be a numeric value.")
   keep <- which(wf == 1)
   p_ <- p[keep]
   time_ <- time[keep]
   if (length(cf) != 1) {
      cf_ <- cf[keep]
   } else {
      cf_ <- cf
   }
   prod_ <- prod[keep,]
   inj_ <- inj[keep,]
   if (aqu_cls == "We") We_ <- We[keep]
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
      Bo[i] <- approx(pvt$p, pvt$Bo, xout = p_[i])$y
      Rs[i] <- approx(pvt$p, pvt$Rs, xout = p_[i])$y
      muo[i] <- approx(pvt$p, pvt$muo, xout = p_[i])$y
      Rv[i] <- approx(pvt$p, pvt$Rv, xout = p_[i])$y
      Bg[i] <- approx(pvt$p, pvt$Bg, xout = p_[i])$y
      mug[i] <- approx(pvt$p, pvt$mug, xout = p_[i])$y
      Bw[i] <- approx(pvt$p, pvt$Bw, xout = p_[i])$y
      muw[i] <- approx(pvt$p, pvt$muw, xout = p_[i])$y
   }
   Np <- prod_$Np
   Gp <- prod_$Gp
   Wp <- prod_$Wp
   Wi <- inj_
   dp <- vector(length = l)
   Bto <- vector(length = l)
   Btg <- vector(length = l)
   Eo <- vector(length = l)
   Eg <- vector(length = l)
   Egwf <- vector(length = l)
   Ew <- vector(length = l)
   Ef <- vector(length = l)
   F_ <- vector(length = l)


   if (cls == "water_drive_M_We_optim_gas") {
      if (aqu_cls == "We") {
         # We <- We
      }
      if (aqu_cls == "uss_rad_edge") {
         We_ <- veh_uss_rad_edge(aquifer$phi, aquifer$perm_h, aquifer$h_a, aquifer$r_a, aquifer$r_R, aquifer$tetha, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "uss_rad_bottom") {
         We_ <- yk_uss_rad_bottom(aquifer$phi, aquifer$perm_h, aquifer$perm_v, aquifer$h_a, aquifer$r_a, aquifer$r_R, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "pss_rad_edge") {
         We_ <- fetkovich_pss_rad_edge(aquifer$phi, aquifer$perm_h, aquifer$h_a, aquifer$r_a, aquifer$r_R, aquifer$tetha, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "uss_lin_edge") {
         We_ <- nb_uss_lin_edge(aquifer$phi, aquifer$perm_h, aquifer$h_a, aquifer$w_a, aquifer$l_a, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "uss_lin_bottom") {
         We_ <- nb_uss_lin_bottom(aquifer$phi, aquifer$perm_v, aquifer$h_a, aquifer$w_a, aquifer$l_a, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "pss_lin_edge") {
         We_ <- fetk_pss_lin_edge(aquifer$phi, aquifer$perm_h, aquifer$h_a, aquifer$w_a, aquifer$l_a, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "pss_lin_bottom") {
         We_ <- fetk_pss_lin_bottom(aquifer$phi, aquifer$perm_v, aquifer$h_a, aquifer$w_a, aquifer$l_a, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "pot") {
         We_ <- pot(aquifer$phi, aquifer$h_a, aquifer$r_a, aquifer$r_R, aquifer$tetha, aquifer$c_water, aquifer$c_rock, aquifer$pressure[keep], aquifer$mult_len)
      }
      dp <- p_[1] - p_
      Bto <- (Bo * (1 - Rs[1] * Rv) + Bg * (Rs[1] - Rs)) / (1 - Rs * Rv)
      Btg <- (Bg * (1 - Rs * Rv[1]) + Bo * (Rv[1] - Rv)) / (1 - Rs * Rv)
      Eo <- Bto - Bto[1]
      Eg <- Btg - Btg[1]
      Ew <- Bw - Bw[1]
      for (i in 1:l) {
         Ef[i] <- -1 * trapz(p_[1:i], cf_[1:i])
      }
      Egwf <- Eg + Bg[1] * ((swi + M) * Ew / Bw[1] + (1 + M) * Ef) / (1 - swi)
      F_ <- Gp * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_
      data <- data.frame(Egwf = Egwf, `F_` = F_)
      mbal_lm <- lm(`F_` ~ 0 + Egwf, data = data)
      Gfgi <- mbal_lm$coefficients
      names(Gfgi) <- NULL
      Nfoi <- 0
      N <- Nfoi + Rv[1] * Gfgi
      G <- Gfgi + Rs[1] * Nfoi
      PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
      BV <- PV / phi
      W <- PV * swi / Bw[1]
      Et <- Gfgi * Egwf + We_ + (Wi - Wp) * Bw
      STDERR <- sqrt(sum((F_ - Gfgi * Egwf) * ((F_ - Gfgi * Egwf))) / (l - 2))
      plot(Egwf, F_, pch = 21, bg = "blue" , xlab = "Egwf", ylab = "F")
      abline(a = 0, b = Gfgi)
      mtext(paste("slope = Gfgi = ", round(Gfgi,0), " , STDERR = ", round(STDERR,2)), side = 3)
      xmin <- min(Egwf, na.rm = TRUE) + 0.25 * diff(range(Egwf, na.rm = TRUE))
      ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
      text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F - W[e] - G[fgi] * E[gwf])^2), n - 2))))
      p_calc <- vector(length = l)
      for (i in 1:l) {
         if (i == 1) {
            p_calc[i] <- p_[i]
         } else {
            fun <- function(x, pd, p_ini, Bg_ini, Btg_ini, Rs_ini, Rv_ini, Bw_ini, Gfgi, swi, p, Ef, Np, Gp, Wp, Wi, We, M, pvt) {
               Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear", ties = max)$y
               Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear", ties = max)$y
               Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear", ties = max)$y
               Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear", ties = max)$y
               Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear", ties = max)$y
               Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear", ties = max)$y
               dp <- p_ini - x
               Btg <- (Bg * (1 - Rs * Rv_ini) + Bo * (Rv_ini - Rv)) / (1 - Rs * Rv)
               Eg <- Btg - Btg_ini
               Ew <- Bw - Bw_ini
               Egwf <- Eg + Bg_ini * ((swi + M) * Ew / Bw_ini + (1 + M) * Ef_internal) / (1 - swi)
               F_lhs_calc <- Gp * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
               F_rhs_calc <- Gfgi * Egwf
               diff <- F_lhs_calc - F_rhs_calc
               return(diff)
            }
            p_calc[i] <- fzero(fun = fun, x = p_[i], pd = pd, p_ini = p_[1], Bg_ini = Bg[1], Btg_ini = Btg[1], Rs_ini = Rs[1], Rv_ini = Rv[1], Bw_ini = Bw[1], Gfgi = Gfgi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], We = We_[i], M = M, pvt = pvt, maxiter = 1000)$x
         }
      }
      STDERR <- sqrt(sum((p_ - p_calc) * ((p_ - p_calc)), na.rm = TRUE) / (l - 2))
      xmin <- min(c(p_,p_calc), na.rm = TRUE)
      xmax <- max(c(p_,p_calc), na.rm = TRUE)
      plot(p_, p_calc, pch = 21, bg = "blue" , xlab = "Pressure", ylab = "Pressure_estimated", xlim = c(xmin, xmax), ylim = c(xmin, xmax))
      abline(a = 0, b = 1)
      mtext(paste("slope = 1.0 " , "STDERR = ", round(STDERR,2)), side = 3)
      xrange <- xmax - xmin
      text(xmin + 0.25 * xrange, xmax - 0.25 * xrange, expression(STDERR == sqrt(frac(sum((p - p[est])^2), n - 2))))

      ymin <- max(0, min((F_[2:l] + We_[2:l]) / Egwf[2:l]) / 1.2)
      ymax <- 1.2 * max((F_[2:l] + We_[2:l]) / Egwf[2:l])
      plot(Gp, (F_ + We_) / Egwf, pch = 21, bg = "blue" , xlab = "Gp", ylab = "F/Egwf", ylim = c(ymin, ymax))

      mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, G = G, phi = phi, swi = swi, pd = pd, p = p, cf = cf, M = M, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sgrw = sgrw, p_est = p_calc)

      class(mbal_lst) <- c("water_drive_gas", "mbal_gas")
      return(mbal_lst)
   }



   if (cls == "water_drive_G_M_optim_gas") {
      dp <- p_[1] - p_
      Bto <- (Bo * (1 - Rs[1] * Rv) + Bg * (Rs[1] - Rs)) / (1 - Rs * Rv)
      Btg <- (Bg * (1 - Rs * Rv[1]) + Bo * (Rv[1] - Rv)) / (1 - Rs * Rv)
      Eo <- Bto - Bto[1]
      Eg <- Btg - Btg[1]
      Ew <- Bw - Bw[1]
      for (i in 1:l) {
         Ef[i] <- -1 * trapz(p_[1:i], cf_[1:i])
      }
      Egwf <- Eg + Bg[1] * ((swi + M) * Ew / Bw[1] + (1 + M) * Ef) / (1 - swi)
      Nfoi <- 0
      Gfgi <- G - Rs[1] * Nfoi
      F_ <- Gp * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw
      We_ <- F_ - Gfgi * Egwf
      if (aqu_cls %in% c("uss_rad_edge", "uss_rad_bottom", "pss_rad_edge")) {
         phi <- aquifer$phi
         perm_h <- aquifer$perm_h
         if (aqu_cls == "uss_rad_bottom") {
            perm_v <- aquifer$perm_v
         }
         h_a <- aquifer$h_a
         r_a <- aquifer$r_a
         r_R <- aquifer$r_R
         if (aqu_cls %in% c("uss_rad_edge", "pss_rad_edge")) {
            tetha <- aquifer$tetha
         }
         mu_water <- aquifer$mu_water
         c_water <- aquifer$c_water
         c_rock <- aquifer$c_rock
         pressure_ <- aquifer$pressure[keep]
         mult_len <- aquifer$mult_len
         if (is.null(mult_len)) {
            if (aqu_cls %in% c("uss_rad_edge", "pss_rad_edge")) {
               mult_len <- c(1, 1)
            }
            if (aqu_cls == "uss_rad_bottom") {
               mult_len <- c(1, 1, 1)
            }
         }
      }
      if (aqu_cls %in% c("uss_lin_edge", "uss_lin_bottom", "pss_lin_edge", "pss_lin_bottom")) {
         phi <- aquifer$phi
         if (aqu_cls %in% c("uss_lin_edge", "pss_lin_edge")) {
            perm_h <- aquifer$perm_h
         }
         if (aqu_cls %in% c("uss_lin_bottom", "pss_lin_bottom")) {
            perm_v <- aquifer$perm_v
         }
         h_a <- aquifer$h_a
         w_a <- aquifer$w_a
         l_a <- aquifer$l_a
         mu_water <- aquifer$mu_water
         c_water <- aquifer$c_water
         c_rock <- aquifer$c_rock
         pressure_ <- aquifer$pressure[keep]
         mult_len <- aquifer$mult_len
         if (is.null(mult_len)) {
            mult_len <- c(1, 1)
         }
      }
      if (aqu_cls == "pot") {
         phi <- aquifer$phi
         h_a <- aquifer$h_a
         r_a <- aquifer$r_a
         r_R <- aquifer$r_R
         tetha <- aquifer$tetha
         c_water <- aquifer$c_water
         c_rock <- aquifer$c_rock
         pressure_ <- aquifer$pressure[keep]
         mult_len <- aquifer$mult_len
         if (is.null(mult_len)) {
            mult_len <- c(1)
         }
      }
      if (aqu_cls == "uss_rad_edge") {
         error_fun_ussre_gas <- function(par, We, phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- veh_uss_rad_edge(phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_ussre_gas, We = We_, phi = phi, perm_h = perm_h, h_a = h_a, r_a = r_a, r_R = r_R,
                           tetha = tetha, mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_,
                           pressure = pressure_, lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "uss_rad_bottom") {
         error_fun_ussrb_gas <- function(par, We, phi, perm_h, perm_v, h_a, r_a, r_R, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- yk_uss_rad_bottom(phi, perm_h, perm_v, h_a, r_a, r_R, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_ussrb_gas, We = We_, phi = phi, perm_h = perm_h, perm_v = perm_v, h_a = h_a, r_a = r_a,
                           r_R = r_R, mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_, pressure = pressure_,
                           lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "pss_rad_edge") {
         error_fun_pssre_gas <- function(par, We, phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- fetkovich_pss_rad_edge(phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_pssre_gas, We = We_, phi = phi, perm_h = perm_h, h_a = h_a, r_a = r_a,
                           r_R = r_R, tetha = tetha, mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_,
                           pressure = pressure_, lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "uss_lin_edge") {
         error_fun_ussle_gas <- function(par, We, phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- nb_uss_lin_edge(phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_ussle_gas, We = We_, phi = phi, perm_h = perm_h, h_a = h_a, w_a = w_a, l_a = l_a,
                           mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_, pressure = pressure_,
                           lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "uss_lin_bottom") {
         error_fun_usslb_gas <- function(par, We, phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- nb_uss_lin_bottom(phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_usslb_gas, We = We_, phi = phi, perm_v = perm_v, h_a = h_a, w_a = w_a, l_a = l_a,
                           mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_, pressure = pressure_,
                           lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "pss_lin_edge") {
         error_fun_pssle_gas <- function(par, We, phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- fetk_pss_lin_edge(phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_pssle_gas, We = We_, phi = phi, perm_h = perm_h, h_a = h_a, w_a = w_a, l_a = l_a,
                           mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_, pressure = pressure_,
                           lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "pss_lin_bottom") {
         error_fun_psslb_gas <- function(par, We, phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- fetk_pss_lin_bottom(phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_psslb_gas, We = We_, phi = phi, perm_v = perm_v, h_a = h_a, w_a = w_a, l_a = l_a,
                           mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_, pressure = pressure_,
                           lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "pot") {
         error_fun_pot_gas <- function(par, We, phi, h_a, r_a, r_R, tetha, c_water, c_rock, time, pressure) {
            We_calc <- pot(phi, h_a, r_a, r_R, tetha, c_water, c_rock, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_pot_gas, We = We_, phi = phi, h_a = h_a, r_a = r_a, r_R = r_R,
                           tetha = tetha, c_water = c_water, c_rock = c_rock, time = time_,
                           pressure = pressure_, lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }

      if (aqu_cls == "uss_rad_edge") {
         aquifer$r_a <- aquifer$r_a * mult_len_est[1]
         r_a <- r_a * mult_len_est[1]
         aquifer$perm_h <- aquifer$perm_h * mult_len_est[2]
         perm_h <- perm_h * mult_len_est[2]
         We_calc <- veh_uss_rad_edge(phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time_, pressure_, c(1,1))
      }
      if (aqu_cls == "uss_rad_bottom") {
         aquifer$r_a <- aquifer$r_a * mult_len_est[1]
         r_a <- r_a * mult_len_est[1]
         aquifer$perm_h <- aquifer$perm_h * mult_len_est[2]
         perm_h <- perm_h * mult_len_est[2]
         aquifer$perm_v <- aquifer$perm_v * mult_len_est[3]
         perm_h <- perm_v * mult_len_est[3]
         We_calc <- yk_uss_rad_bottom(phi, perm_h, perm_v, h_a, r_a, r_R, mu_water, c_water, c_rock, time_, pressure_, c(1,1,1))
      }
      if (aqu_cls == "pss_rad_edge") {
         aquifer$r_a <- aquifer$r_a * mult_len_est[1]
         r_a <- r_a * mult_len_est[1]
         aquifer$perm_h <- aquifer$perm_h * mult_len_est[2]
         perm_h <- perm_h * mult_len_est[2]
         We_calc <- fetkovich_pss_rad_edge(phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time_, pressure_, c(1,1))
      }
      if (aqu_cls == "uss_lin_edge") {
         aquifer$l_a <- aquifer$l_a * mult_len_est[1]
         l_a <- l_a * mult_len_est[1]
         aquifer$perm_h <- aquifer$perm_h * mult_len_est[2]
         perm_h <- perm_h * mult_len_est[2]
         We_calc <- nb_uss_lin_edge(phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time_, pressure_, c(1,1))
      }
      if (aqu_cls == "uss_lin_bottom") {
         aquifer$h_a <- aquifer$h_a * mult_len_est[1]
         h_a <- h_a * mult_len_est[1]
         aquifer$perm_v <- aquifer$perm_v * mult_len_est[2]
         perm_v <- perm_v * mult_len_est[2]
         We_calc <- nb_uss_lin_bottom(phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time_, pressure_, c(1,1))
      }
      if (aqu_cls == "pss_lin_edge") {
         aquifer$l_a <- aquifer$l_a * mult_len_est[1]
         l_a <- l_a * mult_len_est[1]
         aquifer$perm_h <- aquifer$perm_h * mult_len_est[2]
         perm_h <- perm_h * mult_len_est[2]
         We_calc <- fetk_pss_lin_edge(phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time_, pressure_, c(1,1))
      }
      if (aqu_cls == "pss_lin_bottom") {
         aquifer$h_a <- aquifer$h_a * mult_len_est[1]
         h_a <- h_a * mult_len_est[1]
         aquifer$perm_v <- aquifer$perm_v * mult_len_est[2]
         perm_v <- perm_v * mult_len_est[2]
         We_calc <- fetk_pss_lin_bottom(phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time_, pressure_, c(1,1))
      }
      if (aqu_cls == "pot") {
         aquifer$r_a <- aquifer$r_a * mult_len_est[1]
         r_a <- r_a * mult_len_est[1]
         We_calc <- pot(phi, h_a, r_a, r_R, tetha, c_water, c_rock, pressure_, c(1))
      }
      if (aqu_cls == "pot") {
         aquifer$mult_len <- c(1)
      } else if (aqu_cls == "uss_rad_bottom") {
         aquifer$mult_len <- c(1,1,1)
      } else {
         aquifer$mult_len <- c(1,1)
      }
      aquifer$lower <- NULL
      aquifer$upper <- NULL
      aquifer$control <- NULL
      F_ <- Gp * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_calc
      N <- Nfoi + Rv[1] * Gfgi
      G <- Gfgi + Rs[1] * Nfoi
      PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
      BV <- PV / phi
      W <- PV * swi / Bw[1]
      Et <- Gfgi * Egwf + We_calc + (Wi - Wp) * Bw
      STDERR <- sqrt(sum((F_ - Gfgi * Egwf) * ((F_ - Gfgi * Egwf))) / (l - 2))
      plot(Egwf, F_, pch = 21, bg = "blue" , xlab = "Egwf", ylab = "F")
      abline(a = 0, b = Gfgi)
      mtext(paste("slope = Gfgi = ", round(Gfgi,0), " , STDERR = ", round(STDERR,2)), side = 3)
      xmin <- min(Egwf, na.rm = TRUE) + 0.25 * diff(range(Egwf, na.rm = TRUE))
      ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
      text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F - W[e] - G[fgi] * E[gwf])^2), n - 2))))
      p_calc <- vector(length = l)
      for (i in 1:l) {
         if (i == 1) {
            p_calc[i] <- p_[i]
         } else {
            fun <- function(x, pd, p_ini, Bg_ini, Btg_ini, Rs_ini, Rv_ini, Bw_ini, Gfgi, swi, p, Ef, Np, Gp, Wp, Wi, We, M, pvt) {
               Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear", ties = max)$y
               Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear", ties = max)$y
               Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear", ties = max)$y
               Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear", ties = max)$y
               Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear", ties = max)$y
               Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear", ties = max)$y
               dp <- p_ini - x
               Btg <- (Bg * (1 - Rs * Rv_ini) + Bo * (Rv_ini - Rv)) / (1 - Rs * Rv)
               Eg <- Btg - Btg_ini
               Ew <- Bw - Bw_ini
               Egwf <- Eg + Bg_ini * ((swi + M) * Ew / Bw_ini + (1 + M) * Ef_internal) / (1 - swi)
               F_lhs_calc <- Gp * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
               F_rhs_calc <- Gfgi * Egwf
               diff <- F_lhs_calc - F_rhs_calc
               return(diff)
            }
            p_calc[i] <- fzero(fun = fun, x = p_[i], pd = pd, p_ini = p_[1], Bg_ini = Bg[1], Btg_ini = Btg[1], Rs_ini = Rs[1], Rv_ini = Rv[1], Bw_ini = Bw[1], Gfgi = Gfgi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], We = We_calc[i], M = M, pvt = pvt, maxiter = 1000)$x
         }
      }
      STDERR <- sqrt(sum((p_ - p_calc) * ((p_ - p_calc)), na.rm = TRUE) / (l - 2))
      xmin <- min(c(p_,p_calc), na.rm = TRUE)
      xmax <- max(c(p_,p_calc), na.rm = TRUE)
      plot(p_, p_calc, pch = 21, bg = "blue" , xlab = "Pressure", ylab = "Pressure_estimated", xlim = c(xmin, xmax), ylim = c(xmin, xmax))
      abline(a = 0, b = 1)
      mtext(paste("slope = 1.0 " , "STDERR = ", round(STDERR,2)), side = 3)
      xrange <- xmax - xmin
      text(xmin + 0.25 * xrange, xmax - 0.25 * xrange, expression(STDERR == sqrt(frac(sum((p - p[est])^2), n - 2))))

      ymin <- max(0, min((F_[2:l] + We_calc[2:l]) / Egwf[2:l]) / 1.2)
      ymax <- 1.2 * max((F_[2:l] + We_calc[2:l]) / Egwf[2:l])
      plot(Gp, (F_ + We_calc)/ Egwf, pch = 21, bg = "blue" , xlab = "Gp", ylab = "F/Egwf", ylim = c(ymin, ymax))

      mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, G = G, phi = phi, swi = swi, pd = pd, p = p, cf = cf, M = M, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sgrw = sgrw, p_est = p_calc)

      class(mbal_lst) <- c("water_drive_gas", "mbal_gas")
      return(mbal_lst)
   }

   if (cls == "water_drive_We_optim_gas") {
      if (aqu_cls == "We") {
         # We <- We
      }
      if (aqu_cls == "uss_rad_edge") {
         We_ <- veh_uss_rad_edge(aquifer$phi, aquifer$perm_h, aquifer$h_a, aquifer$r_a, aquifer$r_R, aquifer$tetha, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "uss_rad_bottom") {
         We_ <- yk_uss_rad_bottom(aquifer$phi, aquifer$perm_h, aquifer$perm_v, aquifer$h_a, aquifer$r_a, aquifer$r_R, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "pss_rad_edge") {
         We_ <- fetkovich_pss_rad_edge(aquifer$phi, aquifer$perm_h, aquifer$h_a, aquifer$r_a, aquifer$r_R, aquifer$tetha, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "uss_lin_edge") {
         We_ <- nb_uss_lin_edge(aquifer$phi, aquifer$perm_h, aquifer$h_a, aquifer$w_a, aquifer$l_a, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "uss_lin_bottom") {
         We_ <- nb_uss_lin_bottom(aquifer$phi, aquifer$perm_v, aquifer$h_a, aquifer$w_a, aquifer$l_a, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "pss_lin_edge") {
         We_ <- fetk_pss_lin_edge(aquifer$phi, aquifer$perm_h, aquifer$h_a, aquifer$w_a, aquifer$l_a, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "pss_lin_bottom") {
         We_ <- fetk_pss_lin_bottom(aquifer$phi, aquifer$perm_v, aquifer$h_a, aquifer$w_a, aquifer$l_a, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "pot") {
         We_ <- pot(aquifer$phi, aquifer$h_a, aquifer$r_a, aquifer$r_R, aquifer$tetha, aquifer$c_water, aquifer$c_rock, aquifer$pressure[keep], aquifer$mult_len)
      }
      dp <- p_[1] - p_
      Bto <- (Bo * (1 - Rs[1] * Rv) + Bg * (Rs[1] - Rs)) / (1 - Rs * Rv)
      Btg <- (Bg * (1 - Rs * Rv[1]) + Bo * (Rv[1] - Rv)) / (1 - Rs * Rv)
      Eo <- Bto - Bto[1]
      Eg <- Btg - Btg[1]
      Ew <- Bw - Bw[1]
      for (i in 1:l) {
         Ef[i] <- -1 * trapz(p_[1:i], cf_[1:i])
      }
      A <- Eg + Bg[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      B <- Bg[1] * (Ew / Bw[1] + Ef) / (1 - swi)
      F_ <- Gp * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_
      data <- data.frame(LHS = F_ / A, RHS = B / A)
      mbal_lm <- lm(LHS ~ RHS, data = data)
      Gfgi <- mbal_lm$coefficients[1]
      names(Gfgi) <- NULL
      M <- mbal_lm$coefficients[2] / Gfgi
      names(M) <- NULL
      Egwf <- Eg + Bg[1] * ((swi + M) * Ew / Bw[1] + (1 + M) * Ef) / (1 - swi)
      Nfoi <- 0
      N <- Nfoi + Rv[1] * Gfgi
      G <- Gfgi + Rs[1] * Nfoi
      PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
      BV <- PV / phi
      W <- PV * swi / Bw[1]
      Et <- Gfgi * Egwf + We_ + (Wi - Wp) * Bw
      STDERR <- sqrt(sum((F_ - Gfgi * Egwf) * ((F_ - Gfgi * Egwf))) / (l - 2))
      plot(Egwf, F_, pch = 21, bg = "blue" , xlab = "Egwf", ylab = "F")
      abline(a = 0, b = Gfgi)
      mtext(paste("slope = Gfgi = ", round(Gfgi,0), ", M = ", round(M,3), " , STDERR = ", round(STDERR,2)), side = 3)
      xmin <- min(Egwf, na.rm = TRUE) + 0.25 * diff(range(Egwf, na.rm = TRUE))
      ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
      text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F - G[fgi] * E[gwf])^2), n - 2))))
      p_calc <- vector(length = l)
      for (i in 1:l) {
         if (i == 1) {
            p_calc[i] <- p_[i]
         } else {
            fun <- function(x, pd, p_ini, Bg_ini, Btg_ini, Rs_ini, Rv_ini, Bw_ini, Gfgi, swi, p, Ef, Np, Gp, Wp, Wi, We, M, pvt) {
               Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear", ties = max)$y
               Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear", ties = max)$y
               Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear", ties = max)$y
               Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear", ties = max)$y
               Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear", ties = max)$y
               Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear", ties = max)$y
               dp <- p_ini - x
               Btg <- (Bg * (1 - Rs * Rv_ini) + Bo * (Rv_ini - Rv)) / (1 - Rs * Rv)
               Eg <- Btg - Btg_ini
               Ew <- Bw - Bw_ini
               Egwf <- Eg + Bg_ini * ((swi + M) * Ew / Bw_ini + (1 + M) * Ef_internal) / (1 - swi)
               F_lhs_calc <- Gp * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
               F_rhs_calc <- Gfgi * Egwf
               diff <- F_lhs_calc - F_rhs_calc
               return(diff)
            }
            p_calc[i] <- fzero(fun = fun, x = p_[i], pd = pd, p_ini = p_[1], Bg_ini = Bg[1], Btg_ini = Btg[1], Rs_ini = Rs[1], Rv_ini = Rv[1], Bw_ini = Bw[1], Gfgi = Gfgi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], We = We_[i], M = M, pvt = pvt, maxiter = 1000)$x
         }
      }
      STDERR <- sqrt(sum((p_ - p_calc) * ((p_ - p_calc)), na.rm = TRUE) / (l - 2))
      xmin <- min(c(p_,p_calc), na.rm = TRUE)
      xmax <- max(c(p_,p_calc), na.rm = TRUE)
      plot(p_, p_calc, pch = 21, bg = "blue" , xlab = "Pressure", ylab = "Pressure_estimated", xlim = c(xmin, xmax), ylim = c(xmin, xmax))
      abline(a = 0, b = 1)
      mtext(paste("slope = 1.0 " , "STDERR = ", round(STDERR,2)), side = 3)
      xrange <- xmax - xmin
      text(xmin + 0.25 * xrange, xmax - 0.25 * xrange, expression(STDERR == sqrt(frac(sum((p - p[est])^2), n - 2))))

      ymin <- max(0, min((F_[2:l] + We_[2:l]) / Egwf[2:l]) / 1.2)
      ymax <- 1.2 * max((F_[2:l] + We_[2:l]) / Egwf[2:l])
      plot(Gp, (F_ + We_)/ Egwf, pch = 21, bg = "blue" , xlab = "Gp", ylab = "F/Egwf", ylim = c(ymin, ymax))

      mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, G = G, phi = phi, swi = swi, pd = pd, p = p, cf = cf, M = M, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sgrw = sgrw, p_est = p_calc)

      class(mbal_lst) <- c("water_drive_gas", "mbal_gas")
      return(mbal_lst)
   }
   if (cls == "water_drive_G_We_optim_gas") {
      if (aqu_cls == "We") {
         # We <- We
      }
      if (aqu_cls == "uss_rad_edge") {
         We_ <- veh_uss_rad_edge(aquifer$phi, aquifer$perm_h, aquifer$h_a, aquifer$r_a, aquifer$r_R, aquifer$tetha, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "uss_rad_bottom") {
         We_ <- yk_uss_rad_bottom(aquifer$phi, aquifer$perm_h, aquifer$perm_v, aquifer$h_a, aquifer$r_a, aquifer$r_R, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "pss_rad_edge") {
         We_ <- fetkovich_pss_rad_edge(aquifer$phi, aquifer$perm_h, aquifer$h_a, aquifer$r_a, aquifer$r_R, aquifer$tetha, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "uss_lin_edge") {
         We_ <- nb_uss_lin_edge(aquifer$phi, aquifer$perm_h, aquifer$h_a, aquifer$w_a, aquifer$l_a, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "uss_lin_bottom") {
         We_ <- nb_uss_lin_bottom(aquifer$phi, aquifer$perm_v, aquifer$h_a, aquifer$w_a, aquifer$l_a, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "pss_lin_edge") {
         We_ <- fetk_pss_lin_edge(aquifer$phi, aquifer$perm_h, aquifer$h_a, aquifer$w_a, aquifer$l_a, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "pss_lin_bottom") {
         We_ <- fetk_pss_lin_bottom(aquifer$phi, aquifer$perm_v, aquifer$h_a, aquifer$w_a, aquifer$l_a, aquifer$mu_water, aquifer$c_water, aquifer$c_rock, time_, aquifer$pressure[keep], aquifer$mult_len)
      }
      if (aqu_cls == "pot") {
         We_ <- pot(aquifer$phi, aquifer$h_a, aquifer$r_a, aquifer$r_R, aquifer$tetha, aquifer$c_water, aquifer$c_rock, aquifer$pressure[keep], aquifer$mult_len)
      }
      dp <- p_[1] - p_
      Bto <- (Bo * (1 - Rs[1] * Rv) + Bg * (Rs[1] - Rs)) / (1 - Rs * Rv)
      Btg <- (Bg * (1 - Rs * Rv[1]) + Bo * (Rv[1] - Rv)) / (1 - Rs * Rv)
      Eo <- Bto - Bto[1]
      Eg <- Btg - Btg[1]
      Ew <- Bw - Bw[1]
      for (i in 1:l) {
         Ef[i] <- -1 * trapz(p_[1:i], cf_[1:i])
      }
      Gfgi <- G
      A <- Eg + Bg[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      B <- Bg[1] * (Ew / Bw[1] + Ef) / (1 - swi)
      F_ <- Gp * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_
      data <- data.frame(LHS = F_ / Gfgi - A, RHS = B)
      mbal_lm <- lm(LHS ~ 0 + RHS, data = data)
      M <- mbal_lm$coefficients
      names(M) <- NULL
      Egwf <- Eg + Bg[1] * ((swi + M) * Ew / Bw[1] + (1 + M) * Ef) / (1 - swi)
      Nfoi <- 0
      N <- Nfoi + Rv[1] * Gfgi
      PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
      BV <- PV / phi
      W <- PV * swi / Bw[1]
      Et <- Gfgi * Egwf + We_ + (Wi - Wp) * Bw
      STDERR <- sqrt(sum((F_ - Gfgi * Egwf) * ((F_ - Gfgi * Egwf))) / (l - 2))
      plot(Egwf, F_, pch = 21, bg = "blue" , xlab = "Egwf", ylab = "F")
      abline(a = 0, b = Gfgi)
      mtext(paste("slope = Gfgi = ", round(Gfgi,0), ", M = ", round(M,3), " , STDERR = ", round(STDERR,2)), side = 3)
      xmin <- min(Egwf, na.rm = TRUE) + 0.25 * diff(range(Egwf, na.rm = TRUE))
      ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
      text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F - G[fgi] * E[gwf])^2), n - 2))))
      p_calc <- vector(length = l)
      for (i in 1:l) {
         if (i == 1) {
            p_calc[i] <- p_[i]
         } else {
            fun <- function(x, pd, p_ini, Bg_ini, Btg_ini, Rs_ini, Rv_ini, Bw_ini, Gfgi, swi, p, Ef, Np, Gp, Wp, Wi, We, M, pvt) {
               Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear", ties = max)$y
               Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear", ties = max)$y
               Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear", ties = max)$y
               Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear", ties = max)$y
               Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear", ties = max)$y
               Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear", ties = max)$y
               dp <- p_ini - x
               Btg <- (Bg * (1 - Rs * Rv_ini) + Bo * (Rv_ini - Rv)) / (1 - Rs * Rv)
               Eg <- Btg - Btg_ini
               Ew <- Bw - Bw_ini
               Egwf <- Eg + Bg_ini * ((swi + M) * Ew / Bw_ini + (1 + M) * Ef_internal) / (1 - swi)
               F_lhs_calc <- Gp * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
               F_rhs_calc <- Gfgi * Egwf
               diff <- F_lhs_calc - F_rhs_calc
               return(diff)
            }
            p_calc[i] <- fzero(fun = fun, x = p_[i], pd = pd, p_ini = p_[1], Bg_ini = Bg[1], Btg_ini = Btg[1], Rs_ini = Rs[1], Rv_ini = Rv[1], Bw_ini = Bw[1], Gfgi = Gfgi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], We = We_[i], M = M, pvt = pvt, maxiter = 1000)$x
         }
      }
      STDERR <- sqrt(sum((p_ - p_calc) * ((p_ - p_calc)), na.rm = TRUE) / (l - 2))
      xmin <- min(c(p_,p_calc), na.rm = TRUE)
      xmax <- max(c(p_,p_calc), na.rm = TRUE)
      plot(p_, p_calc, pch = 21, bg = "blue" , xlab = "Pressure", ylab = "Pressure_estimated", xlim = c(xmin, xmax), ylim = c(xmin, xmax))
      abline(a = 0, b = 1)
      mtext(paste("slope = 1.0 " , "STDERR = ", round(STDERR,2)), side = 3)
      xrange <- xmax - xmin
      text(xmin + 0.25 * xrange, xmax - 0.25 * xrange, expression(STDERR == sqrt(frac(sum((p - p[est])^2), n - 2))))
      mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, G = G, phi = phi, swi = swi, pd = pd, p = p, cf = cf, M = M, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sgrw = sgrw, p_est = p_calc)

      ymin <- max(0, min((F_[2:l] + We_[2:l]) / Egwf[2:l]) / 1.2)
      ymax <- 1.2 * max((F_[2:l] + We_[2:l]) / Egwf[2:l])
      plot(Gp, (F_ + We_)/ Egwf, pch = 21, bg = "blue" , xlab = "Gp", ylab = "F/Egwf", ylim = c(ymin, ymax))

      class(mbal_lst) <- c("water_drive_gas", "mbal_gas")
      return(mbal_lst)
   }
}




