
#' A list object of class 'optimization_oil' for material balance analysis
#'
#' Create an object of class 'optimization_oil'
#'
#' @param input_unit a unit system for parameters, only the character string 'Field' is accepted
#' @param output_unit a unit system for properties, only the character string 'Field' is accepted
#' @param unknown_param a character string showing the unknown parameter(s). One of the following options: 'N', 'm', 'We', or 'N_m'
#' @param aquifer_model defaulted to `NULL`, otherwise must be a character string, one of the following eight options: 'uss_rad_edge', 'uss_rad_bottom', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom', 'pot'. For further information about each model, please see 'Raquifer' package reference manual (https://cran.r-project.org/web/packages/Raquifer/index.html)
#' @param N original oil in place, STB. If unknown, a `NULL` value must be assigned
#' @param m ratio of original gas cap volume to original oil leg volume, a numeric. If unknown, a `NULL` value must be assigned
#' @param phi reservoir porosity, a numeric fraction
#' @param swi initial water saturation in the reservoir, a numeric fraction
#' @param Np cumulative oil production, STB
#' @param Rp ratio of cumulative produced gas to cumulative produced oil
#' @param Wp cumulative water production, STB
#' @param Gi cumulative gas injection, SCF
#' @param Wi cumulative water injection, STB
#' @param We cumulative aquifer water influx, BBL. If unknown, a `NULL` value must be assigned
#' @param pb bubble point pressure, a numeric value, psi
#' @param p reservoir pressure, a numeric vector, psi
#' @param pvt a data frame of PVT properties including pressure 'p' in 'psi', oil formation volume factor 'Bo' in 'bbl/stb', solution gas-oil ratio 'Rs' in 'scf/stb', oil viscosity 'muo' in 'cp', volatilized oil-gas ratio 'Rv' in 'stb/scf', gas formation volume factor 'Bg' in 'bbl/scf', gas viscosity 'mug' in 'cp', water formation volume factor 'Bw' in 'bbl/stb', and water viscosity 'muw' in 'cp'
#' @param cf formation compressibility, a numeric value or vector, 1/psi
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
#' @param sorg residual oil saturation in gas invaded zone (gas cap expansion or gas injection), a numeric fraction
#' @param sorw residual oil saturation in water invaded zone (aquifer encroachment or water injection), a numeric fraction
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
#' @return a list of class 'mbal_oil' with all the required parameters for the mbal_perform_oil() S3 methods
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
#' p <- c(3330, 3150, 3000, 2850, 2700, 2550, 2400)
#'
#' We <- rep(0, length.out = length(p))
#'
#' Np <- c(0, 3.295, 5.903, 8.852, 11.503, 14.513, 17.730) * 1e6
#'
#' Rp <- c(0, 1050, 1060, 1160, 1235, 1265, 1300)
#'
#' Wp <- rep(0, length.out = length(p))
#'
#' Wi <- rep(0, length.out = length(p))
#'
#' Gi <- rep(0, length.out = length(p))
#'
#' wf <- c(1, 1, 1, 0, 1, 0, 1)
#'
#' mbal_optim_oil_lst <- mbal_optim_param_oil(input_unit = "Field",
#' output_unit = "Field", unknown_param = "N_m", aquifer_model = NULL,
#' phi = 0.2, swi = 0.2, Np = Np,Rp = Rp, Wp = Wp, Gi = Gi, Wi = Wi,
#' We = We, pb = 3330, p = p, pvt = pvt_table, cf = 0, wf = wf,
#' sorg = 0.2, sorw = 0)
#'
#' dplyr::glimpse(mbal_optim_oil_lst)

mbal_optim_param_oil <- function(input_unit = "Field", output_unit = "Field", unknown_param = NULL, aquifer_model = NULL, N = NULL, m = NULL, phi = NULL, swi = NULL, Np = NULL, Rp = NULL, Wp = NULL, Gi = NULL, Wi = NULL, We = NULL, pb = NULL, p = NULL, pvt = NULL, cf = NULL, phi_a = NULL, perm_h_a = NULL, perm_v_a = NULL, h_a = NULL, r_a = NULL, r_R = NULL, w_a = NULL, l_a = NULL, tetha = NULL, muw_a = NULL, cw_a = NULL, cf_a = NULL, wf = NULL, sorg = NULL, sorw = NULL, mult_len = NULL, lower = NULL, upper = NULL, control = NULL) {

   if (!is.character(input_unit)) stop("'input_unit' must be the character string 'Field'.")
   if (input_unit != "Field") stop("'input_unit' must be the character string 'Field'.")
   if (!is.character(output_unit)) stop("'output_unit' must be the character string 'Field'.")
   if (output_unit != "Field") stop("'output_unit' must be the character string 'Field'.")
   if (!is.character(unknown_param))  stop("'unknown_param' must be one of the character strings from the following combinations: 'N', 'm', 'We', or 'N_m'.")
   if (is.null(unknown_param)) stop("One of the following combinations must be selected: 'N', 'm', 'We', or 'N_m'.")
   if (!(unknown_param %in% c('N', 'm', 'We', 'N_m'))) stop("One of the following combinations must be selected: 'N', 'm', 'We', or 'N_m'.")

   if (unknown_param == "N") {
      if (!is.null(N)) stop("'N' must be a numeric for the 'unknown_param = N' case.")
      if (is.null(m)) stop("'m' must be a numeric for the 'unknown_param = N' case.")
      if (!is.numeric(m)) stop("'m' must be a numeric for the 'unknown_param = N' case.")
      if (m < 0) stop("'m' must be a numeric equal to or greater than zero for the 'unknown_param = N' case.")
      if (is.null(aquifer_model)) {
         if (is.null(We)) stop("Either 'We' or 'aquifer_model' must be known in advance for the unknown_param = 'N' case.")
      }
      if ((!is.null(aquifer_model)) & (!is.null(We))) {
         stop("Either 'We' or 'aquifer_model' must be NULL for the unknown_param = 'N' case.")
      }
   }
   if (unknown_param == "m") {
      if (!is.null(m)) stop("'m' must be a numeric for the 'unknown_param = m' case.")
      if (is.null(N)) stop("'N' must be a numeric for the 'unknown_param = m' case.")
      if (!is.numeric(N)) stop("'N' must be a numeric for the 'unknown_param = m' case.")
      if (N < 0) stop("'N' must be a numeric equal to or greater than zero for the 'unknown_param = m' case.")
      if (is.null(aquifer_model)) {
         if (is.null(We)) stop("Either 'We' or 'aquifer_model' must be known in advance for the unknown_param = 'm' case.")
      }
      if ((!is.null(aquifer_model)) & (!is.null(We))) {
         stop("Either 'We' or 'aquifer_model' must be NULL for the unknown_param = 'm' case.")
      }
   }
   if (unknown_param == "We") {
      if (is.null(aquifer_model)) stop("'aquifer_model' must be known for the 'unknown_param = We' case.")
      if (!is.null(We)) stop("'We' must be NULL for the 'unknown_param = We' case.")
      if (is.null(N)) stop("'N' must be a numeric for the 'unknown_param = We' case.")
      if (!is.numeric(N)) stop("'N' must be a numeric for the 'unknown_param = We' case.")
      if (N < 0) stop("'N' must be a numeric equal to or greater than zero for the 'unknown_param = We' case.")
      if (is.null(m)) stop("'m' must be a numeric for the 'unknown_param = We' case.")
      if (!is.numeric(m)) stop("'m' must be a numeric for the 'unknown_param = We' case.")
      if (m < 0) stop("'m' must be a numeric equal to or greater than zero for the 'unknown_param = We' case.")
   }
   if (unknown_param == "N_m") {
      if (!is.null(N)) stop("'N' must be a numeric for the 'unknown_param = N_m' case.")
      if (!is.null(m)) stop("'m' must be a numeric for the 'unknown_param = N_m' case.")
      if (is.null(aquifer_model)) {
         if (is.null(We)) stop("Either 'We' or 'aquifer_model' must be known in advance for the unknown_param = 'N_m' case.")
      }
      if ((!is.null(aquifer_model)) & (!is.null(We))) {
         stop("Either 'We' or 'aquifer_model' must be NULL for the unknown_param = 'N_m' case.")
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
   if (is.null(Rp)) stop("'Rp' must be a numeric vector.")
   if (!is.numeric(Rp)) stop("'Rp' must be a numeric vector.")
   if (Rp[1] != 0) stop("First reported 'Rp' value must be zero.")
   if (is.null(pb)) stop("'pb' must be a numeric value.")
   if (!is.numeric(pb)) stop("'pb' must be a numeric value.")
   if (length(pb) != 1) stop("'pb' must be a numeric value.")
   if (is.null(p)) stop("'p' must be a numeric vector.")
   if (!is.numeric(p)) stop("'p' must be a numeric vector.")
   if (p[1] < pb) stop("Initial reservoir pressure must be equal to or greater than 'pb'.")
   l <- length(p)
   if (length(Np) != l) stop("'p' and 'Np' vectors must have the same length.")
   if (length(Rp) != l) stop("'p' and 'Rp' vectors must have the same length.")
   if (is.null(pvt)) stop("'pvt' must be a data frame with columns 'p', 'Bo', 'Rs', 'muo', 'Rv', 'Bg', 'mug', 'Bw', and 'muw'.")
   if (!is.data.frame(pvt)) stop("'pvt' must be a data frame with columns 'p', 'Bo', 'Rs', 'muo', 'Rv', 'Bg', 'mug', 'Bw', and 'muw'.")
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
      Rs_pvt_max <- pvt[pvt$p == max(pvt$p),]$Rs
      for (i in 1:l) {
         if (i == 1) {
            Rp[i] <- Rs_pvt_max
         } else {
            if (p[i] >= pb) {
               diff <- abs(Rs_pvt_max - Rp[i])
               if (diff > 1e-9) stop("'Rp' values are not equal to the 'Rs' values in the 'pvt' data frame at pressures equal to or greater than 'pb'.")
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
   if (is.null(Gi)) {
      Gi <- rep(0, l)
   }
   if (!is.null(Gi)) {
      if (!is.numeric(Gi)) stop("'Gi' must be a numeric vector.")
      if (is.numeric(Gi)) {
         if (length(Gi) != l) stop("'p' and 'Gi' vectors must have the same length.")
      }
      if (Gi[1] != 0) stop("First reported 'Gi' value must be zero.")
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
   prod <- data.frame(Np = Np, Rp = Rp, Wp = Wp)
   inj <- data.frame(Gi = Gi, Wi = Wi)
   if (unknown_param == "N") {
      if (m == 0) {
         if (is.null(aquifer_model)) {
            if (all(We == 0)) res_type <- "volumetric"
            if (any(We != 0)) res_type <- "water_drive"
         } else {
            res_type <- "water_drive"
         }
      } else {
         if (is.null(aquifer_model)) {
            if (all(We == 0)) res_type <- "gas_cap"
            if (any(We != 0)) res_type <- "combination"
         } else {
            res_type <- "combination"
         }
      }
   }
   if (unknown_param == "m") {
      if (is.null(aquifer_model)) {
         if (all(We == 0)) res_type <- "gas_cap"
         if (any(We != 0)) res_type <- "combination"
      } else {
         res_type <- "combination"
      }
   }
   if (unknown_param == "We") {
      if (m == 0) {
         res_type <- "water_drive"
      } else {
         res_type <- "combination"
      }
   }
   if (unknown_param == "N_m") {
      if (is.null(aquifer_model)) {
         if (all(We == 0)) res_type <- "gas_cap"
         if (any(We != 0)) res_type <- "combination"
      } else {
         res_type <- "combination"
      }
   }
   if (res_type %in% c("gas_cap", "combination")) {
      if (max(p) != pb) {
         stop("Initial reservoir pressure must be equal to 'pb' in reservoirs with an associated gas cap.")
      }
      if (wf[1] == 0) {
         stop("Information at 'pb' cannot be removed from the calculations. Change the corresponding 'wf' value to one.")
      }
   }
   if (res_type == "volumetric") {
      optim_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sorg = sorg, sorw = sorw)
      class(optim_lst) <- c("volumetric_optim_oil", "optimization_oil")
   }
   if (res_type == "gas_cap") {
      if (is.null(m) & is.null(N)) {
         optim_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sorg = sorg, sorw = sorw)
         class(optim_lst) <- c("gas_cap_optim_oil", "optimization_oil")
      }
      if (!(is.null(m)) & is.null(N)) {
         optim_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sorg = sorg, sorw = sorw)
         class(optim_lst) <- c("gas_cap_m_optim_oil", "gas_cap_optim_oil", "optimization_oil")
      }
      if (is.null(m) & !(is.null(N))) {
         optim_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sorg = sorg, sorw = sorw)
         class(optim_lst) <- c("gas_cap_N_optim_oil", "gas_cap_optim_oil", "optimization_oil")
      }
   }
   if (res_type == "water_drive") {
      if (is.null(N)) {
         optim_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sorg = sorg, sorw = sorw)
         class(optim_lst) <- c("water_drive_We_optim_oil", "water_drive_optim_oil", "optimization_oil")
      }
      if (!is.null(N)) {
         optim_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sorg = sorg, sorw = sorw)
         class(optim_lst) <- c("water_drive_N_optim_oil", "water_drive_optim_oil", "optimization_oil")
      }
   }
   if (res_type == "combination") {
      if (!is.null(We)) {
         if (is.null(N) & is.null(m)) {
            optim_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sorg = sorg, sorw = sorw)
            class(optim_lst) <- c("combination_We_optim_oil", "combination_optim_oil", "optimization_oil")
         } else {
            if (is.null(m)) {
               optim_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sorg = sorg, sorw = sorw)
               class(optim_lst) <- c("combination_We_N_optim_oil", "combination_optim_oil", "optimization_oil")
            }
            if (is.null(N)) {
               optim_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sorg = sorg, sorw = sorw)
               class(optim_lst) <- c("combination_We_m_optim_oil", "combination_optim_oil", "optimization_oil")
            }
         }
      }
      if (!is.null(aquifer_model)) {
         if (is.null(N) & is.null(m)) {
            optim_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sorg = sorg, sorw = sorw)
            class(optim_lst) <- c("combination_We_optim_oil", "combination_optim_oil", "optimization_oil")
         } else if (is.null(m)) {
            optim_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sorg = sorg, sorw = sorw)
            class(optim_lst) <- c("combination_We_N_optim_oil", "combination_optim_oil", "optimization_oil")
         } else if (is.null(N)) {
            optim_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sorg = sorg, sorw = sorw)
            class(optim_lst) <- c("combination_We_m_optim_oil", "combination_optim_oil", "optimization_oil")
         } else {
            optim_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer_lst, wf = wf, sorg = sorg, sorw = sorw)
            class(optim_lst) <- c("combination_N_m_optim_oil", "combination_optim_oil", "optimization_oil")
         }
      }
   }
   return(optim_lst)
}



# ******************************************************************************

#' Generic function for predicting unknown parameters of a material balance model
#'
#' Generate a list of class 'mbal_oil' with estimates for the unknown parameters of the material balance model according to the class of 'optim_lst' and 'time_lst' objects
#'
#' @param optim_lst a list object of class 'optimization_oil'
#' @param time_lst a list object of class 'time/date'
#'
#' @return a list of class 'mbal_oil' with estimates for the unknown parameters of the material balance model according to the class of 'optim_lst' and 'time_lst' objects
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
#' p <- c(3330, 3150, 3000, 2850, 2700, 2550, 2400)
#'
#' We <- rep(0, length.out = length(p))
#'
#' Np <- c(0, 3.295, 5.903, 8.852, 11.503, 14.513, 17.730) * 1e6
#'
#' Rp <- c(0, 1050, 1060, 1160, 1235, 1265, 1300)
#'
#' Wp <- rep(0, length.out = length(p))
#'
#' Wi <- rep(0, length.out = length(p))
#'
#' Gi <- rep(0, length.out = length(p))
#'
#' wf <- c(1, 1, 1, 0, 1, 0, 1)
#'
#' mbal_optim_oil_lst <- mbal_optim_param_oil(input_unit = "Field",
#' output_unit = "Field", unknown_param = "N_m", aquifer_model = NULL,
#' phi = 0.2, swi = 0.2, Np = Np,Rp = Rp, Wp = Wp, Gi = Gi, Wi = Wi,
#' We = We, pb = 3330, p = p, pvt = pvt_table, cf = 0, wf = wf,
#' sorg = 0.2, sorw = 0)
#'
#' time_lst <- mbal_time(c(0, 365, 730, 1095, 1460, 1825, 2190), "day")
#'
#' optim_results <- mbal_optim_oil(mbal_optim_oil_lst, time_lst)
#'
#' dplyr::glimpse(optim_results)

mbal_optim_oil <- function(optim_lst, time_lst) {

   if (inherits(optim_lst, "optimization_oil") == TRUE & inherits(time_lst, "time")) {
      UseMethod("mbal_optim_oil")
   } else {
      if (!inherits(optim_lst, "optimization_oil")) {
         stop("A class of 'optimization_oil' must be assigned to the 'optim_lst' parameter of the mbal_optim_oil() function.")
      }
      if (!inherits(time_lst, "time")) {
         stop("A class of 'time' must be assigned to the 'time_lst' parameter of the mbal_optim_oil() function.")
      }
   }
}



# ******************************************************************************

#' S3 method for class 'mbal_optim_oil'
#'
#' Generate a list of class 'mbal_oil' with estimates for the unknown parameters of a volumetric oil reservoir
#'
#' @param optim_lst a list object of class 'optimization_oil'
#' @param time_lst a list object of class 'time'
#'
#' @return a list of class 'mbal_oil' with estimates for the unknown parameters of a volumetric oil reservoir
#' @export
mbal_optim_oil.volumetric_optim_oil <- function(optim_lst, time_lst) {
   volumetric_optim_oil(optim_lst, time_lst)
}


# ******************************************************************************

#' S3 method for class 'mbal_optim_oil'
#'
#' Generate a list of class 'mbal_oil' with estimates for the unknown parameters of a gas_cap_drive oil reservoir
#'
#' @param optim_lst a list object of class 'optimization_oil'
#' @param time_lst a list object of class 'time'
#'
#' @return a list of class 'mbal_oil' with estimates for the unknown parameters of a gas_cap_drive oil reservoir
#' @export
mbal_optim_oil.gas_cap_optim_oil <- function(optim_lst, time_lst) {
   gas_cap_optim_oil(optim_lst, time_lst)
}


# ******************************************************************************

#' S3 method for class 'mbal_optim_oil'
#'
#' Generate a list of class 'mbal_oil' with estimates for the unknown parameters of a water_drive oil reservoir
#'
#' @param optim_lst a list object of class 'optimization_oil'
#' @param time_lst a list object of class 'time'
#'
#' @return a list of class 'mbal_oil' with estimates for the unknown parameters of a water_drive oil reservoir
#' @export
mbal_optim_oil.water_drive_optim_oil <- function(optim_lst, time_lst) {
   water_drive_optim_oil(optim_lst, time_lst)
}

# ******************************************************************************

#' S3 method for class 'mbal_optim_oil'
#'
#' Generate a list of class 'mbal_oil' with estimates for the unknown parameters of a combination_drive oil reservoir
#'
#' @param optim_lst a list object of class 'optimization_oil'
#' @param time_lst a list object of class 'time'
#'
#' @return a list of class 'mbal_oil' with estimates for the unknown parameters of a combination_drive oil reservoir
#' @export
mbal_optim_oil.combination_optim_oil <- function(optim_lst, time_lst) {
   combination_optim_oil(optim_lst, time_lst)
}





# *************** Volumetric **************

volumetric_optim_oil <- function(optim_lst, time_lst) {

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
   N <- optim_lst$N
   m <- optim_lst$m
   phi <- optim_lst$phi
   swi <- optim_lst$swi
   pb <- optim_lst$pb
   p <- optim_lst$p
   l <- length(p)
   if (length(time_lst$t) != length(p)) stop("Lenght of 'time_lst$t' vector must be equal to the length of 'p' vector.")
   cf <- optim_lst$cf
   pvt <- optim_lst$pvt
   prod <- optim_lst$prod
   inj <- optim_lst$inj
   We <- optim_lst$We
   aquifer <- optim_lst$aquifer
   cls <- class(optim_lst)[1]
   aqu_cls <- class(aquifer)[1]
   if (aqu_cls == "NoA") We <- aquifer$We
   wf <- optim_lst$wf
   sorw <- optim_lst$sorw
   sorg <- optim_lst$sorg
   if (is.null(sorw)) stop("'sorw' must be a numeric value.")
   if (is.null(sorg)) stop("'sorg' must be a numeric value.")
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
   Rp <- prod_$Rp
   Wp <- prod_$Wp
   Gi <- inj_$Gi
   Wi <- inj_$Wi
   Gp <- Np * Rp
   dp <- vector(length = l)
   Bto <- vector(length = l)
   Eo <- vector(length = l)
   Eowf <- vector(length = l)
   Eg <- vector(length = l)
   Egwf <- vector(length = l)
   Ew <- vector(length = l)
   Ef <- vector(length = l)
   F_ <- vector(length = l)
   dp <- p_[1] - p_
   Bto <- (Bo * (1 - Rs[1] * Rv) + Bg * (Rs[1] - Rs)) / (1 - Rs * Rv)
   Btg <- (Bg * (1 - Rs * Rv[1]) + Bo * (Rv[1] - Rv)) / (1 - Rs * Rv)
   Eo <- Bto - Bto[1]
   Eg <- Btg - Btg[1]
   Ew <- Bw - Bw[1]
   for (i in 1:l) {
      Ef[i] <- -1 * trapz(p_[1:i], cf_[1:i])
   }
   Eowf <- Eo + Bo[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
   Egwf <- Eg + Bg[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
   F_ <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_
   data <- data.frame(Eowf = Eowf, `F_` = F_)
   mbal_lm <- lm(`F_` ~ 0 + Eowf, data = data)
   Nfoi <- mbal_lm$coefficients
   names(Nfoi) <- NULL
   Gfgi <- 0
   N <- Nfoi + Rv[1] * Gfgi
   G <- Gfgi + Rs[1] * Nfoi
   m <- Gfgi * Bg[1] / (Nfoi * Bo[1])
   PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
   BV <- PV / phi
   W <- PV * swi / Bw[1]
   Et <- Gfgi * Eg + Nfoi * Eo + W * Ew + PV * Ef + We_ + (Wi - Wp) * Bw
   STDERR <- sqrt(sum((F_ - Nfoi * Eowf) * ((F_ - Nfoi * Eowf))) / (l - 2))
   plot(Eowf, F_, pch = 21, bg = "blue" , xlab = "Eowf", ylab = "F")
   abline(a = 0, b = Nfoi)
   mtext(paste("slope = Nfoi = ", round(Nfoi,0), " , STDERR = ", round(STDERR,2)), side = 3)
   xmin <- min(Eowf, na.rm = TRUE) + 0.25 * diff(range(Eowf, na.rm = TRUE))
   ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
   text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F - N[foi] * E[owf])^2), n - 2))))
   p_calc <- vector(length = l)
   for (i in 1:l) {
      if (i == 1) {
         p_calc[i] <- p_[i]
      } else {
         fun <- function(x, pb, p_ini, Bo_ini, Bto_ini, Rs_ini, Bw_ini, Nfoi, swi, p, Ef, Np, Rp, Gp, Wp, Wi, Gi, We, pvt) {
            Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear")$y
            Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear")$y
            Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear")$y
            Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear")$y
            Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear")$y
            Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear")$y
            dp <- p_ini - x
            Bto <- (Bo * (1 - Rs_ini * Rv) + Bg * (Rs_ini - Rs)) / (1 - Rs * Rv)
            Eo <- Bto - Bto_ini
            Ew <- Bw - Bw_ini
            Eowf <- Eo + Bo_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
            F_lhs_calc <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
            F_rhs_calc <- Nfoi * Eowf
            diff <- F_lhs_calc - F_rhs_calc
            return(diff)
         }
         p_calc[i] <- fzero(fun = fun, x = p_[i], pb = pb, p_ini = p_[1], Bo_ini = Bo[1], Bto_ini = Bto[1], Rs_ini = Rs[1], Bw_ini = Bw[1], Nfoi = Nfoi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Rp = Rp[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], Gi = Gi[i], We = We_[i], pvt = pvt, maxiter = 1000)$x
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

   ymin <- max(0, min((F_[2:l] + We_[2:l]) / Eowf[2:l]) / 1.2)
   ymax <- 1.2 * max((F_[2:l] + We_[2:l]) / Eowf[2:l])
   plot(Np, (F_ + We_) / Eowf, pch = 21, bg = "blue" , xlab = "Np", ylab = "F/Eowf", ylim = c(ymin, ymax))

   mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sorw = sorw, sorg = sorg, p_est = p_calc)
   class(mbal_lst) <- c("volumetric_oil", "mbal_oil")
   return(mbal_lst)
}




# *************** Gas Cap **************

gas_cap_optim_oil <- function(optim_lst, time_lst) {

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
   N <- optim_lst$N
   m <- optim_lst$m
   phi <- optim_lst$phi
   swi <- optim_lst$swi
   pb <- optim_lst$pb
   p <- optim_lst$p
   l <- length(p)
   if (length(time_lst$t) != length(p)) stop("Lenght of 'time_lst$t' vector must be equal to the length of 'p' vector.")
   cf <- optim_lst$cf
   pvt <- optim_lst$pvt
   prod <- optim_lst$prod
   inj <- optim_lst$inj
   We <- optim_lst$We
   aquifer <- optim_lst$aquifer
   cls <- class(optim_lst)[1]
   aqu_cls <- class(aquifer)[1]
   if (aqu_cls == "NoA") We <- aquifer$We
   wf <- optim_lst$wf
   sorw <- optim_lst$sorw
   sorg <- optim_lst$sorg
   if (is.null(sorw)) stop("'sorw' must be a numeric value.")
   if (is.null(sorg)) stop("'sorg' must be a numeric value.")
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
   Rp <- prod_$Rp
   Wp <- prod_$Wp
   Gi <- inj_$Gi
   Wi <- inj_$Wi
   Gp <- Np * Rp
   dp <- vector(length = l)
   Bto <- vector(length = l)
   Eo <- vector(length = l)
   Eowf <- vector(length = l)
   Eg <- vector(length = l)
   Egwf <- vector(length = l)
   Ew <- vector(length = l)
   Ef <- vector(length = l)
   Et <- vector(length = l)
   F_ <- vector(length = l)

   if (cls == "gas_cap_optim_oil") {
      dp <- p_[1] - p_
      Bto <- (Bo * (1 - Rs[1] * Rv) + Bg * (Rs[1] - Rs)) / (1 - Rs * Rv)
      Btg <- (Bg * (1 - Rs * Rv[1]) + Bo * (Rv[1] - Rv)) / (1 - Rs * Rv)
      Eo <- Bto - Bto[1]
      Eg <- Btg - Btg[1]
      Ew <- Bw - Bw[1]
      for (i in 1:l) {
         Ef[i] <- -1 * trapz(p_[1:i], cf_[1:i])
      }
      Eowf <- Eo + Bo[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      Egwf <- Eg + Bg[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      F_ <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_
      data <- data.frame(Eowf = Eowf, Egwf = Egwf, `F_` = F_)
      mbal_lm <- lm(`F_` ~ 0 + Eowf + Egwf, data = data)
      Nfoi <- mbal_lm$coefficients[1]
      names(Nfoi) <- NULL
      Gfgi <- mbal_lm$coefficients[2]
      names(Gfgi) <- NULL
      N <- Nfoi + Rv[1] * Gfgi
      G <- Gfgi + Rs[1] * Nfoi
      m <- Gfgi * Bg[1] / (Nfoi * Bo[1])
      PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
      BV <- PV / phi
      W <- PV * swi / Bw[1]
      Et <- Eowf + m * Egwf * (Bo[1] / Bg[1])
      STDERR <- sqrt(sum((F_ - Nfoi * Eowf - Gfgi * Egwf) * ((F_ - Nfoi * Eowf - Gfgi * Egwf))) / (l - 2))
      plot(Et, F_, pch = 21, bg = "blue" , xlab = "Et", ylab = "F")
      abline(a = 0, b = Nfoi)
      mtext(paste("slope = Nfoi = ", round(Nfoi,0), ", m = ", round(m,3), " , STDERR = ", round(STDERR,2)), side = 3)
      xmin <- min(Et, na.rm = TRUE) + 0.25 * diff(range(Et, na.rm = TRUE))
      ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
      text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F -  N[foi] * E[owf] - G[fgi] * E[gwf])^2), n - 2))))
      p_calc <- vector(length = l)
      for (i in 1:l) {
         if (i == 1) {
            p_calc[i] <- p_[i]
         } else {
            fun <- function(x, pb, p_ini, Bo_ini, Bto_ini, Rs_ini, Bg_ini, Btg_ini, Rv_ini, Bw_ini, Nfoi, Gfgi, swi, p, Ef, Np, Rp, Gp, Wp, Wi, Gi, We, pvt) {
               Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear")$y
               Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear")$y
               Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear")$y
               Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear")$y
               Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear")$y
               Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear")$y
               dp <- p_ini - x
               Bto <- (Bo * (1 - Rs_ini * Rv) + Bg * (Rs_ini - Rs)) / (1 - Rs * Rv)
               Btg <- (Bg * (1 - Rs * Rv_ini) + Bo * (Rv_ini - Rv)) / (1 - Rs * Rv)
               Eo <- Bto - Bto_ini
               Eg <- Btg - Btg_ini
               Ew <- Bw - Bw_ini
               Eowf <- Eo + Bo_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
               Egwf <- Eg + Bg_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
               F_lhs_calc <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
               F_rhs_calc <- Nfoi * Eowf + Gfgi * Egwf
               diff <- F_lhs_calc - F_rhs_calc
               return(diff)
            }
            p_calc[i] <- fzero(fun = fun, x = p_[i], pb = pb, p_ini = p_[1], Bo_ini = Bo[1], Bto_ini = Bto[1], Rs_ini = Rs[1], Bg_ini = Bg[1], Btg_ini = Btg[1], Rv_ini = Rv[1], Bw_ini = Bw[1], Nfoi = Nfoi, Gfgi = Gfgi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Rp = Rp[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], Gi = Gi[i], We = We_[i], pvt = pvt, maxiter = 1000)$x
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

      ymin <- max(0, min((F_[2:l] + We_[2:l]) / Eowf[2:l]) / 1.2)
      ymax <- 1.2 * max((F_[2:l] + We_[2:l]) / Eowf[2:l])
      plot(Np, (F_ + We_) / Eowf, pch = 21, bg = "blue" , xlab = "Np", ylab = "F/Eowf", ylim = c(ymin, ymax))

      mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sorw = sorw, sorg = sorg, p_est = p_calc)
      class(mbal_lst) <- c("gas_cap_oil", "mbal_oil")
      return(mbal_lst)
   }


   if (cls == "gas_cap_m_optim_oil") {
      dp <- p_[1] - p_
      Bto <- (Bo * (1 - Rs[1] * Rv) + Bg * (Rs[1] - Rs)) / (1 - Rs * Rv)
      Btg <- (Bg * (1 - Rs * Rv[1]) + Bo * (Rv[1] - Rv)) / (1 - Rs * Rv)
      Eo <- Bto - Bto[1]
      Eg <- Btg - Btg[1]
      Ew <- Bw - Bw[1]
      for (i in 1:l) {
         Ef[i] <- -1 * trapz(p_[1:i], cf_[1:i])
      }
      Eowf <- Eo + Bo[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      Egwf <- Eg + Bg[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      Et <- Eowf + m * Egwf * (Bo[1] / Bg[1])
      F_ <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_
      data <- data.frame(Et = Et, `F_` = F_)
      mbal_lm <- lm(`F_` ~ 0 + Et, data = data)
      Nfoi <- mbal_lm$coefficients[1]
      names(Nfoi) <- NULL
      Gfgi <- m * Nfoi * Bo[1] / Bg[1]
      N <- Nfoi + Rv[1] * Gfgi
      G <- Gfgi + Rs[1] * Nfoi
      PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
      BV <- PV / phi
      W <- PV * swi / Bw[1]
      STDERR <- sqrt(sum((F_ - Nfoi * Eowf - Gfgi * Egwf) * ((F_ - Nfoi * Eowf - Gfgi * Egwf))) / (l - 2))
      plot(Et, F_, pch = 21, bg = "blue" , xlab = "Et", ylab = "F")
      abline(a = 0, b = Nfoi)
      mtext(paste("slope = Nfoi = ", round(Nfoi,0), ", m = ", round(m,3), " , STDERR = ", round(STDERR,2)), side = 3)
      xmin <- min(Et, na.rm = TRUE) + 0.25 * diff(range(Et, na.rm = TRUE))
      ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
      text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F -  N[foi] * E[owf] - G[fgi] * E[gwf])^2), n - 2))))
      p_calc <- vector(length = l)
      for (i in 1:l) {
         if (i == 1) {
            p_calc[i] <- p_[i]
         } else {
            fun <- function(x, pb, p_ini, Bo_ini, Bto_ini, Rs_ini, Bg_ini, Btg_ini, Rv_ini, Bw_ini, Nfoi, Gfgi, swi, p, Ef, Np, Rp, Gp, Wp, Wi, Gi, We, pvt) {
               Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear")$y
               Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear")$y
               Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear")$y
               Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear")$y
               Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear")$y
               Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear")$y
               dp <- p_ini - x
               Bto <- (Bo * (1 - Rs_ini * Rv) + Bg * (Rs_ini - Rs)) / (1 - Rs * Rv)
               Btg <- (Bg * (1 - Rs * Rv_ini) + Bo * (Rv_ini - Rv)) / (1 - Rs * Rv)
               Eo <- Bto - Bto_ini
               Eg <- Btg - Btg_ini
               Ew <- Bw - Bw_ini
               Eowf <- Eo + Bo_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
               Egwf <- Eg + Bg_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
               F_lhs_calc <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
               F_rhs_calc <- Nfoi * Eowf + Gfgi * Egwf
               diff <- F_lhs_calc - F_rhs_calc
               return(diff)
            }
            p_calc[i] <- fzero(fun = fun, x = p_[i], pb = pb, p_ini = p_[1], Bo_ini = Bo[1], Bto_ini = Bto[1], Rs_ini = Rs[1], Bg_ini = Bg[1], Btg_ini = Btg[1], Rv_ini = Rv[1], Bw_ini = Bw[1], Nfoi = Nfoi, Gfgi = Gfgi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Rp = Rp[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], Gi = Gi[i], We = We_[i], pvt = pvt, maxiter = 1000)$x
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

      ymin <- max(0, min((F_[2:l] + We_[2:l]) / Eowf[2:l]) / 1.2)
      ymax <- 1.2 * max((F_[2:l] + We_[2:l]) / Eowf[2:l])
      plot(Np, (F_ + We_) / Eowf, pch = 21, bg = "blue" , xlab = "Np", ylab = "F/Eowf", ylim = c(ymin, ymax))

      mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sorw = sorw, sorg = sorg, p_est = p_calc)
      class(mbal_lst) <- c("gas_cap_oil", "mbal_oil")
      return(mbal_lst)
   }

   if (cls == "gas_cap_N_optim_oil") {
      dp <- p_[1] - p_
      Bto <- (Bo * (1 - Rs[1] * Rv) + Bg * (Rs[1] - Rs)) / (1 - Rs * Rv)
      Btg <- (Bg * (1 - Rs * Rv[1]) + Bo * (Rv[1] - Rv)) / (1 - Rs * Rv)
      Eo <- Bto - Bto[1]
      Eg <- Btg - Btg[1]
      Ew <- Bw - Bw[1]
      for (i in 1:l) {
         Ef[i] <- -1 * trapz(p_[1:i], cf_[1:i])
      }
      for (i in 1:l) {
         Ef[i] <- -1 * trapz(p_[1:i], cf_[1:i])
      }
      Eowf <- Eo + Bo[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      Egwf <- Eg + Bg[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      F_ <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_
      data <- data.frame(En = (Egwf - Rv[1] * Eowf) / Eowf, `F_Eowf_N` = F_ / Eowf - N)
      mbal_lm <- lm(`F_Eowf_N` ~ 0 + En, data = data)
      Gfgi <- mbal_lm$coefficients[1]
      names(Gfgi) <- NULL
      Nfoi <- (N - Gfgi * Rv[1])
      m <- Gfgi * Bg[1] / (Nfoi * Bo[1])
      G <- Gfgi + Rs[1] * Nfoi
      PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
      BV <- PV / phi
      W <- PV * swi / Bw[1]
      Et <- Eowf + m * Egwf * (Bo[1] / Bg[1])
      STDERR <- sqrt(sum((F_ - Nfoi * Eowf - Gfgi * Egwf) * ((F_ - Nfoi * Eowf - Gfgi * Egwf))) / (l - 2))
      plot(Et, F_, pch = 21, bg = "blue" , xlab = "Et", ylab = "F")
      abline(a = 0, b = Nfoi)
      mtext(paste("slope = Nfoi = ", round(Nfoi,0), ", m = ", round(m,3), " , STDERR = ", round(STDERR,2)), side = 3)
      xmin <- min(Et, na.rm = TRUE) + 0.25 * diff(range(Et, na.rm = TRUE))
      ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
      text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F -  N[foi] * E[owf] - G[fgi] * E[gwf])^2), n - 2))))
      p_calc <- vector(length = l)
      for (i in 1:l) {
         if (i == 1) {
            p_calc[i] <- p_[i]
         } else {
            fun <- function(x, pb, p_ini, Bo_ini, Bto_ini, Rs_ini, Bg_ini, Btg_ini, Rv_ini, Bw_ini, Nfoi, Gfgi, swi, p, Ef, Np, Rp, Gp, Wp, Wi, Gi, We, pvt) {
               Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear")$y
               Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear")$y
               Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear")$y
               Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear")$y
               Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear")$y
               Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear")$y
               dp <- p_ini - x
               Bto <- (Bo * (1 - Rs_ini * Rv) + Bg * (Rs_ini - Rs)) / (1 - Rs * Rv)
               Btg <- (Bg * (1 - Rs * Rv_ini) + Bo * (Rv_ini - Rv)) / (1 - Rs * Rv)
               Eo <- Bto - Bto_ini
               Eg <- Btg - Btg_ini
               Ew <- Bw - Bw_ini
               Eowf <- Eo + Bo_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
               Egwf <- Eg + Bg_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
               F_lhs_calc <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
               F_rhs_calc <- Nfoi * Eowf + Gfgi * Egwf
               diff <- F_lhs_calc - F_rhs_calc
               return(diff)
            }
            p_calc[i] <- fzero(fun = fun, x = p_[i], pb = pb, p_ini = p_[1], Bo_ini = Bo[1], Bto_ini = Bto[1], Rs_ini = Rs[1], Bg_ini = Bg[1], Btg_ini = Btg[1], Rv_ini = Rv[1], Bw_ini = Bw[1], Nfoi = Nfoi, Gfgi = Gfgi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Rp = Rp[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], Gi = Gi[i], We = We_[i], pvt = pvt, maxiter = 1000)$x
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

      ymin <- max(0, min((F_[2:l] + We_[2:l]) / Eowf[2:l]) / 1.2)
      ymax <- 1.2 * max((F_[2:l] + We_[2:l]) / Eowf[2:l])
      plot(Np, (F_ + We_) / Eowf, pch = 21, bg = "blue" , xlab = "Np", ylab = "F/Eowf", ylim = c(ymin, ymax))

      mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sorw = sorw, sorg = sorg, p_est = p_calc)
      class(mbal_lst) <- c("gas_cap_oil", "mbal_oil")
      return(mbal_lst)
   }
}





# *************** Water Drive **************

water_drive_optim_oil <- function(optim_lst, time_lst) {

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
   N <- optim_lst$N
   m <- optim_lst$m
   phi <- optim_lst$phi
   swi <- optim_lst$swi
   pb <- optim_lst$pb
   p <- optim_lst$p
   l <- length(p)
   if (length(time_lst$t) != length(p)) stop("Lenght of 'time_lst$t' vector must be equal to the length of 'p' vector.")
   cf <- optim_lst$cf
   pvt <- optim_lst$pvt
   prod <- optim_lst$prod
   inj <- optim_lst$inj
   We <- optim_lst$We
   aquifer <- optim_lst$aquifer
   cls <- class(optim_lst)[1]
   aqu_cls <- class(aquifer)[1]
   if (aqu_cls == "We") We <- aquifer$We
   wf <- optim_lst$wf
   sorw <- optim_lst$sorw
   sorg <- optim_lst$sorg
   if (is.null(sorw)) stop("'sorw' must be a numeric value.")
   if (is.null(sorg)) stop("'sorg' must be a numeric value.")
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
   Rp <- prod_$Rp
   Wp <- prod_$Wp
   Gi <- inj_$Gi
   Wi <- inj_$Wi
   Gp <- Np * Rp
   dp <- vector(length = l)
   Bto <- vector(length = l)
   Eo <- vector(length = l)
   Eowf <- vector(length = l)
   Eg <- vector(length = l)
   Egwf <- vector(length = l)
   Ew <- vector(length = l)
   Ef <- vector(length = l)
   Et <- vector(length = l)
   F_ <- vector(length = l)


   if (cls == "water_drive_We_optim_oil") {
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
      Eowf <- Eo + Bo[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      Egwf <- Eg + Bg[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      F_ <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_
      data <- data.frame(Eowf = Eowf, `F_` = F_)
      mbal_lm <- lm(`F_` ~ 0 + Eowf, data = data)
      Nfoi <- mbal_lm$coefficients
      names(Nfoi) <- NULL
      Gfgi <- 0
      N <- Nfoi + Rv[1] * Gfgi
      G <- Gfgi + Rs[1] * Nfoi
      m <- Gfgi * Bg[1] / (Nfoi * Bo[1])
      PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
      BV <- PV / phi
      W <- PV * swi / Bw[1]
      Et <- Gfgi * Eg + Nfoi * Eo + W * Ew + PV * Ef + We_ + (Wi - Wp) * Bw
      STDERR <- sqrt(sum((F_ - Nfoi * Eowf) * ((F_ - Nfoi * Eowf))) / (l - 2))
      plot(Eowf, F_, pch = 21, bg = "blue" , xlab = "Eowf", ylab = "F")
      abline(a = 0, b = Nfoi)
      mtext(paste("slope = Nfoi = ", round(Nfoi,0), " , STDERR = ", round(STDERR,2)), side = 3)
      xmin <- min(Eowf, na.rm = TRUE) + 0.25 * diff(range(Eowf, na.rm = TRUE))
      ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
      text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F - W[e] - N[foi] * E[owf])^2), n - 2))))
      p_calc <- vector(length = l)
      for (i in 1:l) {
         if (i == 1) {
            p_calc[i] <- p_[i]
         } else {
            fun <- function(x, pb, p_ini, Bo_ini, Bto_ini, Rs_ini, Bw_ini, Nfoi, swi, p, Ef, Np, Rp, Gp, Wp, Wi, Gi, We, pvt) {
               Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear")$y
               Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear")$y
               Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear")$y
               Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear")$y
               Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear")$y
               Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear")$y
               dp <- p_ini - x
               Bto <- (Bo * (1 - Rs_ini * Rv) + Bg * (Rs_ini - Rs)) / (1 - Rs * Rv)
               Eo <- Bto - Bto_ini
               Ew <- Bw - Bw_ini
               Eowf <- Eo + Bo_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
               F_lhs_calc <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
               F_rhs_calc <- Nfoi * Eowf
               diff <- F_lhs_calc - F_rhs_calc
               return(diff)
            }
            p_calc[i] <- fzero(fun = fun, x = p_[i], pb = pb, p_ini = p_[1], Bo_ini = Bo[1], Bto_ini = Bto[1], Rs_ini = Rs[1], Bw_ini = Bw[1], Nfoi = Nfoi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Rp = Rp[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], Gi = Gi[i], We = We_[i], pvt = pvt, maxiter = 1000)$x
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

      ymin <- max(0, min((F_[2:l] + We_[2:l]) / Eowf[2:l]) / 1.2)
      ymax <- 1.2 * max((F_[2:l] + We_[2:l]) / Eowf[2:l])
      plot(Np, (F_ + We_) / Eowf, pch = 21, bg = "blue" , xlab = "Np", ylab = "F/Eowf", ylim = c(ymin, ymax))

      mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sorw = sorw, sorg = sorg, p_est = p_calc)
      class(mbal_lst) <- c("water_drive_oil", "mbal_oil")
      return(mbal_lst)
   }


   if (cls == "water_drive_N_optim_oil") {
      dp <- p_[1] - p_
      Bto <- (Bo * (1 - Rs[1] * Rv) + Bg * (Rs[1] - Rs)) / (1 - Rs * Rv)
      Btg <- (Bg * (1 - Rs * Rv[1]) + Bo * (Rv[1] - Rv)) / (1 - Rs * Rv)
      Eo <- Bto - Bto[1]
      Eg <- Btg - Btg[1]
      Ew <- Bw - Bw[1]
      for (i in 1:l) {
         Ef[i] <- -1 * trapz(p_[1:i], cf_[1:i])
      }
      Eowf <- Eo + Bo[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      Egwf <- Eg + Bg[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      Gfgi <- 0
      Nfoi <- N - Rv[1] * Gfgi
      F_ <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw
      We_ <- F_ - Nfoi * Eowf - Gfgi * Egwf
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
         error_fun_ussre_oil <- function(par, We, phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- veh_uss_rad_edge(phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_ussre_oil, We = We_, phi = phi, perm_h = perm_h, h_a = h_a, r_a = r_a, r_R = r_R,
                           tetha = tetha, mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_,
                           pressure = pressure_, lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "uss_rad_bottom") {
         error_fun_ussrb_oil <- function(par, We, phi, perm_h, perm_v, h_a, r_a, r_R, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- yk_uss_rad_bottom(phi, perm_h, perm_v, h_a, r_a, r_R, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_ussrb_oil, We = We_, phi = phi, perm_h = perm_h, perm_v = perm_v, h_a = h_a, r_a = r_a,
                           r_R = r_R, mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_, pressure = pressure_,
                           lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "pss_rad_edge") {
         error_fun_pssre_oil <- function(par, We, phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- fetkovich_pss_rad_edge(phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_pssre_oil, We = We_, phi = phi, perm_h = perm_h, h_a = h_a, r_a = r_a,
                           r_R = r_R, tetha = tetha, mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_,
                           pressure = pressure_, lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "uss_lin_edge") {
         error_fun_ussle_oil <- function(par, We, phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- nb_uss_lin_edge(phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_ussle_oil, We = We_, phi = phi, perm_h = perm_h, h_a = h_a, w_a = w_a, l_a = l_a,
                           mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_, pressure = pressure_,
                           lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "uss_lin_bottom") {
         error_fun_usslb_oil <- function(par, We, phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- nb_uss_lin_bottom(phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_usslb_oil, We = We_, phi = phi, perm_v = perm_v, h_a = h_a, w_a = w_a, l_a = l_a,
                           mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_, pressure = pressure_,
                           lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "pss_lin_edge") {
         error_fun_pssle_oil <- function(par, We, phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- fetk_pss_lin_edge(phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_pssle_oil, We = We_, phi = phi, perm_h = perm_h, h_a = h_a, w_a = w_a, l_a = l_a,
                           mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_, pressure = pressure_,
                           lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "pss_lin_bottom") {
         error_fun_psslb_oil <- function(par, We, phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- fetk_pss_lin_bottom(phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_psslb_oil, We = We_, phi = phi, perm_v = perm_v, h_a = h_a, w_a = w_a, l_a = l_a,
                           mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_, pressure = pressure_,
                           lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "pot") {
         error_fun_pot_oil <- function(par, We, phi, h_a, r_a, r_R, tetha, c_water, c_rock, time, pressure) {
            We_calc <- pot(phi, h_a, r_a, r_R, tetha, c_water, c_rock, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_pot_oil, We = We_, phi = phi, h_a = h_a, r_a = r_a, r_R = r_R,
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
      F_ <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_calc
      Gfgi <- 0
      N <- Nfoi + Rv[1] * Gfgi
      G <- Gfgi + Rs[1] * Nfoi
      m <- Gfgi * Bg[1] / (Nfoi * Bo[1])
      PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
      BV <- PV / phi
      W <- PV * swi / Bw[1]
      Et <- Gfgi * Eg + Nfoi * Eo + W * Ew + PV * Ef + We_calc + (Wi - Wp) * Bw
      STDERR <- sqrt(sum((F_ - Nfoi * Eowf) * ((F_ - Nfoi * Eowf))) / (l - 2))
      plot(Eowf, F_, pch = 21, bg = "blue" , xlab = "Eowf", ylab = "F")
      abline(a = 0, b = Nfoi)
      mtext(paste("slope = Nfoi = ", round(Nfoi,0), " , STDERR = ", round(STDERR,2)), side = 3)
      xmin <- min(Eowf, na.rm = TRUE) + 0.25 * diff(range(Eowf, na.rm = TRUE))
      ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
      text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F - W[e] - N[foi] * E[owf])^2), n - 2))))
      p_calc <- vector(length = l)
      for (i in 1:l) {
         if (i == 1) {
            p_calc[i] <- p_[i]
         } else {
            fun <- function(x, pb, p_ini, Bo_ini, Bto_ini, Rs_ini, Bw_ini, Nfoi, swi, p, Ef, Np, Rp, Gp, Wp, Wi, Gi, We, pvt) {
               Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear")$y
               Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear")$y
               Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear")$y
               Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear")$y
               Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear")$y
               Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear")$y
               dp <- p_ini - x
               Bto <- (Bo * (1 - Rs_ini * Rv) + Bg * (Rs_ini - Rs)) / (1 - Rs * Rv)
               Eo <- Bto - Bto_ini
               Ew <- Bw - Bw_ini
               Eowf <- Eo + Bo_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
               F_lhs_calc <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
               F_rhs_calc <- Nfoi * Eowf
               diff <- F_lhs_calc - F_rhs_calc
               return(diff)
            }
            p_calc[i] <- fzero(fun = fun, x = p_[i], pb = pb, p_ini = p_[1], Bo_ini = Bo[1], Bto_ini = Bto[1], Rs_ini = Rs[1], Bw_ini = Bw[1], Nfoi = Nfoi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Rp = Rp[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], Gi = Gi[i], We = We_calc[i], pvt = pvt, maxiter = 1000)$x
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

      ymin <- max(0, min((F_[2:l] + We_calc[2:l]) / Eowf[2:l]) / 1.2)
      ymax <- 1.2 * max((F_[2:l] + We_calc[2:l]) / Eowf[2:l])
      plot(Np, (F_ + We_calc) / Eowf, pch = 21, bg = "blue" , xlab = "Np", ylab = "F/Eowf", ylim = c(ymin, ymax))

      mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sorw = sorw, sorg = sorg, p_est = p_calc)
      class(mbal_lst) <- c("water_drive_oil", "mbal_oil")
      return(mbal_lst)
   }
}







# *************** Combination **************

combination_optim_oil <- function(optim_lst, time_lst) {

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
   N <- optim_lst$N
   m <- optim_lst$m
   phi <- optim_lst$phi
   swi <- optim_lst$swi
   pb <- optim_lst$pb
   p <- optim_lst$p
   l <- length(p)
   if (length(time_lst$t) != length(p)) stop("Lenght of 'time_lst$t' vector must be equal to the length of 'p' vector.")
   cf <- optim_lst$cf
   pvt <- optim_lst$pvt
   prod <- optim_lst$prod
   inj <- optim_lst$inj
   We <- optim_lst$We
   aquifer <- optim_lst$aquifer
   cls <- class(optim_lst)[1]
   aqu_cls <- class(aquifer)[1]
   if (aqu_cls == "We") We <- aquifer$We
   wf <- optim_lst$wf
   sorw <- optim_lst$sorw
   sorg <- optim_lst$sorg
   if (is.null(sorw)) stop("'sorw' must be a numeric value.")
   if (is.null(sorg)) stop("'sorg' must be a numeric value.")
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
   Rp <- prod_$Rp
   Wp <- prod_$Wp
   Gi <- inj_$Gi
   Wi <- inj_$Wi
   Gp <- Np * Rp
   dp <- vector(length = l)
   Bto <- vector(length = l)
   Eo <- vector(length = l)
   Eowf <- vector(length = l)
   Eg <- vector(length = l)
   Egwf <- vector(length = l)
   Ew <- vector(length = l)
   Ef <- vector(length = l)
   Et <- vector(length = l)
   F_ <- vector(length = l)

   if (cls == "combination_We_N_optim_oil") {
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
      Eowf <- Eo + Bo[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      Egwf <- Eg + Bg[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      F_ <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_
      data <- data.frame(En = (Egwf - Rv[1] * Eowf) / Eowf, `F_Eowf_N` = F_ / Eowf - N)
      mbal_lm <- lm(`F_Eowf_N` ~ 0 + En, data = data)
      Gfgi <- mbal_lm$coefficients[1]
      names(Gfgi) <- NULL
      Nfoi <- (N - Gfgi * Rv[1])
      m <- Gfgi * Bg[1] / (Nfoi * Bo[1])
      G <- Gfgi + Rs[1] * Nfoi
      PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
      BV <- PV / phi
      W <- PV * swi / Bw[1]
      Et <- Eowf + m * Egwf * (Bo[1] / Bg[1])
      STDERR <- sqrt(sum((F_ - Nfoi * Eowf - Gfgi * Egwf) * ((F_ - Nfoi * Eowf - Gfgi * Egwf))) / (l - 2))
      plot(Et, F_, pch = 21, bg = "blue" , xlab = "Et", ylab = "F")
      abline(a = 0, b = Nfoi)
      mtext(paste("slope = Nfoi = ", round(Nfoi,0), ", m = ", round(m,3), " , STDERR = ", round(STDERR,2)), side = 3)
      xmin <- min(Et, na.rm = TRUE) + 0.25 * diff(range(Et, na.rm = TRUE))
      ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
      text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F - W[e] - N[foi] * E[owf] - G[fgi] * E[gwf])^2), n - 2))))
      p_calc <- vector(length = l)
      for (i in 1:l) {
         if (i == 1) {
            p_calc[i] <- p_[i]
         } else {
            fun <- function(x, pb, p_ini, Bo_ini, Bto_ini, Rs_ini, Bg_ini, Btg_ini, Rv_ini, Bw_ini, Nfoi, Gfgi, swi, p, Ef, Np, Rp, Gp, Wp, Wi, Gi, We, pvt) {
               Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear")$y
               Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear")$y
               Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear")$y
               Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear")$y
               Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear")$y
               Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear")$y
               dp <- p_ini - x
               Bto <- (Bo * (1 - Rs_ini * Rv) + Bg * (Rs_ini - Rs)) / (1 - Rs * Rv)
               Btg <- (Bg * (1 - Rs * Rv_ini) + Bo * (Rv_ini - Rv)) / (1 - Rs * Rv)
               Eo <- Bto - Bto_ini
               Eg <- Btg - Btg_ini
               Ew <- Bw - Bw_ini
               Eowf <- Eo + Bo_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
               Egwf <- Eg + Bg_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
               F_lhs_calc <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
               F_rhs_calc <- Nfoi * Eowf + Gfgi * Egwf
               diff <- F_lhs_calc - F_rhs_calc
               return(diff)
            }
            p_calc[i] <- fzero(fun = fun, x = p_[i], pb = pb, p_ini = p_[1], Bo_ini = Bo[1], Bto_ini = Bto[1], Rs_ini = Rs[1], Bg_ini = Bg[1], Btg_ini = Btg[1], Rv_ini = Rv[1], Bw_ini = Bw[1], Nfoi = Nfoi, Gfgi = Gfgi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Rp = Rp[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], Gi = Gi[i], We = We_[i], pvt = pvt, maxiter = 1000)$x
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

      ymin <- max(0, min((F_[2:l] + We_[2:l]) / Eowf[2:l]) / 1.2)
      ymax <- 1.2 * max((F_[2:l] + We_[2:l]) / Eowf[2:l])
      plot(Np, (F_ + We_) / Eowf, pch = 21, bg = "blue" , xlab = "Np", ylab = "F/Eowf", ylim = c(ymin, ymax))

      mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sorw = sorw, sorg = sorg, p_est = p_calc)
      class(mbal_lst) <- c("combination_oil", "mbal_oil")
      return(mbal_lst)
   }
   if (cls == "combination_We_m_optim_oil") {
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
      Eowf <- Eo + Bo[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      Egwf <- Eg + Bg[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      Et <- Eowf + m * Egwf * (Bo[1] / Bg[1])
      F_ <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_
      data <- data.frame(Et = Et, `F_` = F_)
      mbal_lm <- lm(`F_` ~ 0 + Et, data = data)
      Nfoi <- mbal_lm$coefficients[1]
      names(Nfoi) <- NULL
      Gfgi <- m * Nfoi * Bo[1] / Bg[1]
      N <- Nfoi + Rv[1] * Gfgi
      G <- Gfgi + Rs[1] * Nfoi
      PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
      BV <- PV / phi
      W <- PV * swi / Bw[1]
      STDERR <- sqrt(sum((F_ - Nfoi * Eowf - Gfgi * Egwf) * ((F_ - Nfoi * Eowf - Gfgi * Egwf))) / (l - 2))
      plot(Et, F_, pch = 21, bg = "blue" , xlab = "Et", ylab = "F")
      abline(a = 0, b = Nfoi)
      mtext(paste("slope = Nfoi = ", round(Nfoi,0), ", m = ", round(m,3), " , STDERR = ", round(STDERR,2)), side = 3)
      xmin <- min(Et, na.rm = TRUE) + 0.25 * diff(range(Et, na.rm = TRUE))
      ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
      text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F - W[e] - N[foi] * E[owf] - G[fgi] * E[gwf])^2), n - 2))))
      p_calc <- vector(length = l)
      for (i in 1:l) {
         if (i == 1) {
            p_calc[i] <- p_[i]
         } else {
            fun <- function(x, pb, p_ini, Bo_ini, Bto_ini, Rs_ini, Bg_ini, Btg_ini, Rv_ini, Bw_ini, Nfoi, Gfgi, swi, p, Ef, Np, Rp, Gp, Wp, Wi, Gi, We, pvt) {
               Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear")$y
               Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear")$y
               Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear")$y
               Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear")$y
               Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear")$y
               Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear")$y
               dp <- p_ini - x
               Bto <- (Bo * (1 - Rs_ini * Rv) + Bg * (Rs_ini - Rs)) / (1 - Rs * Rv)
               Btg <- (Bg * (1 - Rs * Rv_ini) + Bo * (Rv_ini - Rv)) / (1 - Rs * Rv)
               Eo <- Bto - Bto_ini
               Eg <- Btg - Btg_ini
               Ew <- Bw - Bw_ini
               Eowf <- Eo + Bo_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
               Egwf <- Eg + Bg_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
               F_lhs_calc <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
               F_rhs_calc <- Nfoi * Eowf + Gfgi * Egwf
               diff <- F_lhs_calc - F_rhs_calc
               return(diff)
            }
            p_calc[i] <- fzero(fun = fun, x = p_[i], pb = pb, p_ini = p_[1], Bo_ini = Bo[1], Bto_ini = Bto[1], Rs_ini = Rs[1], Bg_ini = Bg[1], Btg_ini = Btg[1], Rv_ini = Rv[1], Bw_ini = Bw[1], Nfoi = Nfoi, Gfgi = Gfgi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Rp = Rp[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], Gi = Gi[i], We = We_[i], pvt = pvt, maxiter = 1000)$x
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

      ymin <- max(0, min((F_[2:l] + We_[2:l]) / Eowf[2:l]) / 1.2)
      ymax <- 1.2 * max((F_[2:l] + We_[2:l]) / Eowf[2:l])
      plot(Np, (F_ + We_) / Eowf, pch = 21, bg = "blue" , xlab = "Np", ylab = "F/Eowf", ylim = c(ymin, ymax))

      mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sorw = sorw, sorg = sorg, p_est = p_calc)
      class(mbal_lst) <- c("combination_oil", "mbal_oil")
      return(mbal_lst)
   }


   if (cls == "combination_We_optim_oil") {
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
      Eowf <- Eo + Bo[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      Egwf <- Eg + Bg[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      F_ <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_
      data <- data.frame(Eowf = Eowf, Egwf = Egwf, `F_` = F_)
      mbal_lm <- lm(`F_` ~ 0 + Eowf + Egwf, data = data)
      Nfoi <- mbal_lm$coefficients[1]
      names(Nfoi) <- NULL
      Gfgi <- mbal_lm$coefficients[2]
      names(Gfgi) <- NULL
      N <- Nfoi + Rv[1] * Gfgi
      G <- Gfgi + Rs[1] * Nfoi
      m <- Gfgi * Bg[1] / (Nfoi * Bo[1])
      PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
      BV <- PV / phi
      W <- PV * swi / Bw[1]
      Et <- Eowf + m * Egwf * (Bo[1] / Bg[1])
      STDERR <- sqrt(sum((F_ - Nfoi * Eowf - Gfgi * Egwf) * ((F_ - Nfoi * Eowf - Gfgi * Egwf))) / (l - 2))
      plot(Et, F_, pch = 21, bg = "blue" , xlab = "Et", ylab = "F")
      abline(a = 0, b = Nfoi)
      mtext(paste("slope = Nfoi = ", round(Nfoi,0), ", m = ", round(m,3), " , STDERR = ", round(STDERR,2)), side = 3)
      xmin <- min(Et, na.rm = TRUE) + 0.25 * diff(range(Et, na.rm = TRUE))
      ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
      text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F - W[e] - N[foi] * E[owf] - G[fgi] * E[gwf])^2), n - 2))))
      p_calc <- vector(length = l)
      for (i in 1:l) {
         if (i == 1) {
            p_calc[i] <- p_[i]
         } else {
            fun <- function(x, pb, p_ini, Bo_ini, Bto_ini, Rs_ini, Bg_ini, Btg_ini, Rv_ini, Bw_ini, Nfoi, Gfgi, swi, p, Ef, Np, Rp, Gp, Wp, Wi, Gi, We, pvt) {
               Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear")$y
               Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear")$y
               Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear")$y
               Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear")$y
               Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear")$y
               Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear")$y
               dp <- p_ini - x
               Bto <- (Bo * (1 - Rs_ini * Rv) + Bg * (Rs_ini - Rs)) / (1 - Rs * Rv)
               Btg <- (Bg * (1 - Rs * Rv_ini) + Bo * (Rv_ini - Rv)) / (1 - Rs * Rv)
               Eo <- Bto - Bto_ini
               Eg <- Btg - Btg_ini
               Ew <- Bw - Bw_ini
               Eowf <- Eo + Bo_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
               Egwf <- Eg + Bg_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
               F_lhs_calc <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
               F_rhs_calc <- Nfoi * Eowf + Gfgi * Egwf
               diff <- F_lhs_calc - F_rhs_calc
               return(diff)
            }
            p_calc[i] <- fzero(fun = fun, x = p_[i], pb = pb, p_ini = p_[1], Bo_ini = Bo[1], Bto_ini = Bto[1], Rs_ini = Rs[1], Bg_ini = Bg[1], Btg_ini = Btg[1], Rv_ini = Rv[1], Bw_ini = Bw[1], Nfoi = Nfoi, Gfgi = Gfgi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Rp = Rp[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], Gi = Gi[i], We = We_[i], pvt = pvt, maxiter = 1000)$x
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

      ymin <- max(0, min((F_[2:l] + We_[2:l]) / Eowf[2:l]) / 1.2)
      ymax <- 1.2 * max((F_[2:l] + We_[2:l]) / Eowf[2:l])
      plot(Np, (F_ + We_) / Eowf, pch = 21, bg = "blue" , xlab = "Np", ylab = "F/Eowf", ylim = c(ymin, ymax))

      mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sorw = sorw, sorg = sorg, p_est = p_calc)
      class(mbal_lst) <- c("combination_oil", "mbal_oil")
      return(mbal_lst)
   }

   if (cls == "combination_N_m_optim_oil") {
      dp <- p_[1] - p_
      Bto <- (Bo * (1 - Rs[1] * Rv) + Bg * (Rs[1] - Rs)) / (1 - Rs * Rv)
      Btg <- (Bg * (1 - Rs * Rv[1]) + Bo * (Rv[1] - Rv)) / (1 - Rs * Rv)
      Eo <- Bto - Bto[1]
      Eg <- Btg - Btg[1]
      Ew <- Bw - Bw[1]
      for (i in 1:l) {
         Ef[i] <- -1 * trapz(p_[1:i], cf_[1:i])
      }
      Eowf <- Eo + Bo[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      Egwf <- Eg + Bg[1] * (swi * Ew / Bw[1] + Ef) / (1 - swi)
      Nfoi <- N / (1 + (m * Bo[1] * Rv[1]) / Bg[1])
      Gfgi <- m * Nfoi * Bo[1] / Bg[1]
      PV <- (Nfoi * Bo[1] + Gfgi * Bg[1]) / (1 - swi)
      BV <- PV / phi
      W <- PV * swi / Bw[1]
      F_ <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw
      We_ <- F_ - Nfoi * Eowf - Gfgi * Egwf
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
         error_fun_ussre_oil <- function(par, We, phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- veh_uss_rad_edge(phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_ussre_oil, We = We_, phi = phi, perm_h = perm_h, h_a = h_a, r_a = r_a, r_R = r_R,
                           tetha = tetha, mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_,
                           pressure = pressure_, lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "uss_rad_bottom") {
         error_fun_ussrb_oil <- function(par, We, phi, perm_h, perm_v, h_a, r_a, r_R, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- yk_uss_rad_bottom(phi, perm_h, perm_v, h_a, r_a, r_R, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_ussrb_oil, We = We_, phi = phi, perm_h = perm_h, perm_v = perm_v, h_a = h_a, r_a = r_a,
                           r_R = r_R, mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_, pressure = pressure_,
                           lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "pss_rad_edge") {
         error_fun_pssre_oil <- function(par, We, phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- fetkovich_pss_rad_edge(phi, perm_h, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_pssre_oil, We = We_, phi = phi, perm_h = perm_h, h_a = h_a, r_a = r_a,
                           r_R = r_R, tetha = tetha, mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_,
                           pressure = pressure_, lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "uss_lin_edge") {
         error_fun_ussle_oil <- function(par, We, phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- nb_uss_lin_edge(phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_ussle_oil, We = We_, phi = phi, perm_h = perm_h, h_a = h_a, w_a = w_a, l_a = l_a,
                           mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_, pressure = pressure_,
                           lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "uss_lin_bottom") {
         error_fun_usslb_oil <- function(par, We, phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- nb_uss_lin_bottom(phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_usslb_oil, We = We_, phi = phi, perm_v = perm_v, h_a = h_a, w_a = w_a, l_a = l_a,
                           mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_, pressure = pressure_,
                           lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "pss_lin_edge") {
         error_fun_pssle_oil <- function(par, We, phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- fetk_pss_lin_edge(phi, perm_h, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_pssle_oil, We = We_, phi = phi, perm_h = perm_h, h_a = h_a, w_a = w_a, l_a = l_a,
                           mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_, pressure = pressure_,
                           lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "pss_lin_bottom") {
         error_fun_psslb_oil <- function(par, We, phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure) {
            We_calc <- fetk_pss_lin_bottom(phi, perm_v, h_a, w_a, l_a, mu_water, c_water, c_rock, time, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_psslb_oil, We = We_, phi = phi, perm_v = perm_v, h_a = h_a, w_a = w_a, l_a = l_a,
                           mu_water = mu_water, c_water = c_water, c_rock = c_rock, time = time_, pressure = pressure_,
                           lower = aquifer$lower, upper = aquifer$upper, control = aquifer$control)
         mult_len_est <- nls.out$par
      }
      if (aqu_cls == "pot") {
         error_fun_pot_oil <- function(par, We, phi, h_a, r_a, r_R, tetha, c_water, c_rock, time, pressure) {
            We_calc <- pot(phi, h_a, r_a, r_R, tetha, c_water, c_rock, pressure, par)
            SSE <- (We_calc - We) ^ 2
            return(SSE)
         }
         nls.out <- nls.lm(par = mult_len, fn = error_fun_pot_oil, We = We_, phi = phi, h_a = h_a, r_a = r_a, r_R = r_R,
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
      F_ <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We_calc
      STDERR <- sqrt(sum((F_ - Nfoi * Eowf) * ((F_ - Nfoi * Eowf))) / (l - 2))
      plot(Nfoi *Eowf + Gfgi * Egwf, F_, pch = 21, bg = "blue" , xlab = "Nfoi * Eowf + Gfgi * Egwf", ylab = "F - We")
      abline(a = 0, b = 1)
      mtext(paste("slope = ", 1.0, " , STDERR = ", round(STDERR,2)), side = 3)
      xmin <- min(Nfoi *Eowf + Gfgi * Egwf, na.rm = TRUE) + 0.25 * diff(range(Nfoi *Eowf + Gfgi * Egwf, na.rm = TRUE))
      ymax <- min(F_, na.rm = TRUE) + 0.75 * diff(range(F_, na.rm = TRUE))
      text(xmin, ymax, expression(STDERR == sqrt(frac(sum((F - W[e] - N[foi] * E[owf] - G[fgi] * E[gwf])^2), n - 2))))
      p_calc <- vector(length = l)
      for (i in 1:l) {
         if (i == 1) {
            p_calc[i] <- p_[i]
         } else {
            fun <- function(x, pb, p_ini, Bo_ini, Bto_ini, Rs_ini, Bg_ini, Btg_ini, Rv_ini, Bw_ini, Nfoi, Gfgi, swi, p, Ef, Np, Rp, Gp, Wp, Wi, Gi, We, pvt) {
               Bo <- approx(pvt$p, pvt$Bo, xout = x, rule = 2, method = "linear")$y
               Rs <- approx(pvt$p, pvt$Rs, xout = x, rule = 2, method = "linear")$y
               Bg <- approx(pvt$p, pvt$Bg, xout = x, rule = 2, method = "linear")$y
               Bw <- approx(pvt$p, pvt$Bw, xout = x, rule = 2, method = "linear")$y
               Rv <- approx(pvt$p, pvt$Rv, xout = x, rule = 2, method = "linear")$y
               Ef_internal <- approx(p, Ef, xout = x, rule = 2, method = "linear")$y
               dp <- p_ini - x
               Bto <- (Bo * (1 - Rs_ini * Rv) + Bg * (Rs_ini - Rs)) / (1 - Rs * Rv)
               Btg <- (Bg * (1 - Rs * Rv_ini) + Bo * (Rv_ini - Rv)) / (1 - Rs * Rv)
               Eo <- Bto - Bto_ini
               Eg <- Btg - Btg_ini
               Ew <- Bw - Bw_ini
               Eowf <- Eo + Bo_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
               Egwf <- Eg + Bg_ini * (swi * Ew / Bw_ini + Ef_internal) / (1 - swi)
               F_lhs_calc <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw - We
               F_rhs_calc <- Nfoi * Eowf + Gfgi * Egwf
               diff <- F_lhs_calc - F_rhs_calc
               return(diff)
            }
            p_calc[i] <- fzero(fun = fun, x = p_[i], pb = pb, p_ini = p_[1], Bo_ini = Bo[1], Bto_ini = Bto[1], Rs_ini = Rs[1], Bg_ini = Bg[1], Btg_ini = Btg[1], Rv_ini = Rv[1], Bw_ini = Bw[1], Nfoi = Nfoi, Gfgi = Gfgi, swi = swi, p = p_, Ef = Ef, Np = Np[i], Rp = Rp[i], Gp = Gp[i], Wp = Wp[i], Wi = Wi[i], Gi = Gi[i], We = We_calc[i], pvt = pvt, maxiter = 1000)$x
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

      ymin <- max(0, min((F_[2:l] + We_calc[2:l]) / Eowf[2:l]) / 1.2)
      ymax <- 1.2 * max((F_[2:l] + We_calc[2:l]) / Eowf[2:l])
      plot(Np, (F_ + We_calc) / Eowf, pch = 21, bg = "blue" , xlab = "Np", ylab = "F/Eowf", ylim = c(ymin, ymax))

      mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, We = We, aquifer = aquifer, wf = wf, sorw = sorw, sorg = sorg, p_est = p_calc)
      class(mbal_lst) <- c("combination_oil", "mbal_oil")
      return(mbal_lst)
   }
}







