
#' A list object of class 'time' for material balance models
#'
#' Create an object of class 'time'
#'
#' @param x a vector of times or a daily sequence of dates
#' @param unit time/date unit of vector x
#'
#' @return a list of class 'time' with all the required parameters for the mbal_perform_oil(), mbal_perform_gas(), mbal_optim_oil(), mbal_optim_gas(), mbal_forecast_oil(), and mbal_forecast_gas() S3 methods
#' @export
#'
#' @examples
#' mbal_time_1 <- mbal_time(c(0:4) * 365, unit = "day")
#'
#' mbal_time_1
#'
#' mbal_time_2 <- mbal_time(c(0:4), unit = "month")
#'
#' mbal_time_2
#'
#' mbal_time_3 <- mbal_time(c(0:4), unit = "year")
#'
#' mbal_time_3
#'
#' mbal_time_4 <- mbal_time(seq(as.Date("2020/1/1"), by = "year",
#' length.out = 5), unit = "date")
#'
#' mbal_time_4

mbal_time <- function(x, unit = "day") {

   if (!is.character(unit)) stop("'unit' must be a character string of type 'day', 'month', 'year', or 'date'.")
   if (!(unit %in% c("day", "month", "year", "date"))) stop("'unit' must be a character string of type 'day', 'month', 'year', or 'date'.")
   if (unit == "date") {
      if (class(x) != "Date") stop("'x' must be a sequence of type 'Date'.")
   }
   if (unit %in% c("day", "month", "year")) {
      if (!is.vector(x)) stop("'x' must be a vector of days, months, or years.")
      if (!is.numeric(x)) stop("'x' must be a numeric vector of days, months, or years.")
      if (any(duplicated(x))) stop("'x' must be a non-duplicated numeric vector of days, months, or years.")
      if (is.unsorted(x)) stop("'x' must be sorted in an ascending order.")
      time_lst <- list(t = x, unit = unit, reference_date = Sys.Date())
      class(time_lst) <- c(unit, "time")
   }
   if (unit == "date") {
      x <- as.Date(x)
      if (any(duplicated(x))) stop("'x' must be a non-duplicated sequence of dates.")
      if (is.unsorted(x)) stop("'x' must be sorted in an ascending order.")
      time_day <- as.numeric(x - x[1])
      time_lst <- list(t = time_day, unit = unit, reference_date = x[1])
      class(time_lst) <- c("day", "time")
   }
   return(time_lst)
}



#***********************mbal_param_oil******************************************

#' A list object of class 'mbal_oil' for material balance analysis
#'
#' Create an object of class 'mbal_oil'
#'
#' @param input_unit a unit system for parameters, only the character string 'Field' is accepted
#' @param output_unit a unit system for properties, only the character string 'Field' is accepted
#' @param aquifer_model defaulted to `NULL`, otherwise must be a character string, one of the following eight options: 'uss_rad_edge', 'uss_rad_bottom', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom', 'pot'. For further information about each model, please see 'Raquifer' package reference manual (https://cran.r-project.org/web/packages/Raquifer/index.html)
#' @param N original oil in place, STB
#' @param m ratio of original gas cap volume to original oil leg volume, a numeric fraction
#' @param phi reservoir porosity, a numeric fraction
#' @param swi initial water saturation in the reservoir, a numeric fraction
#' @param Np cumulative oil production, STB
#' @param Rp ratio of cumulative produced gas to cumulative produced oil
#' @param Wp cumulative water production, STB
#' @param Gi cumulative gas injection, SCF
#' @param Wi cumulative water injection, STB
#' @param We cumulative aquifer water influx, BBL
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
#'
#' @return a list of class ’mbal_oil’ with all the required parameters for the mbal_perform_oil() S3 methods
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
#' mbal_param_oil_lst <- mbal_perform_param_oil(input_unit = "Field", output_unit = "Field",
#' aquifer_model = NULL, N = 1.37e8, m = 0.377, phi = 0.2, swi = 0.2, Np = Np,
#' Rp = Rp, Wp = Wp, Gi = Gi, Wi = Wi, We = We, pb = 3330, p = p, pvt = pvt_table,
#' cf = 0, wf = wf, sorg = 0.2, sorw = 0)
#'
#' dplyr::glimpse(mbal_param_oil_lst)

mbal_perform_param_oil <- function(input_unit = "Field", output_unit = "Field",  aquifer_model = NULL, N = NULL, m = NULL, phi = NULL, swi = NULL, Np = NULL, Rp = NULL, Wp = NULL, Gi = NULL, Wi = NULL, We = NULL, pb = NULL, p = NULL, pvt = NULL, cf = NULL, phi_a = NULL, perm_h_a = NULL, perm_v_a = NULL, h_a = NULL, r_a = NULL, r_R = NULL, w_a = NULL, l_a = NULL, tetha = NULL, muw_a = NULL, cw_a = NULL, cf_a = NULL, wf = NULL, sorg = NULL, sorw = NULL) {

   if (!is.character(input_unit)) stop("'input_unit' must be the character string 'Field'.")
   if (input_unit != "Field") stop("'input_unit' must be the character string 'Field'.")
   if (!is.character(output_unit)) stop("'output_unit' must be the character string 'Field'.")
   if (output_unit != "Field") stop("'output_unit' must be the character string 'Field'.")
   if (!is.null(aquifer_model)) {
      if (!is.character(aquifer_model)) stop("'aquifer_model' must be a character string, one of the following eight options: 'uss_rad_edge', 'uss_rad_bottom', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom', 'pot'.")
      if (!(aquifer_model %in% c('uss_rad_edge', 'uss_rad_bottom', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom', 'pot'))) stop("'aquifer_model' must be one of the following eight options: 'uss_rad_edge', 'uss_rad_bottom', 'uss_lin_edge', 'uss_lin_bottom', 'pss_rad_edge', 'pss_lin_edge', 'pss_lin_bottom', 'pot'.")
      if (length(aquifer_model) > 1) stop("Only one 'aquifer-model' is acceptable for the material balance analysis.")
   }
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
   if (m > 0) {
      if (max(p) != pb) {
         stop("Initial reservoir pressure must be equal to 'pb' in reservoirs with an associated gas cap.")
      }
   }
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
   if (m > 0) {
      if (class(aquifer_lst)[1] == "NoA") {
         if (is.null(sorg)) stop("A numeric value must be assigned to 'sorg' for a 'gas_cap' drive reservoir (m > 0).")
         if (is.null(sorw)) {
            if (any(Wi > 0)) stop("A numeric value must be assigned to 'sorw' for a water injection case.")
         }
      }
      if (class(aquifer_lst)[1] != "NoA") {
         if (is.null(sorw)) stop("A numeric value must be assigned to 'sorw' for a 'combination' drive reservoir.")
         if (is.null(sorg)) stop("A numeric value must be assigned to 'sorg' for a 'combination' drive reservoir.")
      }
   }
   if (m == 0) {
      if (class(aquifer_lst)[1] == "NoA") {
         if (is.null(sorw) & is.null(sorg)) {
            if (any(Gi > 0)) stop("A numeric value must be assigned to 'sorg' for a gas injection case.")
            if (any(Wi > 0)) stop("A numeric value must be assigned to 'sorw' for a water injection case.")
         } else if (is.null(sorw)) {
            if (any(Wi > 0)) stop("A numeric value must be assigned to 'sorw' for a water injection case.")
         } else if (is.null(sorg)) {
            if (any(Gi > 0)) stop("A numeric value must be assigned to 'sorg' for a gas injection case.")
         } else {
         }
      }
      if (class(aquifer_lst)[1] != "NoA") {
         if (is.null(sorw)) stop("A numeric value must be assigned to 'sorw' for a 'water_drive' reservoir.")
         if (is.null(sorg)) {
            if (any(Gi > 0)) stop("A numeric value must be assigned to 'sorg' for a gas injection case.")
         }
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
   pvt <- pvt
   prod <- data.frame(Np = Np, Rp = Rp, Wp = Wp)
   inj <- data.frame(Gi = Gi, Wi = Wi)
   mbal_lst <- list(input_unit = input_unit, output_unit = output_unit, N = N, m = m, phi = phi, swi = swi, pb = pb, p = p, cf = cf, pvt = pvt, prod = prod, inj = inj, aquifer = aquifer_lst, wf = wf, sorw = sorw, sorg = sorg)

   if (m > 0) {
      if (class(aquifer_lst)[1] == "NoA") {
         class(mbal_lst) <- c("gas_cap_oil", "mbal_oil")
      }
      if (class(aquifer_lst)[1] != "NoA") {
         class(mbal_lst) <- c("combination_oil", "mbal_oil")
      }
   }
   if (m == 0) {
      if (class(aquifer_lst)[1] == "NoA") {
         class(mbal_lst) <- c("volumetric_oil", "mbal_oil")
      }
      if (class(aquifer_lst)[1] != "NoA") {
         class(mbal_lst) <- c("water_drive_oil", "mbal_oil")
      }
   }
   return(mbal_lst)
}


# ******************************************************************************

#' Generic function for performance predictions for an oil reservoir
#'
#' Generate a data frame of reservoir performance data according to the class of 'mbal_lst' and 'time_lst' objects
#'
#' @param mbal_lst a list object of class 'mbal_oil'
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
#' mbal_param_oil_lst <- mbal_perform_param_oil(input_unit = "Field", output_unit = "Field",
#' aquifer_model = NULL, N = 1.37e8, m = 0.377, phi = 0.2, swi = 0.2, Np = Np,
#' Rp = Rp, Wp = Wp, Gi = Gi, Wi = Wi, We = We, pb = 3330, p = p, pvt = pvt_table,
#' cf = 0, wf = wf, sorg = 0.2, sorw = 0)
#'
#' time_lst <- mbal_time(c(0, 365, 730, 1095, 1460, 1825, 2190), "day")
#'
#' mbal_results <- mbal_perform_oil(mbal_param_oil_lst, time_lst)
#'
#' dplyr::glimpse(mbal_results)

mbal_perform_oil <- function(mbal_lst, time_lst) {

   if (inherits(mbal_lst, "mbal_oil") == TRUE & inherits(time_lst, "time")) {
      UseMethod("mbal_perform_oil")
   } else {
      if (!inherits(mbal_lst, "mbal_oil")) {
         stop("A class of 'mbal_oil' must be assigned to the 'mbal_lst' parameter of the mbal_perform_oil() function.")
      }
      if (!inherits(time_lst, "time")) {
         stop("A class of 'time' must be assigned to the 'time_lst' parameter of the mbal_perform_oil() function.")
      }
   }
}



# ******************************************************************************

#' S3 method for class 'mbal_perform_oil'
#'
#' Return a data frame with estimates for fluids saturation, drive indices, production rates, and gas-oil ratios over the pressure history of a volumetric oil reservoir
#'
#' @param mbal_lst a list object of class 'mbal_oil'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame with estimates for fluids saturation, drive indices, production rates, and gas-oil ratios over the pressure history of a volumetric oil reservoir
#' @export
#'
mbal_perform_oil.volumetric_oil <- function(mbal_lst, time_lst) {
   volumetric_oil(mbal_lst, time_lst)
}


# ******************************************************************************

#' S3 method for class 'mbal_perform_oil'
#'
#' Return a data frame with estimates for fluids saturation, drive indices, production rates, and gas-oil ratios over the pressure history of a gas_cap_drive oil reservoir
#'
#' @param mbal_lst a list object of class 'mbal_oil'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame with estimates for fluids saturation, drive indices, production rates, and gas-oil ratios over the pressure history of a gas_cap_drive oil reservoir
#' @export
mbal_perform_oil.gas_cap_oil <- function(mbal_lst, time_lst) {
   gas_cap_oil(mbal_lst, time_lst)
}


# ******************************************************************************

#' S3 method for class 'mbal_perform_oil'
#'
#' Return a data frame with estimates for fluids saturation, drive indices, production rates, and gas-oil ratios over the pressure history of a water_drive oil reservoir
#'
#' @param mbal_lst a list object of class 'mbal_oil'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame with estimates for fluids saturation, drive indices, production rates, and gas-oil ratios over the pressure history of a water_drive oil reservoir
#' @export
mbal_perform_oil.water_drive_oil <- function(mbal_lst, time_lst) {
   water_drive_oil(mbal_lst, time_lst)
}


# ******************************************************************************

#' S3 method for class 'mbal_perform_oil'
#'
#' Return a data frame with estimates for fluids saturation, drive indices, production rates, and gas-oil ratios over the pressure history of a combination_drive oil reservoir
#'
#' @param mbal_lst a list object of class 'mbal_oil'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame with estimates for fluids saturation, drive indices, production rates, and gas-oil ratios over the pressure history of a combination_drive oil reservoir
#' @export
mbal_perform_oil.combination_oil <- function(mbal_lst, time_lst) {
   combination_oil(mbal_lst, time_lst)
}


# ******************************************************************************


volumetric_oil <- function(mbal_lst, time_lst) {

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
   N <- mbal_lst$N
   m <- mbal_lst$m
   phi <- mbal_lst$phi
   swi <- mbal_lst$swi
   pb <- mbal_lst$pb
   p <- mbal_lst$p
   l <- length(p)
   if (length(time_lst$t) != length(p)) stop("Lenght of 'time_lst$t' vector must be equal to the length of 'p' vector.")
   cf <- mbal_lst$cf
   pvt <- mbal_lst$pvt
   prod <- mbal_lst$prod
   inj <- mbal_lst$inj
   aquifer <- mbal_lst$aquifer
   aqu_cls <- class(aquifer)[1]
   if (aqu_cls == "NoA") We <- aquifer$We
   wf <- mbal_lst$wf
   sorw <- mbal_lst$sorw
   sorg <- mbal_lst$sorg
   if (is.null(sorw)) stop("'sorw' must be a numeric value.")
   if (is.null(sorg)) stop("'sorg' must be a numeric value.")
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
   Rp <- prod$Rp
   Wp <- prod$Wp
   Gi <- inj$Gi
   Wi <- inj$Wi
   Gp <- Np * Rp
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
   PVgas <- vector(length = l)
   PVwater <- vector(length = l)
   sw <- vector(length = l)
   so <- vector(length = l)
   sg <- vector(length = l)
   sw_t <- vector(length = l)
   so_t <- vector(length = l)
   sg_t <- vector(length = l)
   Nfoi <- N / (1 + (m * Bo[1] * Rv[1]) / Bg[1])
   Gfgi <- m * Nfoi * Bo[1] / Bg[1]
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
   Et <- Gfgi * Eg + Nfoi * Eo + W * Ew + PV * Ef + We + (Wi - Wp) * Bw
   F_ <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw
   RF_oil <- (((Wi - Wp) * Bw + We) / Nfoi + Eowf + (m * Bo[1] * Egwf) / Bg[1]) / ((1 + m * Bo[1] * Rv[1] / Bg[1]) * ((Bo * (1 - Rv * Rp) +  Bg * (Rp - Rs)) / (1 - Rv * Rs)))
   RF_gas <- Rp * (RF_oil) * ((Bg[1] + m * Bo[1] * Rv[1]) / (m * Bo[1] + Bg[1] * Rs[1]))
   Igd <- ifelse(dplyr::row_number(Igd) == 1, NA, Gfgi * Eg / Et)
   Isd <- ifelse(dplyr::row_number(Isd) == 1, NA, Nfoi * Eo / Et)
   Inwd <- ifelse(dplyr::row_number(Inwd) == 1, NA, We / Et)
   Ifwd <- ifelse(dplyr::row_number(Ifwd) == 1, NA, (W * Ew + PV * Ef) / Et)
   Iawd <- ifelse(dplyr::row_number(Iawd) == 1, NA, (Wi - Wp) * Bw / Et)
   Itot <- Igd + Isd + Inwd + Ifwd + Iawd
   soi <- 1- swi
   sgi <- 1- swi - soi
   PVgasi <- (Gfgi * Bg[1] * soi - Nfoi * Bo[1] * sgi) / ((1 - sorg - swi) * soi - sorg * sgi)
   PVoili <- PV - PVgasi
   xp <- 1
   PVgas_r <- PVgasi / PV + PVgasi / PV * xp * ((sorg * Rs[1] / Bo[1] + (1 - sorg - swi) / Bg[1]) / (sorg * Rs[1] / Bo + (1 - sorg - swi) / Bg)- 1) + Gi / ((sorg * Rs / Bo) + ((1 - sorg - swi) / Bg)) / PV
   PVwater_r <- ((Wi - Wp) * Bw + We) / (1 - sorw - swi) / PV
   swi_t <- swi
   soi_t <- soi
   sgi_t <- 1 - swi_t - soi_t
   for (i in 1:l) {
      if (p[i] >= pb) {
         sw[i] <- swi
         so[i] <- 1 - sw[i]
         sg[i] <- 1 - sw[i] - so[i]
         sw_t[i] <- sw[i] * (1 - PVwater_r[i]) + (1 - sorw) * PVwater_r[i]
         so_t[i] <- so[i] * (1 - PVgas_r[i] - PVwater_r[i]) + sorg * PVgas_r[i] + sorw * PVwater_r[i]
         sg_t[i] <- 1 - so_t[i] - sw_t[i]
      } else {
         sw_t[i] <- swi * (1 - PVwater_r[i]) + (1 - sorw) * PVwater_r[i]
         so_t[i] <- ((1 - Np[i] / N) * (soi_t / Bo[1] + sgi_t * Rv[1] / Bg[1]) - (1 - sw_t[i]) * Rv[i] / Bg[i]) / (1 / Bo[i] - Rv[i] / Bg[i])
         sg_t[i] <- 1 - so_t[i] - sw_t[i]
         sw[i] <- swi
         so[i] <- (so_t[i] - sorg * PVgas_r[i] - sorw * PVwater_r[i]) / (1 - PVgas_r[i] - PVwater_r[i])
         sg[i] <- 1 - sw[i] - so[i]
      }
   }
   qo <- vector(length = l)
   qw <- vector(length = l)
   qg <- vector(length = l)
   qg_insitu <- vector(length = l)
   fw <- vector(length = l)
   fg <- vector(length = l)
   gor <- vector(length = l)
   krg_kro <- vector(length = l)
   qo <- ifelse(dplyr::row_number(Np) == 1, 0, (Np - dplyr::lag(Np)) / (time - dplyr::lag(time)))
   qg <- ifelse(dplyr::row_number(Np) == 1, 0, (Np * Rp - dplyr::lag(Np) * dplyr::lag(Rp)) / (time - dplyr::lag(time)))
   qg_insitu <- ifelse(dplyr::row_number(Np) == 1, 0, (Np * (Rp - Rs) - dplyr::lag(Np) * (dplyr::lag(Rp) - dplyr::lag(Rs))) / (time - dplyr::lag(time)))
   qw <- ifelse(dplyr::row_number(Wp) == 1, 0, (Wp - dplyr::lag(Wp)) / (time - dplyr::lag(time)))
   fg <- ifelse(dplyr::row_number(qg_insitu) == 1, 0, qg_insitu * Bg / (qo * Bo + qg_insitu * Bg + qw * Bw))
   fw <- ifelse(dplyr::row_number(qw) == 1, 0, qw * Bw / (qo * Bo + qg_insitu * Bg + qw * Bw))
   gor <- ifelse(p >= pb, Rs, qg / qo)
   krg_kro <- ifelse(p >= pb, 0, ((Rs - gor) / muo / Bo) / ((Rv * gor - 1) / mug / Bg))
   krg_kro <- ifelse(krg_kro < 0, 0, krg_kro)

   names <- c("P (psia)", "Eo (bbl/STB)", "Eg (bbl/SCF)", "Ew (bbl/STB)", "Ef (bbl/bbl)", "Eowf (bbl/STB)", "Egwf (bbl/SCF)", "Et (bbl)", "F (bbl)", "We", "Igd", "Isd", "Inwd", "Ifwd", "Iawd", "Itot", "RF_oil", "RF_gas", "SOo", "SGo", "SWo", "SOT", "SGT", "SWT", "qo (STB/day)", "qg (SCF/day)", "qw (STB/day)", "fg", "fw", "GOR (SCF/STB)", "krg/kro")
   results <- data.frame(p = p, Eo = Eo, Eg = Eg, Ew = Ew, Ef = Ef, Eowf = Eowf, Egwf = Egwf, Et = Et, F_ = F_, We = We, Igd = Igd, Isd = Isd, Inwd = Inwd, Ifwd = Ifwd, Iawd = Iawd, Itot = Itot, RF_oil = RF_oil, RF_gas = RF_gas, SO = so, SG = sg, SW = sw, SOT = so_t, SGT = sg_t, SWT = sw_t, qo = qo, qg = qg, qw = qw, fg = fg, fw = fw, gor = gor, krg_kro = krg_kro)
   colnames(results) <- names
   return(results)
}




gas_cap_oil <- function(mbal_lst, time_lst) {

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
   N <- mbal_lst$N
   m <- mbal_lst$m
   phi <- mbal_lst$phi
   swi <- mbal_lst$swi
   pb <- mbal_lst$pb
   p <- mbal_lst$p
   l <- length(p)
   if (length(time_lst$t) != length(p)) stop("Lenght of 'time_lst$t' vector must be equal to the length of 'p' vector.")
   cf <- mbal_lst$cf
   pvt <- mbal_lst$pvt
   prod <- mbal_lst$prod
   inj <- mbal_lst$inj
   aquifer <- mbal_lst$aquifer
   aqu_cls <- class(aquifer)[1]
   if (aqu_cls == "NoA") We <- aquifer$We
   wf <- mbal_lst$wf
   sorw <- mbal_lst$sorw
   sorg <- mbal_lst$sorg
   if (is.null(sorw)) stop("'sorw' must be a numeric value.")
   if (is.null(sorg)) stop("'sorg' must be a numeric value.")
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
   Rp <- prod$Rp
   Wp <- prod$Wp
   Gi <- inj$Gi
   Wi <- inj$Wi
   Gp <- Np * Rp
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
   sw <- vector(length = l)
   so <- vector(length = l)
   sg <- vector(length = l)
   sw_t <- vector(length = l)
   so_t <- vector(length = l)
   sg_t <- vector(length = l)
   Nfoi <- N / (1 + (m * Bo[1] * Rv[1]) / Bg[1])
   Gfgi <- m * Nfoi * Bo[1] / Bg[1]
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
   Et <- Gfgi * Eg + Nfoi * Eo + W * Ew + PV * Ef + We + (Wi - Wp) * Bw
   F_ <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw
   RF_oil <- (((Wi - Wp) * Bw + We) / Nfoi + Eowf + (m * Bo[1] * Egwf) / Bg[1]) / ((1 + m * Bo[1] * Rv[1] / Bg[1]) * ((Bo * (1 - Rv * Rp) +  Bg * (Rp - Rs)) / (1 - Rv * Rs)))
   RF_gas <- Rp * (RF_oil) * ((Bg[1] + m * Bo[1] * Rv[1]) / (m * Bo[1] + Bg[1] * Rs[1]))
   Igd <- ifelse(dplyr::row_number(Igd) == 1, NA, Gfgi * Eg / Et)
   Isd <- ifelse(dplyr::row_number(Isd) == 1, NA, Nfoi * Eo / Et)
   Inwd <- ifelse(dplyr::row_number(Inwd) == 1, NA, We / Et)
   Ifwd <- ifelse(dplyr::row_number(Ifwd) == 1, NA, (W * Ew + PV * Ef) / Et)
   Iawd <- ifelse(dplyr::row_number(Iawd) == 1, NA, (Wi - Wp) * Bw / Et)
   Itot <- Igd + Isd + Inwd + Ifwd + Iawd
   soi <- 1 - swi
   sgi <- 1 - swi - soi
   PVgasi <- (Gfgi * Bg[1] * soi - Nfoi * Bo[1] * sgi) / ((1 - sorg - swi) * soi - sorg * sgi)
   PVoili <- PV - PVgasi
   xp <- 1
   PVgas_r <- PVgasi / PV + PVgasi / PV * xp * ((sorg * Rs[1] / Bo[1] + (1 - sorg - swi) / Bg[1]) / (sorg * Rs[1] / Bo + (1 - sorg - swi) / Bg)- 1) + Gi / ((sorg * Rs / Bo) + ((1 - sorg - swi) / Bg)) / PV
   PVwater_r <- ((Wi - Wp) * Bw + We) / (1 - sorw - swi) / PV
   swi_t <- swi
   soi_t <- soi * (1 - PVgas_r[1]) + sorg * PVgas_r[1]
   sgi_t <- 1 - swi_t - soi_t
   for (i in 1:l) {
      if (p[i] >= pb) {
         sw[i] <- swi
         so[i] <- 1 - sw[i]
         sg[i] <- 1 - sw[i] - so[i]
         sw_t[i] <- sw[i] * (1 - PVwater_r[i]) + (1 - sorw) * PVwater_r[i]
         so_t[i] <- so[i] * (1 - PVgas_r[i] - PVwater_r[i]) + sorg * PVgas_r[i] + sorw * PVwater_r[i]
         sg_t[i] <- 1 - so_t[i] - sw_t[i]
      } else {
         sw_t[i] <- swi * (1 - PVwater_r[i]) + (1 - sorw) * PVwater_r[i]
         so_t[i] <- ((1 - Np[i] / N) * (soi_t / Bo[1] + sgi_t * Rv[1] / Bg[1]) - (1 - sw_t[i]) * Rv[i] / Bg[i]) / (1 / Bo[i] - Rv[i] / Bg[i])
         sg_t[i] <- 1 - so_t[i] - sw_t[i]
         sw[i] <- swi
         so[i] <- (so_t[i] - sorg * PVgas_r[i] - sorw * PVwater_r[i]) / (1 - PVgas_r[i] - PVwater_r[i])
         sg[i] <- 1 - sw[i] - so[i]
      }
   }
   qo <- vector(length = l)
   qw <- vector(length = l)
   qg <- vector(length = l)
   qg_insitu <- vector(length = l)
   fw <- vector(length = l)
   fg <- vector(length = l)
   gor <- vector(length = l)
   krg_kro <- vector(length = l)
   qo <- ifelse(dplyr::row_number(Np) == 1, 0, (Np - dplyr::lag(Np)) / (time - dplyr::lag(time)))
   qg <- ifelse(dplyr::row_number(Np) == 1, 0, (Np * Rp - dplyr::lag(Np) * dplyr::lag(Rp)) / (time - dplyr::lag(time)))
   qg_insitu <- ifelse(dplyr::row_number(Np) == 1, 0, (Np * (Rp - Rs) - dplyr::lag(Np) * (dplyr::lag(Rp) - dplyr::lag(Rs))) / (time - dplyr::lag(time)))
   qw <- ifelse(dplyr::row_number(Wp) == 1, 0, (Wp - dplyr::lag(Wp)) / (time - dplyr::lag(time)))
   fg <- ifelse(dplyr::row_number(qg_insitu) == 1, 0, qg_insitu * Bg / (qo * Bo + qg_insitu * Bg + qw * Bw))
   fw <- ifelse(dplyr::row_number(qw) == 1, 0, qw * Bw / (qo * Bo + qg_insitu * Bg + qw * Bw))
   gor <- ifelse(p >= pb, Rs, qg / qo)
   krg_kro <- ifelse(p >= pb, 0, ((Rs - gor) / muo / Bo) / ((Rv * gor - 1) / mug / Bg))
   krg_kro <- ifelse(krg_kro < 0, 0, krg_kro)
   names <- c("P (psia)", "Eo (bbl/STB)", "Eg (bbl/SCF)", "Ew (bbl/STB)", "Ef (bbl/bbl)", "Eowf (bbl/STB)", "Egwf (bbl/SCF)", "Et (bbl)", "F (bbl)", "We", "Igd", "Isd", "Inwd", "Ifwd", "Iawd", "Itot", "RF_oil", "RF_gas", "SOo", "SGo", "SWo", "SOT", "SGT", "SWT", "qo (STB/day)", "qg (SCF/day)", "qw (STB/day)", "fg", "fw", "GOR (SCF/STB)", "krg/kro")
   results <- data.frame(p = p, Eo = Eo, Eg = Eg, Ew = Ew, Ef = Ef, Eowf = Eowf, Egwf = Egwf, Et = Et, F_ = F_, We = We, Igd = Igd, Isd = Isd, Inwd = Inwd, Ifwd = Ifwd, Iawd = Iawd, Itot = Itot, RF_oil = RF_oil, RF_gas = RF_gas, SO = so, SG = sg, SW = sw, SOT = so_t, SGT = sg_t, SWT = sw_t, qo = qo, qg = qg, qw = qw, fg = fg, fw = fw, gor = gor, krg_kro = krg_kro)
   colnames(results) <- names
   return(results)
}



water_drive_oil <- function(mbal_lst, time_lst) {

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
   N <- mbal_lst$N
   m <- mbal_lst$m
   phi <- mbal_lst$phi
   swi <- mbal_lst$swi
   pb <- mbal_lst$pb
   p <- mbal_lst$p
   l <- length(p)
   if (length(time_lst$t) != length(p)) stop("Lenght of 'time_lst$t' vector must be equal to the length of 'p' vector.")
   cf <- mbal_lst$cf
   pvt <- mbal_lst$pvt
   prod <- mbal_lst$prod
   inj <- mbal_lst$inj
   aquifer <- mbal_lst$aquifer
   aqu_cls <- class(aquifer)[1]
   if (aqu_cls == "We") We <- aquifer$We
   wf <- mbal_lst$wf
   sorw <- mbal_lst$sorw
   sorg <- mbal_lst$sorg
   if (is.null(sorw)) stop("'sorw' must be a numeric value.")
   if (is.null(sorg)) stop("'sorg' must be a numeric value.")
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
   Rp <- prod$Rp
   Wp <- prod$Wp
   Gi <- inj$Gi
   Wi <- inj$Wi
   Gp <- Np * Rp
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
   Isd <- vector(length = l)
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
   Nfoi <- N / (1 + (m * Bo[1] * Rv[1]) / Bg[1])
   Gfgi <- m * Nfoi * Bo[1] / Bg[1]
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
   Et <- Gfgi * Eg + Nfoi * Eo + W * Ew + PV * Ef + We + (Wi - Wp) * Bw
   F_ <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw
   RF_oil <- (((Wi - Wp) * Bw + We) / Nfoi + Eowf + (m * Bo[1] * Egwf) / Bg[1]) / ((1 + m * Bo[1] * Rv[1] / Bg[1]) * ((Bo * (1 - Rv * Rp) +  Bg * (Rp - Rs)) / (1 - Rv * Rs)))
   RF_gas <- Rp * (RF_oil) * ((Bg[1] + m * Bo[1] * Rv[1]) / (m * Bo[1] + Bg[1] * Rs[1]))
   Igd <- ifelse(dplyr::row_number(Igd) == 1, NA, Gfgi * Eg / Et)
   Isd <- ifelse(dplyr::row_number(Isd) == 1, NA, Nfoi * Eo / Et)
   Inwd <- ifelse(dplyr::row_number(Inwd) == 1, NA, We / Et)
   Ifwd <- ifelse(dplyr::row_number(Ifwd) == 1, NA, (W * Ew + PV * Ef) / Et)
   Iawd <- ifelse(dplyr::row_number(Iawd) == 1, NA, (Wi - Wp) * Bw / Et)
   Itot <- Igd + Isd + Inwd + Ifwd + Iawd
   soi <- 1- swi
   sgi <- 1- swi - soi
   PVgasi <- (Gfgi * Bg[1] * soi - Nfoi * Bo[1] * sgi) / ((1 - sorg - swi) * soi - sorg * sgi)
   PVoili <- PV - PVgasi
   xp <- 1
   PVgas_r <- PVgasi / PV + PVgasi / PV * xp * ((sorg * Rs[1] / Bo[1] + (1 - sorg - swi) / Bg[1]) / (sorg * Rs[1] / Bo + (1 - sorg - swi) / Bg)- 1) + Gi / ((sorg * Rs / Bo) + ((1 - sorg - swi) / Bg)) / PV
   PVwater_r <- ((Wi - Wp) * Bw + We) / (1 - sorw - swi) / PV
   swi_t <- swi
   soi_t <- soi
   sgi_t <- 1 - swi_t - soi_t
   for (i in 1:l) {
      if (p[i] >= pb) {
         sw[i] <- swi
         so[i] <- 1 - sw[i]
         sg[i] <- 1 - sw[i] - so[i]
         sw_t[i] <- sw[i] * (1 - PVwater_r[i]) + (1 - sorw) * PVwater_r[i]
         so_t[i] <- so[i] * (1 - PVgas_r[i] - PVwater_r[i]) + sorg * PVgas_r[i] + sorw * PVwater_r[i]
         sg_t[i] <- 1 - so_t[i] - sw_t[i]
      } else {
         sw_t[i] <- swi * (1 - PVwater_r[i]) + (1 - sorw) * PVwater_r[i]
         so_t[i] <- ((1 - Np[i] / N) * (soi_t / Bo[1] + sgi_t * Rv[1] / Bg[1]) - (1 - sw_t[i]) * Rv[i] / Bg[i]) / (1 / Bo[i] - Rv[i] / Bg[i])
         sg_t[i] <- 1 - so_t[i] - sw_t[i]
         sw[i] <- swi
         so[i] <- (so_t[i] - sorg * PVgas_r[i] - sorw * PVwater_r[i]) / (1 - PVgas_r[i] - PVwater_r[i])
         sg[i] <- 1 - sw[i] - so[i]
      }
   }
   qo <- vector(length = l)
   qw <- vector(length = l)
   qg <- vector(length = l)
   qg_insitu <- vector(length = l)
   fw <- vector(length = l)
   fg <- vector(length = l)
   gor <- vector(length = l)
   krg_kro <- vector(length = l)
   qo <- ifelse(dplyr::row_number(Np) == 1, 0, (Np - dplyr::lag(Np)) / (time - dplyr::lag(time)))
   qg <- ifelse(dplyr::row_number(Np) == 1, 0, (Np * Rp - dplyr::lag(Np) * dplyr::lag(Rp)) / (time - dplyr::lag(time)))
   qg_insitu <- ifelse(dplyr::row_number(Np) == 1, 0, (Np * (Rp - Rs) - dplyr::lag(Np) * (dplyr::lag(Rp) - dplyr::lag(Rs))) / (time - dplyr::lag(time)))
   qw <- ifelse(dplyr::row_number(Wp) == 1, 0, (Wp - dplyr::lag(Wp)) / (time - dplyr::lag(time)))
   fg <- ifelse(dplyr::row_number(qg_insitu) == 1, 0, qg_insitu * Bg / (qo * Bo + qg_insitu * Bg + qw * Bw))
   fw <- ifelse(dplyr::row_number(qw) == 1, 0, qw * Bw / (qo * Bo + qg_insitu * Bg + qw * Bw))
   gor <- ifelse(p >= pb, Rs, qg / qo)
   krg_kro <- ifelse(p >= pb, 0, ((Rs - gor) / muo / Bo) / ((Rv * gor - 1) / mug / Bg))
   krg_kro <- ifelse(krg_kro < 0, 0, krg_kro)
   names <- c("P (psia)", "Eo (bbl/STB)", "Eg (bbl/SCF)", "Ew (bbl/STB)", "Ef (bbl/bbl)", "Eowf (bbl/STB)", "Egwf (bbl/SCF)", "Et (bbl)", "F (bbl)", "We", "Igd", "Isd", "Inwd", "Ifwd", "Iawd", "Itot", "RF_oil", "RF_gas", "SOo", "SGo", "SWo", "SOT", "SGT", "SWT", "qo (STB/day)", "qg (SCF/day)", "qw (STB/day)", "fg", "fw", "GOR (SCF/STB)", "krg/kro")
   results <- data.frame(p = p, Eo = Eo, Eg = Eg, Ew = Ew, Ef = Ef, Eowf = Eowf, Egwf = Egwf, Et = Et, F_ = F_, We = We, Igd = Igd, Isd = Isd, Inwd = Inwd, Ifwd = Ifwd, Iawd = Iawd, Itot = Itot, RF_oil = RF_oil, RF_gas = RF_gas, SO = so, SG = sg, SW = sw, SOT = so_t, SGT = sg_t, SWT = sw_t, qo = qo, qg = qg, qw = qw, fg = fg, fw = fw, gor = gor, krg_kro = krg_kro)
   colnames(results) <- names
   return(results)
}





combination_oil <- function(mbal_lst, time_lst) {

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
   N <- mbal_lst$N
   m <- mbal_lst$m
   phi <- mbal_lst$phi
   swi <- mbal_lst$swi
   pb <- mbal_lst$pb
   p <- mbal_lst$p
   l <- length(p)
   if (length(time_lst$t) != length(p)) stop("Lenght of 'time_lst$t' vector must be equal to the length of 'p' vector.")
   cf <- mbal_lst$cf
   pvt <- mbal_lst$pvt
   prod <- mbal_lst$prod
   inj <- mbal_lst$inj
   aquifer <- mbal_lst$aquifer
   aqu_cls <- class(aquifer)[1]
   if (aqu_cls == "We") We <- aquifer$We
   wf <- mbal_lst$wf
   sorw <- mbal_lst$sorw
   sorg <- mbal_lst$sorg
   if (is.null(sorw)) stop("'sorw' must be a numeric value.")
   if (is.null(sorg)) stop("'sorg' must be a numeric value.")
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
   Rp <- prod$Rp
   Wp <- prod$Wp
   Gi <- inj$Gi
   Wi <- inj$Wi
   Gp <- Np * Rp
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
   Isd <- vector(length = l)
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
   Nfoi <- N / (1 + (m * Bo[1] * Rv[1]) / Bg[1])
   Gfgi <- m * Nfoi * Bo[1] / Bg[1]
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
   Et <- Gfgi * Eg + Nfoi * Eo + W * Ew + PV * Ef + We + (Wi - Wp) * Bw
   F_ <- (Gp - Gi) * ((Bg - Bo * Rv) / (1 - Rv * Rs)) + Np * ((Bo - Bg * Rs) / (1 - Rv * Rs)) + (Wp - Wi) * Bw
   RF_oil <- (((Wi - Wp) * Bw + We) / Nfoi + Eowf + (m * Bo[1] * Egwf) / Bg[1]) / ((1 + m * Bo[1] * Rv[1] / Bg[1]) * ((Bo * (1 - Rv * Rp) +  Bg * (Rp - Rs)) / (1 - Rv * Rs)))
   RF_gas <- Rp * (RF_oil) * ((Bg[1] + m * Bo[1] * Rv[1]) / (m * Bo[1] + Bg[1] * Rs[1]))
   Igd <- ifelse(dplyr::row_number(Igd) == 1, NA, Gfgi * Eg / Et)
   Isd <- ifelse(dplyr::row_number(Isd) == 1, NA, Nfoi * Eo / Et)
   Inwd <- ifelse(dplyr::row_number(Inwd) == 1, NA, We / Et)
   Ifwd <- ifelse(dplyr::row_number(Ifwd) == 1, NA, (W * Ew + PV * Ef) / Et)
   Iawd <- ifelse(dplyr::row_number(Iawd) == 1, NA, (Wi - Wp) * Bw / Et)
   Itot <- Igd + Isd + Inwd + Ifwd + Iawd
   soi <- 1 - swi
   sgi <- 1 - swi - soi
   PVgasi <- (Gfgi * Bg[1] * soi - Nfoi * Bo[1] * sgi) / ((1 - sorg - swi) * soi - sorg * sgi)
   PVoili <- PV - PVgasi
   xp <- 1
   PVgas_r <- PVgasi / PV + PVgasi / PV * xp * ((sorg * Rs[1] / Bo[1] + (1 - sorg - swi) / Bg[1]) / (sorg * Rs[1] / Bo + (1 - sorg - swi) / Bg)- 1) + Gi / ((sorg * Rs / Bo) + ((1 - sorg - swi) / Bg)) / PV
   PVwater_r <- ((Wi - Wp) * Bw + We) / (1 - sorw - swi) / PV


   swi_t <- swi
   soi_t <- soi * soi * (1 - PVgas_r[1]) + sorg * PVgas_r[1]
   sgi_t <- 1 - swi_t - soi_t
   for (i in 1:l) {
      if (p[i] >= pb) {
         sw[i] <- swi
         so[i] <- 1 - sw[i]
         sg[i] <- 1 - sw[i] - so[i]
         sw_t[i] <- sw[i] * (1 - PVwater_r[i]) + (1 - sorw) * PVwater_r[i]
         so_t[i] <- so[i] * (1 - PVgas_r[i] - PVwater_r[i]) + sorg * PVgas_r[i] + sorw * PVwater_r[i]
         sg_t[i] <- 1 - so_t[i] - sw_t[i]
      } else {
         sw_t[i] <- swi * (1 - PVwater_r[i]) + (1 - sorw) * PVwater_r[i]
         so_t[i] <- ((1 - Np[i] / N) * (soi_t / Bo[1] + sgi_t * Rv[1] / Bg[1]) - (1 - sw_t[i]) * Rv[i] / Bg[i]) / (1 / Bo[i] - Rv[i] / Bg[i])
         sg_t[i] <- 1 - so_t[i] - sw_t[i]
         sw[i] <- swi
         so[i] <- (so_t[i] - sorg * PVgas_r[i] - sorw * PVwater_r[i]) / (1 - PVgas_r[i] - PVwater_r[i])
         sg[i] <- 1 - sw[i] - so[i]
      }
   }
   qo <- vector(length = l)
   qw <- vector(length = l)
   qg <- vector(length = l)
   qg_insitu <- vector(length = l)
   fw <- vector(length = l)
   fg <- vector(length = l)
   gor <- vector(length = l)
   krg_kro <- vector(length = l)
   qo <- ifelse(dplyr::row_number(Np) == 1, 0, (Np - dplyr::lag(Np)) / (time - dplyr::lag(time)))
   qg <- ifelse(dplyr::row_number(Np) == 1, 0, (Np * Rp - dplyr::lag(Np) * dplyr::lag(Rp)) / (time - dplyr::lag(time)))
   qg_insitu <- ifelse(dplyr::row_number(Np) == 1, 0, (Np * (Rp - Rs) - dplyr::lag(Np) * (dplyr::lag(Rp) - dplyr::lag(Rs))) / (time - dplyr::lag(time)))
   qw <- ifelse(dplyr::row_number(Wp) == 1, 0, (Wp - dplyr::lag(Wp)) / (time - dplyr::lag(time)))
   fg <- ifelse(dplyr::row_number(qg_insitu) == 1, 0, qg_insitu * Bg / (qo * Bo + qg_insitu * Bg + qw * Bw))
   fw <- ifelse(dplyr::row_number(qw) == 1, 0, qw * Bw / (qo * Bo + qg_insitu * Bg + qw * Bw))
   gor <- ifelse(p >= pb, Rs, qg / qo)
   krg_kro <- ifelse(p >= pb, 0, ((Rs - gor) / muo / Bo) / ((Rv * gor - 1) / mug / Bg))
   krg_kro <- ifelse(krg_kro < 0, 0, krg_kro)
   names <- c("P (psia)", "Eo (bbl/STB)", "Eg (bbl/SCF)", "Ew (bbl/STB)", "Ef (bbl/bbl)", "Eowf (bbl/STB)", "Egwf (bbl/SCF)", "Et (bbl)", "F (bbl)", "We", "Igd", "Isd", "Inwd", "Ifwd", "Iawd", "Itot", "RF_oil", "RF_gas", "SOo", "SGo", "SWo", "SOT", "SGT", "SWT", "qo (STB/day)", "qg (SCF/day)", "qw (STB/day)", "fg", "fw", "GOR (SCF/STB)", "krg/kro")
   results <- data.frame(p = p, Eo = Eo, Eg = Eg, Ew = Ew, Ef = Ef, Eowf = Eowf, Egwf = Egwf, Et = Et, F_ = F_, We = We, Igd = Igd, Isd = Isd, Inwd = Inwd, Ifwd = Ifwd, Iawd = Iawd, Itot = Itot, RF_oil = RF_oil, RF_gas = RF_gas, SO = so, SG = sg, SW = sw, SOT = so_t, SGT = sg_t, SWT = sw_t, qo = qo, qg = qg, qw = qw, fg = fg, fw = fw, gor = gor, krg_kro = krg_kro)
   colnames(results) <- names
   return(results)
}


