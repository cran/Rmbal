
#' @importFrom pracma fzero
#' @importFrom stats approx lag lm

#****************************************************************************************************
# Stehfest Algorithm

stehfest <- function(n) {
   # Coefficients for the Gaver-Stehfest algorithm
   v <- vector(length = n)
   m <- 0.5 * n
   for (i in 1:n) {
      sum <- 0
      j <- floor(0.5 * (i + 1))
      l <- min(i, m)
      for (k in j:l) {
         sum <- sum + (k ^ m * factorial(2 * k)) / (factorial(m - k) * factorial(k) * factorial(k - 1) * factorial(i - k) * factorial(2 * k - i))
      }
      v[i] = (-1) ^ (i + m) * sum
   }
   return(v)
}


#****************************************************************************************************
# Van Everdingen and Hurst Unsteady State Model - Radial Flow - Edge Water - Laplace Transform

QDs_inf <- function(s) {

   a <- bessel_K0(sqrt(s)) / bessel_K1(sqrt(s)) / (s ^ 1.5)
   sol <- (1 / a) / (s * s * s)
   return(sol)
}

QD_inf <- function(tD, stehfest_size) {

   v <- stehfest(stehfest_size)
   sum <- 0
   for (i in 1:stehfest_size) {
      s = i * log(2) / tD
      sum <- sum + v[i] * QDs_inf(s)
   }
   sol <- sum * log(2) / tD
   return(sol)
}

QDs <- function(s, reD) {

   if (sqrt(s) > 742) {
      ratio <- 1.0
      f <- ratio / s / sqrt(s)
      g <- 1 / f / (s * s * s)
      return(g)
   } else {
      if ((reD * sqrt(s)) > 713) {
         if (sqrt(s) > 742) {
            ratio <- 1.0
            f <- ratio / s / sqrt(s)
            g <- 1 / f / (s * s * s)
            return(g)
         } else {
            ratio <- bessel_K0(sqrt(s)) / bessel_K1(sqrt(s))
            f <- ratio / s / sqrt(s)
            g <- 1 / f / (s * s * s)
            return(g)
         }
      } else {
         a <- bessel_K0(sqrt(s)) / bessel_I0(sqrt(s))
         b <- bessel_K1(sqrt(s)) / bessel_I0(sqrt(s))
         c <- bessel_I1(sqrt(s)) / bessel_I0(sqrt(s))
         d <- bessel_K1(reD * sqrt(s)) / bessel_I1(reD * sqrt(s))
         ratio <- (a + d) / (b - d * c)
         f <- ratio / s / sqrt(s)
         g <- 1 / f / (s * s * s)
         return(g)
      }
   }
}

QD <- function(tD, reD, stehfest_size) {

   v <- stehfest(stehfest_size)
   sum <- 0
   for (i in 1:stehfest_size) {
      s = i * log(2) / tD
      sum <- sum + v[i] * QDs(s, reD)
   }
   sol <- sum * log(2) / tD
   return(sol)
}

veh_uss_rad_edge_WeD <- function(tD, reD) {

   stehfest_size <- 8
   if (reD > 10000) {
      vec <- QD_inf(tD, stehfest_size)
   } else {
      vec <- QD(tD, reD, stehfest_size)
   }
   return(vec)
}

veh_uss_rad_edge <- function(phi, perm, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure, mult_len) {

   # time in days
   len_t <- length(time)
   len_p <- length(pressure)
   if (len_t != len_p) {
      message("'time' and 'pressure' vectors must have the same length.")
      return(0)
   } else {
      tD <- vector(length = len_t)
      P_avg <- vector(length = len_p)
      dP <- vector(length = len_p)
      We <- vector(length = len_t)
      c_total <- c_water + c_rock
      f <- tetha / 360
      perm <- mult_len[2] * perm
      r_a <- mult_len[1] * r_a
      B <- 1.119 * phi * c_total * r_R * r_R * h_a * f
      tD <- 6.328e-03 * perm * time / phi / mu_water / c_total / r_R / r_R
      raD <- r_a / r_R
      for (i in 1:len_p) {
         if (i == 1) {
            P_avg[i] <- pressure[i]
            dP[i] <- 0
         } else {
            P_avg[i] <- 0.5 * (pressure[i] + pressure[i - 1])
            dP[i] <- P_avg[i - 1] - P_avg[i]
         }
      }
      for (i in seq_along(time)) {
         if (i == 1) {
            We[i] <- 0
         } else {
            dt_subset <- tD[1:i]
            dp_subset <- dP[2:i]
            len_t <- length(dt_subset)
            tD_temp <- vector(length = len_t)
            WeD <- vector(length = (len_t - 1))
            tD_temp <- (dt_subset[len_t] - dt_subset)[1:(len_t - 1)]
            for (j in 1:(len_t - 1)) {
               WeD[j] <- veh_uss_rad_edge_WeD(tD_temp[j], raD)
            }
            We[i] <- sum(B * dp_subset * WeD)
         }
      }
      return(We)
   }
}


#*******************************************************************************

# Yildiz and Khosravi Unsteady State Model - Radial Flow - Bottom Water

QDs_yk <- function(s, raD, haD) {

   sai_1 <- coth(sqrt(s) * haD) / s ^ 1.5
   k <- 100
   zeros <- bessel_zero_J1(s = c(1:k))
   beta <- zeros / raD
   sai_2 <- vector(length = k)
   sum_sai_2 <- 0
   for (i in 1:k) {
      sai_2[i] <- (bessel_J1(beta[i]) ^ 2 *
                      coth(sqrt(beta[i] * beta[i] + s) * haD)) /
         (beta[i] * beta[i] * sqrt(beta[i] * beta[i] + s) *
             bessel_J0(beta[i] * raD) ^ 2)
      sum_sai_2 <- sum_sai_2 + sai_2[i]
   }
   pDs <- sai_1 / raD / raD + 4 * sum_sai_2 / raD / raD / s
   QDs <- 1 / pDs / s / s / s / 2 / haD
   return(QDs)
}

QD_yk <- function(tD, raD, haD, stehfest_size) {

   v <- stehfest(stehfest_size)
   sum <- 0
   for (i in 1:stehfest_size) {
      s <- i * log(2) / tD
      sum <- sum + v[i] * QDs_yk(s, raD, haD)
   }
   sol <- sum * log(2) / tD
   return(sol)
}

yk_uss_rad_bottom_WeD <- function(tD, raD, haD) {

   stehfest_size <- 8
   vec <- QD_yk(tD, raD, haD, stehfest_size)
   return(vec)
}


# Yildiz and Khosravi Unsteady State Model - Radial Flow - Bottom Water

yk_uss_rad_bottom <- function(phi, perm_h, perm_v, h_a, r_a, r_R, mu_water, c_water, c_rock, time, pressure, mult_len) {

   len_t <- length(time)
   len_p <- length(pressure)
   if (len_t != len_p) {
      message("'time' and 'pressure' vectors must have the same length.")
      return(0)
   } else {
      tD <- vector(length = len_t)
      P_avg <- vector(length = len_p)
      dP <- vector(length = len_p)
      WeD <- vector(length = len_t)
      We <- vector(length = len_t)
      c_total <- c_water + c_rock
      r_a <- mult_len[1] * r_a
      perm_h <- mult_len[2] * perm_h
      perm_v <- mult_len[3] * perm_v
      Fk <- perm_v / perm_h
      B <- 1.119 * phi * c_total * r_R * r_R * h_a
      tD <- 6.328e-03 * perm_h * time / phi / mu_water / c_total / r_R / r_R
      raD <- r_a / r_R
      zD <- h_a / r_R / sqrt(Fk)
      for (i in 1:len_p) {
         if (i == 1) {
            P_avg[i] <- pressure[i]
            dP[i] <- 0
         } else {
            P_avg[i] <- 0.5 * (pressure[i] + pressure[i - 1])
            dP[i] <- P_avg[i - 1] - P_avg[i]
         }
      }
      for (i in seq_along(time)) {
         if (i == 1) {
            We[i] <- 0
         } else {
            dt_subset <- tD[1:i]
            dp_subset <- dP[2:i]
            len_t <- length(dt_subset)
            tD_temp <- vector(length = len_t)
            WeD <- vector(length = (len_t - 1))
            tD_temp <- (dt_subset[len_t] - dt_subset)[1:(len_t - 1)]
            for (j in 1:(len_t - 1)) {
               WeD[j] <- yk_uss_rad_bottom_WeD(tD_temp[j], raD, zD)
            }
            We[i] <- sum(B * dp_subset * WeD)
         }
      }
      return(We)
   }
}


#*******************************************************************************

# Fetkovich Pseudo-Steady State Model - Radial Flow - Edge Water

fetkovich_pss_rad_edge <- function(phi, perm, h_a, r_a, r_R, tetha, mu_water, c_water, c_rock, time, pressure, mult_len) {

   # time in days
   len_t <- length(time)
   len_p <- length(pressure)
   if (len_t != len_p) {
      message("'time' and 'pressure' vectors must have the same length.")
      return(0)
   } else {
      dt <- vector(length = len_t)
      p_r <- vector(length = len_t)
      p_a <- vector(length = len_t)
      dp_ar <- vector(length = len_t)
      dWe <- vector(length = len_t)
      We <- vector(length = len_t)
      c_total <- c_water + c_rock
      f <- tetha / 360
      r_a <- mult_len[1] * r_a
      perm <- mult_len[2] * perm
      raD <- r_a / r_R
      Wi <- pi * (r_a * r_a - r_R * r_R) * h_a * phi / 5.615
      Wei <- c_total * Wi * pressure[1] * f
      for (i in seq_along(time)) {
         if (i == 1) {
            dt[i] <- 0
            p_r[i] <- pressure[1]
            p_a[i] <- pressure[1]
            dp_ar[i] <- 0
            dWe[i] <- 0
            We[i] <- 0
         } else {
            J <- 0.00708 * perm * h_a * f / mu_water / (log(raD) - 0.75)
            dt[i] <- time[i] - time[i - 1]
            p_r[i] <- 0.5 * (pressure[i] + pressure[i - 1])
            dp_ar[i] <- p_a[i - 1] - p_r[i]
            dWe[i] <- Wei * dp_ar[i] * (1 - exp(-J * pressure[1] * dt[i] / Wei)) / pressure[1]
            We[i] <- We[i - 1] + dWe[i]
            p_a[i] <- pressure[1] * (1 - We[i] / Wei)
         }
      }
      return(We)
   }
}




#*******************************************************************************

# Nabor and Barham Unsteady State Model - Linear Flow - Edge Water - Exact Solution

nb_uss_lin_edge_WeD <- function(tD) {

   if (tD >= 10) {
      WeD <- 1
      return(WeD)
   } else {
      sum <- 0
      for (i in seq(1, 101, by = 2)) {
         sum <- sum + (1 / i / i) * exp(-1 * i * i * pi * pi * tD / 4)
      }
      WeD <- 1 - 8 * sum / pi / pi
      return(WeD)
   }
}


# Nabor and Barham Unsteady State Model - Linear Flow - Edge Water

nb_uss_lin_edge <- function(phi, perm, h_a, w_a, L_a, mu_water, c_water, c_rock, time, pressure, mult_len) {

   # time in days
   len_t <- length(time)
   len_p <- length(pressure)
   if (len_t != len_p) {
      message("'time' and 'pressure' vectors must have the same length.")
      return(0)
   } else {
      tD <- vector(length = len_t)
      P_avg <- vector(length = len_p)
      dP <- vector(length = len_p)
      WeD <- vector(length = len_t)
      We <- vector(length = len_t)
      c_total <- c_water + c_rock
      L_a <- mult_len[1] * L_a
      perm <- mult_len[2] * perm
      V_a <- h_a * L_a * w_a * phi
      B <- V_a * c_total / 5.615
      tD <- 6.328e-03 * perm * time / phi / mu_water / c_total / L_a / L_a
      for (i in 1:len_p) {
         if (i == 1) {
            P_avg[i] <- pressure[i]
            dP[i] <- 0
         } else {
            P_avg[i] <- 0.5 * (pressure[i] + pressure[i - 1])
            dP[i] <- P_avg[i - 1] - P_avg[i]
         }
      }
      for (i in seq_along(time)) {
         if (i == 1) {
            We[i] <- 0
         } else {
            dt_subset <- tD[1:i]
            dp_subset <- dP[2:i]
            len_t <- length(dt_subset)
            tD_temp <- vector(length = len_t)
            WeD <- vector(length = (len_t - 1))
            tD_temp <- (dt_subset[len_t] - dt_subset)[1:(len_t - 1)]
            for (j in 1:(len_t - 1)) {
               WeD[j] <- nb_uss_lin_edge_WeD(tD_temp[j])
            }
            We[i] <- sum(B * dp_subset * WeD)
         }
      }
      return(We)
   }
}




#*******************************************************************************
# Nabor and Barham Unsteady State Model - Linear Flow - Bottom Water - Exact Solution

nb_uss_lin_bottom_WeD <- function(tD) {

   if (tD >= 10) {
      WeD <- 1
      return(WeD)
   } else {
      sum <- 0
      for (i in seq(1, 101, by = 2)) {
         sum <- sum + (1 / i / i) * exp(-1 * i * i * pi * pi * tD / 4)
      }
      WeD <- 1 - 8 * sum / pi / pi
      return(WeD)
   }
}


# Nabor and Barham Unsteady State Model - Linear Flow - Bottom Water

nb_uss_lin_bottom <- function(phi, perm, h_a, w_a, L_a, mu_water, c_water, c_rock, time, pressure, mult_len) {

   # time in days
   len_t <- length(time)
   len_p <- length(pressure)
   if (len_t != len_p) {
      message("'time' and 'pressure' vectors must have the same length.")
      return(0)
   } else {
      tD <- vector(length = len_t)
      P_avg <- vector(length = len_p)
      dP <- vector(length = len_p)
      WeD <- vector(length = len_t)
      We <- vector(length = len_t)
      c_total <- c_water + c_rock
      h_a <- mult_len[1] * h_a
      perm <- mult_len[2] * perm
      V_a <- h_a * L_a * w_a * phi
      B <- V_a * c_total / 5.615
      tD <- 6.328e-03 * perm * time / phi / mu_water / c_total / h_a / h_a
      for (i in 1:len_p) {
         if (i == 1) {
            P_avg[i] <- pressure[i]
            dP[i] <- 0
         } else {
            P_avg[i] <- 0.5 * (pressure[i] + pressure[i - 1])
            dP[i] <- P_avg[i - 1] - P_avg[i]
         }
      }
      for (i in seq_along(time)) {
         if (i == 1) {
            We[i] <- 0
         } else {
            dt_subset <- tD[1:i]
            dp_subset <- dP[2:i]
            len_t <- length(dt_subset)
            tD_temp <- vector(length = len_t)
            WeD <- vector(length = (len_t - 1))
            tD_temp <- (dt_subset[len_t] - dt_subset)[1:(len_t - 1)]
            for (j in 1:(len_t - 1)) {
               WeD[j] <- nb_uss_lin_bottom_WeD(tD_temp[j])
            }
            We[i] <- sum(B * dp_subset * WeD)
         }
      }
      return(We)
   }
}




#*******************************************************************************

# Fetkovich Pseudo-Steady State Model - Linear Flow - Edge Water

fetk_pss_lin_edge <- function(phi, perm, h_a, w_a, L_a, mu_water, c_water, c_rock, time, pressure, mult_len) {

   # time in days
   len_t <- length(time)
   len_p <- length(pressure)
   if (len_t != len_p) {
      message("'time' and 'pressure' vectors must have the same length.")
      return(0)
   } else {
      dt <- vector(length = len_t)
      p_r <- vector(length = len_t)
      p_a <- vector(length = len_t)
      dp_ar <- vector(length = len_t)
      dWe <- vector(length = len_t)
      We <- vector(length = len_t)
      c_total <- c_water + c_rock
      L_a <- mult_len[1] * L_a
      perm <- mult_len[2] * perm
      Wi <- h_a * L_a * w_a * phi / 5.615
      Wei <- c_total * Wi * pressure[1]
      for (i in seq_along(time)) {
         if (i == 1) {
            dt[i] <- 0
            p_r[i] <- pressure[1]
            p_a[i] <- pressure[1]
            dp_ar[i] <- 0
            dWe[i] <- 0
            We[i] <- 0
         } else {
            J <- 0.003381 * perm * h_a * w_a / mu_water / L_a
            dt[i] <- time[i] - time[i - 1]
            p_r[i] <- 0.5 * (pressure[i] + pressure[i - 1])
            dp_ar[i] <- p_a[i - 1] - p_r[i]
            dWe[i] <- Wei * dp_ar[i] * (1 - exp(-J * pressure[1] * dt[i] / Wei)) / pressure[1]
            We[i] <- We[i - 1] + dWe[i]
            p_a[i] <- pressure[1] * (1 - We[i] / Wei)
         }
      }
      return(We)
   }
}



#*******************************************************************************

# Fetkovich Pseudo-Steady State Model - Linear Flow - Bottom Water

fetk_pss_lin_bottom <- function(phi, perm, h_a, w_a, L_a, mu_water, c_water, c_rock, time, pressure, mult_len) {

   # time in days
   len_t <- length(time)
   len_p <- length(pressure)
   if (len_t != len_p) {
      message("'time' and 'pressure' vectors must have the same length.")
      return(0)
   } else {
      dt <- vector(length = len_t)
      p_r <- vector(length = len_t)
      p_a <- vector(length = len_t)
      dp_ar <- vector(length = len_t)
      dWe <- vector(length = len_t)
      We <- vector(length = len_t)
      c_total <- c_water + c_rock
      h_a <- mult_len[1] * h_a
      perm <- mult_len[2] * perm
      Wi <- h_a * L_a * w_a * phi / 5.615
      Wei <- c_total * Wi * pressure[1]
      for (i in seq_along(time)) {
         if (i == 1) {
            dt[i] <- 0
            p_r[i] <- pressure[1]
            p_a[i] <- pressure[1]
            dp_ar[i] <- 0
            dWe[i] <- 0
            We[i] <- 0
         } else {
            J <- 0.003381 * perm * w_a * L_a / mu_water / h_a
            dt[i] <- time[i] - time[i - 1]
            p_r[i] <- 0.5 * (pressure[i] + pressure[i - 1])
            dp_ar[i] <- p_a[i - 1] - p_r[i]
            dWe[i] <- Wei * dp_ar[i] * (1 - exp(-J * pressure[1] * dt[i] / Wei)) / pressure[1]
            We[i] <- We[i - 1] + dWe[i]
            p_a[i] <- pressure[1] * (1 - We[i] / Wei)
         }
      }
      return(We)
   }
}



#*******************************************************************************

# Pot aquifer Model - Edge Water

pot <- function(phi, h_a, r_a, r_R, tetha, c_water, c_rock, pressure, mult_len) {

   len_p <- length(pressure)
   dp <- vector(length = len_p)
   We <- vector(length = len_p)
   c_total <- c_water + c_rock
   f <- tetha / 360
   r_a <- mult_len[1] * r_a
   Wi <- pi * (r_a * r_a - r_R * r_R) * h_a * phi / 5.615
   dp <- pressure[1] - pressure
   We <- c_total * Wi * f * dp
   return(We)
}
