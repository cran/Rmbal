% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mbal_gas_forecast.R
\name{mbal_forecast_gas.volumetric_forecast_gas}
\alias{mbal_forecast_gas.volumetric_forecast_gas}
\title{S3 method for class 'mbal_forecast_gas'}
\usage{
\method{mbal_forecast_gas}{volumetric_forecast_gas}(forecast_lst, time_lst)
}
\arguments{
\item{forecast_lst}{a list object of class 'forecast_gas'}

\item{time_lst}{a list object of class 'time'}
}
\value{
a data frame with estimates for saturation of fluids, liquid dropout, gas-oil ratio, recovery factor, and drive indices over a range of given pressures for a volumetric reservoir
}
\description{
Return a data frame with estimates for saturation of fluids, liquid dropout, gas-oil ratio, recovery factor, and drive indices over a range of given pressures for a volumetric reservoir
}
