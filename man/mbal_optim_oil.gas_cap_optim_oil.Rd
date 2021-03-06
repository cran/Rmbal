% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mbal_oil_optim.R
\name{mbal_optim_oil.gas_cap_optim_oil}
\alias{mbal_optim_oil.gas_cap_optim_oil}
\title{S3 method for class 'mbal_optim_oil'}
\usage{
\method{mbal_optim_oil}{gas_cap_optim_oil}(optim_lst, time_lst)
}
\arguments{
\item{optim_lst}{a list object of class 'optimization_oil'}

\item{time_lst}{a list object of class 'time'}
}
\value{
a list of class 'mbal_oil' with estimates for the unknown parameters of a gas_cap_drive oil reservoir
}
\description{
Generate a list of class 'mbal_oil' with estimates for the unknown parameters of a gas_cap_drive oil reservoir
}
