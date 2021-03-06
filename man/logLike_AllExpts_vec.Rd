% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCR.R
\name{logLike_AllExpts_vec}
\alias{logLike_AllExpts_vec}
\title{MCR model log likelihood for all experiments with vector of parameters (interface to C)}
\usage{
logLike_AllExpts_vec(
  theta,
  GFPp_ym_1_F,
  GFPp_yp_1_F,
  GFPm_ym_1_F,
  GFPm_yp_1_F,
  GFPp_ym_1_M,
  GFPp_yp_1_M,
  GFPm_ym_1_M,
  GFPm_yp_1_M,
  GFPp_ym_2_F,
  GFPp_yp_2_F,
  GFPm_ym_2_F,
  GFPm_yp_2_F,
  GFPp_ym_2_M,
  GFPp_yp_2_M,
  GFPm_ym_2_M,
  GFPm_yp_2_M,
  GFPp_ym_3_F,
  GFPp_yp_3_F,
  GFPm_ym_3_F,
  GFPm_yp_3_F,
  GFPp_ym_3_M,
  GFPp_yp_3_M,
  GFPm_ym_3_M,
  GFPm_yp_3_M,
  GFPp_ym_4_F,
  GFPp_yp_4_F,
  GFPm_ym_4_F,
  GFPm_yp_4_F,
  GFPp_ym_4_M,
  GFPp_yp_4_M,
  GFPm_ym_4_M,
  GFPm_yp_4_M,
  gens,
  neg = FALSE
)
}
\arguments{
\item{neg}{if \code{TRUE}, return negative log-likelihood (for minimization algorithms); otherwise return log-likelihood (for maximization)}
}
\description{
MCR model log likelihood for all experiments with vector of parameters (interface to C)
}
