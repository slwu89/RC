/* ##################################################
#   Truncated Gaussian Distribution
#   Code Citation:
#
#     Olaf Mersmann, Heike Trautmann,
#     Detlef Steuer and Björn Bornkamp
#     (2018). truncnorm: Truncated Normal Distribution.
#     R package version 1.0-8.
#     https://CRAN.R-project.org/package=truncnorm
#
#   Sean Wu
################################################## */

#include "rv_truncnorm.h"

/* Exponential rejection sampling (a,inf) */
double ers_a_inf(double a) {

  const double ainv = 1.0 / a;
  double x;
  double rho;

  do {
    x = rexp(ainv) + a; /* rexp works with 1/lambda */
    rho = exp(-0.5 * pow((x - a), 2));
  } while (runif(0, 1) > rho);
  return x;
}

/* Exponential rejection sampling (a,b) */
double ers_a_b(double a, double b) {
  const double ainv = 1.0 / a;
  double x, rho;
  do {
    x = rexp(ainv) + a; /* rexp works with 1/lambda */
    rho = exp(-0.5 * pow((x - a), 2));
  } while (runif(0, 1) > rho || x > b);
  return x;
}

/* Normal rejection sampling (a,b) */
double nrs_a_b(double a, double b) {
  double x = -DBL_MAX;
  while (x < a || x > b) {
    x = rnorm(0, 1);
  }
  return x;
}

/* Normal rejection sampling (a,inf) */
double nrs_a_inf(double a) {
  double x = -DBL_MAX;
  while (x < a) {
    x = rnorm(0, 1);
  }
  return x;
}

/* Half-normal rejection sampling */
double hnrs_a_b(double a, double b) {
  double x = a - 1.0;
  while (x < a || x > b) {
    x = rnorm(0, 1);
    x = fabs(x);
  }
  return x;
}

/* Uniform rejection sampling */
double urs_a_b(double a, double b) {
  const double phi_a = dnorm(a, 0.0, 1.0, FALSE);
  double x = 0.0;

  /* Upper bound of normal density on [a, b] */
  const double ub = a < 0 && b > 0 ? M_1_SQRT_2PI : phi_a;
  do {
    x = runif(a, b);
  } while (runif(0, 1) * ub > dnorm(x, 0, 1, 0));
  return x;
}

/* Previously this was refered to as type 1 sampling: */
double r_lefttruncnorm(double a, double mean, double sd) {
  const double alpha = (a - mean) / sd;
  double x;
  if (alpha < t4) {
    x = mean + sd * nrs_a_inf(alpha);
    return x;
  } else {
    x = mean + sd * ers_a_inf(alpha);
    return x;
  }
}

double r_righttruncnorm(double b, double mean, double sd) {
  const double beta = (b - mean) / sd;
  /* Exploit symmetry: */
  return mean - sd * r_lefttruncnorm(-beta, 0.0, 1.0);
}

/* sampler for univariate truncated Gaussian */
double r_truncnorm(double a, double b, double mean, double sd){

  const double alpha = (a - mean) / sd;
  const double beta = (b - mean) / sd;
  const double phi_a = dnorm(alpha, 0.0, 1.0, FALSE);
  const double phi_b = dnorm(beta, 0.0, 1.0, FALSE);

  double out;
  if(beta <= alpha){
    out = NAN;
    return out;
  } else if (alpha <= 0 && 0 <= beta) {
    if (phi_a <= t1 || phi_b <= t1) {
      return mean + sd * nrs_a_b(alpha, beta);
    } else {
      return mean + sd * urs_a_b(alpha, beta);
    }
  } else if (alpha > 0) {
    if (phi_a / phi_b <= t2) {
      return mean + sd * urs_a_b(alpha, beta);
    } else {
      if (alpha < t3) {
        return mean + sd * hnrs_a_b(alpha, beta);
      } else {
        return mean + sd * ers_a_b(alpha, beta);
      }
    }
  } else {
    if (phi_b / phi_a <= t2) {
      return mean - sd * urs_a_b(-beta, -alpha);
    } else {
      if (beta > -t3) {
        return mean - sd * hnrs_a_b(-beta, -alpha);
      } else {
        return mean - sd * ers_a_b(-beta, -alpha);
      }
    }
  }
};
