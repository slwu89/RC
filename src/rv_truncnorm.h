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

#ifndef TRUNCNORM_H
#define TRUNCNORM_H

/* C headers */
#include <math.h>
#include <float.h>

 /* R's C API */
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>

/* constants */
static const double t1 = 0.15;
static const double t2 = 2.18;
static const double t3 = 0.725;
static const double t4 = 0.45;

/* Exponential rejection sampling (a,inf) */
double ers_a_inf(double a);

/* Exponential rejection sampling (a,b) */
double ers_a_b(double a, double b);

/* Normal rejection sampling (a,b) */
double nrs_a_b(double a, double b);

/* Normal rejection sampling (a,inf) */
double nrs_a_inf(double a);

/* Half-normal rejection sampling */
double hnrs_a_b(double a, double b);

/* Uniform rejection sampling */
double urs_a_b(double a, double b);

/* Previously this was refered to as type 1 sampling: */
double r_lefttruncnorm(double a, double mean, double sd);
double r_righttruncnorm(double b, double mean, double sd);

/* sampler for univariate truncated Gaussian */
double r_truncnorm(double a, double b, double mean, double sd);

#endif
