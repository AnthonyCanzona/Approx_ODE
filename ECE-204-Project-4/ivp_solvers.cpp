#include "ivp_solvers.hpp"
#include <cassert>
#include <cmath>
#include <algorithm>

// Euler's method
//
// f      A function point representing the right-hand side of the
//        ordinary differential equation y'(t) = f( t, y(t) )
// t_rng  A pair of entries, [t0, tf], where t0 is the initial
//        t-value and tf is the final t-value.
// y0     The initial condition y(t0) = y0
// n      The number of intervals the domain [t0, tf] is to be
//        broken into.
//
//                tf - t0
// The value h = --------- is the step size.
//                   n
//
// The next approximation is approximated using the formula
//
//       y[k + 1] = y[k] + h f( t[k], y[k] )
//
// The error is O(h^2) for one step, but only O(h) for multiple
// steps, thus, doubling 'n' (i.e., halving 'h') reduces the error
// by approximately two.

std::tuple<double *, double *, double *> euler(
  std::function<double( double t, double y )> f,
  std::pair<double, double> t_rng, double y0,
  unsigned int n
) {
  assert( n > 0 );

  double h{ (t_rng.second - t_rng.first)/n };

  double  *ts{ new double[n + 1] };
  double  *ys{ new double[n + 1] };
  double *dys{ new double[n + 1] };

  ts[0] = t_rng.first;
  ys[0] = y0;
  dys[0] = f( ts[0], ys[0] );

  for ( unsigned int k{0}; k < n; ++k ) {
     ts[k + 1] = ts[0] + h*(k + 1);
     ys[k + 1] = ys[k] + h*dys[k];
    dys[k + 1] = f( ts[k + 1], ys[k + 1] );
  }

  return std::make_tuple( ts, ys, dys );
}

// Heun's method
//
// f      A function point representing the right-hand side of the
//        ordinary differential equation y'(t) = f( t, y(t) )
// t_rng  A pair of entries, [t0, tf], where t0 is the initial
//        t-value and tf is the final t-value.
// y0     The initial condition y(t0) = y0
// n      The number of intervals the domain [t0, tf] is to be
//        broken into.
//
//                tf - t0
// The value h = --------- is the step size.
//                   n
//
// The next approximation is approximated using the formula:
//
//      s0 = f( t[k],     y[k]        )
//      s1 = f( t[k] + h, y[k] + h*s0 )
//
//                            s0 + s1
//       y[k + 1] = y[k] + h ---------
//                               2
//
// The value 's0' is the slope at (t[k], y[k]), while the value
// 's1' is an approximation of the slope at time t[k+1] = t[k] + h
// It is an approximation because we are using an approximation
// for the y value; specifically, we are approximating the y-value
// using Euler's approximation.
//
// The error is O(h^3) for one step, but only O(h^2) for multiple
// steps, thus, doubling 'n' (i.e., halving 'h') reduces the error
// by approximately four.

std::tuple<double *, double *, double *> heun(
  std::function<double( double t, double y )> f,
  std::pair<double, double> t_rng, double y0,
  unsigned int n
) {
  assert( n > 0 );

  double h{ (t_rng.second - t_rng.first)/n };

  double  *ts{ new double[n + 1] };
  double  *ys{ new double[n + 1] };
  double *dys{ new double[n + 1] };

   ts[0] = t_rng.first;
   ys[0] = y0;
  dys[0] = f( ts[0], ys[0] );

    for ( unsigned int k{0}; k < n; ++k ) {
     ts[k + 1] = ts[0] + h*(k + 1);
    double s0{ dys[k] };
    double s1{ f( ts[k + 1], ys[k] + h*dys[k] ) };
     ys[k + 1] = ys[k] + h*(s0 + s1)/2.0;
    dys[k + 1] = f( ts[k + 1], ys[k + 1] );
  }

  return std::make_tuple( ts, ys, dys );
}

// 4th-order Runge-Kutta method
//
// Note: The term 'Runge-Kutta methods' describes an entire class
//       of algorithms. This is one of the more common due to its
//       simplicitly and clarity.
//
// f      A function point representing the right-hand side of the
//        ordinary differential equation y'(t) = f( t, y(t) )
// t_rng  A pair of entries, [t0, tf], where t0 is the initial
//        t-value and tf is the final t-value.
// y0     The initial condition y(t0) = y0
// n      The number of intervals the domain [t0, tf] is to be
//        broken into.
//
//                tf - t0
// The value h = --------- is the step size.
//                   n
//
// The next approximation is approximated using the formula:
//
//      s0 = f( t[k],       y[k]          )
//      s1 = f( t[k] + h/2, y[k] + h/2*s0 )
//      s2 = f( t[k] + h/2, y[k] + h/2*s1 )
//      s3 = f( t[k] + h,   y[k] +   h*s2 )
//
//                            s0 + 2s1 + 2s2 + s3
//       y[k + 1] = y[k] + h ---------------------
//                                     6
//
// The error is O(h^5) for one step, but only O(h^4) for multiple
// steps, thus, doubling 'n' (i.e., halving 'h') reduces the error
// by approximately sixteen.

std::tuple<double *, double *, double *> rk4(
  std::function<double( double t, double y )> f,
  std::pair<double, double> t_rng, double y0,
  unsigned int n
) {
  assert( n > 0 );

  double h{ (t_rng.second - t_rng.first)/n };

  double  *ts{ new double[n + 1] };
  double  *ys{ new double[n + 1] };
  double *dys{ new double[n + 1] };

  ts[0] = t_rng.first;
  ys[0] = y0;
  dys[0] = f( ts[0], ys[0] );

  for ( unsigned int k{0}; k < n; ++k ) {
    ts[k + 1] = ts[0] + h*(k + 1);
    double s0{ dys[k] };
    double s1{ f( ts[k] + h/2.0, ys[k] + h*s0/2.0 ) };
    double s2{ f( ts[k] + h/2.0, ys[k] + h*s1/2.0 ) };
    double s3{ f( ts[k + 1],     ys[k] + h*s2 ) };

    ys[k + 1] = ys[k] + h*(
      s0 + 2.0*s1 + 2.0*s2 + s3
    )/6.0;

    dys[k + 1] = f( ts[k + 1], ys[k + 1] );
  }

  return std::make_tuple( ts, ys, dys );
}

// RMSE Using Simpsons Rule
/* 
  The function takes in the approximated solution of a
  already existing method implementation and the range    of the function. It then takes the actual solution      given to about 31 digits of precision and using the     root-mean squared method of determining error it        calculates the error, approximating the integral        using Simpsons Rule.
*/
double rmse_simpson(
  std::tuple<double *, double *, double *> aprx_solu, std::pair<double,        double> t_rng, std::function<double( double t )> a, unsigned int n
){
  assert( n > 0 || n%2 != 1 );
  double h { (t_rng.second - t_rng.first)/n };
  double *ts {std::get<0>(aprx_solu)};
  double *ys {std::get<1>(aprx_solu)};
  double integ {0};

  for ( unsigned int k{0}; k < n/2 ; k += 2 ) {
    double s0 {a(ts[k]) - ys[k]};
    double s1 {a(ts[k+1])-ys[k+1]};
    double s2 {a(ts[k+2])-ys[k+2]};
    integ += (2*h/6)*(s0*s0 + 4*s1*s1 + s2*s2);
  }
  return std::sqrt((1/n+1)*integ);
}

// RMSE Using Trapizodial Method 
/* 
  The function takes in the approximated solution of a
  already existing method implementation and the range    of the function. It then takes the actual solution      given to about 31 digits of precision and using the     root-mean squared method of determining error it        calculates the error, approximating the integral        using trapizoidal method.
*/
double rmse_trapiz(
  std::tuple<double *, double *, double *> aprx_solu, std::pair<double,        double> t_rng, std::function<double( double t )> a, unsigned int n
){
  assert( n > 0 );
  double h { (t_rng.second - t_rng.first)/n };
  double *ts {std::get<0>(aprx_solu)};
  double *ys {std::get<1>(aprx_solu)};
  double integ {0};

  for ( unsigned int k{0}; k < n; ++k ) {
    double s0 {a(ts[k])-ys[k]};
    double s1 {a(ts[k+1])-ys[k+1]};
    integ += (h/2)*(s0*s0 + s1*s1);
  }
  return std::sqrt((1/n+1)*integ);
}