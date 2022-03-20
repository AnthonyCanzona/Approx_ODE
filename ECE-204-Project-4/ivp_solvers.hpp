#pragma once
#include <tuple>
#include <utility>
#include <functional>

// Detailed algorithm descriptions are in the source file.

  /////////////////////////////////////////////////////////
 // The various IVP solving algorithms covered in class //
/////////////////////////////////////////////////////////

// Euler's method
std::tuple<double *, double *, double *> euler(
  std::function<double( double t, double y )> f,
  std::pair<double, double> t_rng, double y0,
  unsigned int n
);

// Huen's method
std::tuple<double *, double *, double *> heun(
  std::function<double( double t, double y )> f,
  std::pair<double, double> t_rng, double y0,
  unsigned int n
);

// 4th-order Runge-Kutta method
std::tuple<double *, double *, double *> rk4(
  std::function<double( double t, double y )> f,
  std::pair<double, double> t_rng, double y0,
  unsigned int n
);

// RMSE Using Trapizodial Method 
double rmse_simpson(
  std::tuple<double *, double *, double *> aprx_solu, std::pair<double,        double> t_rng, std::function<double( double t )> a, unsigned int n
);

// RMSE Using Trapizodial Method 
double rmse_trapiz(
  std::tuple<double *, double *, double *> aprx_solu, std::pair<double,        double> t_rng, std::function<double( double t )> a, unsigned int n  
);