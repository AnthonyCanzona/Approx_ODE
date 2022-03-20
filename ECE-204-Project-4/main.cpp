#include <iostream>
#include <cmath>
#include <tuple>
#include "ivp_solvers.hpp"

double f1( double t, double y );
double f2( double t, double y );
double y1( double t );
double y2( double t );
int main();

int main() {
  double t0{  0.0 };
  double tf{ 10.0 };  
  unsigned int N{ 32 };

  std::cout<< "euler:" << std::endl;
  std::cout<< "f1:" << std::endl;
  
  auto result_eu_1{ euler( f1, std::make_pair( 0.0, 10.0 ), 1.0, N ) };
  //double h{ (tf - t0)/N };

    auto result_eu_2{ euler( f1, std::make_pair( 0.0, 10.0 ), 1.0, 2*N ) };
  /*double h{ (tf - t0)/N };

  for ( unsigned int k{ 0 }; k <= N; ++k ) {
    std::cout << std::get<1>( result )[k] << " "
              << y1( k*h ) << std::endl;
  }*/

std::cout << "simp: " << rmse_simpson(result_eu_1, std::make_pair( 0.0, 10.0 ), y1, N) <<  " trapiz: "  << rmse_trapiz(result_eu_1, std::make_pair( 0.0, 10.0 ), y1, N) << std::endl;
  
  std::cout << "simp: " << rmse_simpson(result_eu_2, std::make_pair( 0.0, 10.0 ), y1, 2*N) <<  " trapiz: "  << rmse_trapiz(result_eu_2, std::make_pair( 0.0, 10.0 ), y1, 2*N) << std::endl;

  delete[] std::get<0>( result_eu_1 );
  delete[] std::get<1>( result_eu_1 );
  delete[] std::get<2>( result_eu_1 );
  std::get<0>( result_eu_1 ) = nullptr;
  std::get<1>( result_eu_1 ) = nullptr;
  std::get<2>( result_eu_1) = nullptr;
  
  delete[] std::get<0>( result_eu_2 );
  delete[] std::get<1>( result_eu_2 );
  delete[] std::get<2>( result_eu_2 );
  std::get<0>( result_eu_2 ) = nullptr;
  std::get<1>( result_eu_2 ) = nullptr;
  std::get<2>( result_eu_2) = nullptr;

  std::cout<< "f2:" << std::endl;
  
  result_eu_1 = euler( f2, std::make_pair( 0.0, 10.0 ), 1.0, N );
  result_eu_2 = euler( f2, std::make_pair( 0.0, 10.0 ), 1.0, 2*N );
  /*h = (tf - t0)/N;

  // This just prints out the approximations at the point t[k]
  // and the exact value when evaluating the solution.
  for ( unsigned int k{ 0 }; k <= N; ++k ) {
    std::cout << std::get<1>( result )[k] << " "
              << y2( k*h ) << std::endl;
  }*/
  
  std::cout << "simp: " << rmse_simpson(result_eu_1, std::make_pair( 0.0, 10.0 ), y2, N) <<  " trapiz: " << rmse_trapiz(result_eu_1, std::make_pair( 0.0, 10.0 ), y2, N) << std::endl;
  
  std::cout << "simp: " << rmse_simpson(result_eu_2, std::make_pair( 0.0, 10.0 ), y2, 2*N) <<  " trapiz: " << rmse_trapiz(result_eu_2, std::make_pair( 0.0, 10.0 ), y2, 2*N) << std::endl << std::endl;
  
  delete[] std::get<0>( result_eu_2 );
  delete[] std::get<1>( result_eu_2 );
  delete[] std::get<2>( result_eu_2 );
  delete[] std::get<0>( result_eu_1 );
  delete[] std::get<1>( result_eu_1 );
  delete[] std::get<2>( result_eu_1 );

    std::cout<< "heun:" << std::endl;
  std::cout<< "f1:" << std::endl;
  
  result_eu_1 = heun( f1, std::make_pair( 0.0, 10.0 ), 1.0, N );
  //double h{ (tf - t0)/N };

  result_eu_2 = heun( f1, std::make_pair( 0.0, 10.0 ), 1.0, 2*N );
  /*double h{ (tf - t0)/N };

  for ( unsigned int k{ 0 }; k <= N; ++k ) {
    std::cout << std::get<1>( result )[k] << " "
              << y1( k*h ) << std::endl;
  }*/

std::cout << "simp: " << rmse_simpson(result_eu_1, std::make_pair( 0.0, 10.0 ), y1, N) <<  " trapiz: "  << rmse_trapiz(result_eu_1, std::make_pair( 0.0, 10.0 ), y1, N) << std::endl;
  
  std::cout << "simp: " << rmse_simpson(result_eu_2, std::make_pair( 0.0, 10.0 ), y1, 2*N) <<  " trapiz: "  << rmse_trapiz(result_eu_2, std::make_pair( 0.0, 10.0 ), y1, 2*N) << std::endl;

  delete[] std::get<0>( result_eu_1 );
  delete[] std::get<1>( result_eu_1 );
  delete[] std::get<2>( result_eu_1 );
  std::get<0>( result_eu_1 ) = nullptr;
  std::get<1>( result_eu_1 ) = nullptr;
  std::get<2>( result_eu_1) = nullptr;
  
  delete[] std::get<0>( result_eu_2 );
  delete[] std::get<1>( result_eu_2 );
  delete[] std::get<2>( result_eu_2 );
  std::get<0>( result_eu_2 ) = nullptr;
  std::get<1>( result_eu_2 ) = nullptr;
  std::get<2>( result_eu_2) = nullptr;

  std::cout<< "f2:" << std::endl;
  
  result_eu_1 = heun( f2, std::make_pair( 0.0, 10.0 ), 1.0, N );
  result_eu_2 = heun( f2, std::make_pair( 0.0, 10.0 ), 1.0, 2*N );
  /*h = (tf - t0)/N;

  // This just prints out the approximations at the point t[k]
  // and the exact value when evaluating the solution.
  for ( unsigned int k{ 0 }; k <= N; ++k ) {
    std::cout << std::get<1>( result )[k] << " "
              << y2( k*h ) << std::endl;
  }*/
  
  std::cout << "simp: " << rmse_simpson(result_eu_1, std::make_pair( 0.0, 10.0 ), y2, N) <<  " trapiz: " << rmse_trapiz(result_eu_1, std::make_pair( 0.0, 10.0 ), y2, N) << std::endl;
  
  std::cout << "simp: " << rmse_simpson(result_eu_2, std::make_pair( 0.0, 10.0 ), y2, 2*N) <<  " trapiz: " << rmse_trapiz(result_eu_2, std::make_pair( 0.0, 10.0 ), y2, 2*N) << std::endl << std::endl;
  
  delete[] std::get<0>( result_eu_2 );
  delete[] std::get<1>( result_eu_2 );
  delete[] std::get<2>( result_eu_2 );
  delete[] std::get<0>( result_eu_1 );
  delete[] std::get<1>( result_eu_1 );
  delete[] std::get<2>( result_eu_1 );

    std::cout<< "rk4:" << std::endl;
  std::cout<< "f1:" << std::endl;
  
  result_eu_1 = rk4( f1, std::make_pair( 0.0, 10.0 ), 1.0, N );
  //double h{ (tf - t0)/N };

   result_eu_2 = rk4( f1, std::make_pair( 0.0, 10.0 ), 1.0, 2*N );
  /*double h{ (tf - t0)/N };

  for ( unsigned int k{ 0 }; k <= N; ++k ) {
    std::cout << std::get<1>( result )[k] << " "
              << y1( k*h ) << std::endl;
  }*/

std::cout << "simp: " << rmse_simpson(result_eu_1, std::make_pair( 0.0, 10.0 ), y1, N) <<  " trapiz: "  << rmse_trapiz(result_eu_1, std::make_pair( 0.0, 10.0 ), y1, N) << std::endl;
  
  std::cout << "simp: " << rmse_simpson(result_eu_2, std::make_pair( 0.0, 10.0 ), y1, 2*N) <<  " trapiz: "  << rmse_trapiz(result_eu_2, std::make_pair( 0.0, 10.0 ), y1, 2*N) << std::endl;

  delete[] std::get<0>( result_eu_1 );
  delete[] std::get<1>( result_eu_1 );
  delete[] std::get<2>( result_eu_1 );
  std::get<0>( result_eu_1 ) = nullptr;
  std::get<1>( result_eu_1 ) = nullptr;
  std::get<2>( result_eu_1) = nullptr;
  
  delete[] std::get<0>( result_eu_2 );
  delete[] std::get<1>( result_eu_2 );
  delete[] std::get<2>( result_eu_2 );
  std::get<0>( result_eu_2 ) = nullptr;
  std::get<1>( result_eu_2 ) = nullptr;
  std::get<2>( result_eu_2) = nullptr;

  std::cout<< "f2:" << std::endl;
  
  result_eu_1 = rk4( f2, std::make_pair( 0.0, 10.0 ), 1.0, N );
  result_eu_2 = rk4( f2, std::make_pair( 0.0, 10.0 ), 1.0, 2*N );
  /*h = (tf - t0)/N;

  // This just prints out the approximations at the point t[k]
  // and the exact value when evaluating the solution.
  for ( unsigned int k{ 0 }; k <= N; ++k ) {
    std::cout << std::get<1>( result )[k] << " "
              << y2( k*h ) << std::endl;
  }*/
  
  std::cout << "simp: " << rmse_simpson(result_eu_1, std::make_pair( 0.0, 10.0 ), y2, N) <<  " trapiz: " << rmse_trapiz(result_eu_1, std::make_pair( 0.0, 10.0 ), y2, N) << std::endl;
  
  std::cout << "simp: " << rmse_simpson(result_eu_2, std::make_pair( 0.0, 10.0 ), y2, 2*N) <<  " trapiz: " << rmse_trapiz(result_eu_2, std::make_pair( 0.0, 10.0 ), y2, 2*N) << std::endl << std::endl;
  
  delete[] std::get<0>( result_eu_2 );
  delete[] std::get<1>( result_eu_2 );
  delete[] std::get<2>( result_eu_2 );
  delete[] std::get<0>( result_eu_1 );
  delete[] std::get<1>( result_eu_1 );
  delete[] std::get<2>( result_eu_1 );

  return 0;
}

double f1( double t, double y ) {
  return -1.843*y;
}

double f2( double t, double y ) {
  return -1.887*y + 1.8784*std::cos(t);
}

double y1( double t ) noexcept {
  return 1.0*std::exp(-1.843*t);
}

double y2( double t ) {
  return 0.22281948504736810831682113257655*std::exp(-1.887*t) + 0.87956723210075657783124414004063*std::cos(t - 0.487313068541216739778297474432);
}