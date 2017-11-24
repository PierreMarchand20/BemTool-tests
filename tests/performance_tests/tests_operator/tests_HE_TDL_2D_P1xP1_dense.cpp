#include "tests_operator_dense.hpp"

int main(int argc, char const *argv[]) {

  // Data
  Real kappa = 0.5;
  Real radius = 1.;
  Real lc = 0.0025;

  // Test
  int test = 0;
  tic();
  Real error = 100*Test_operator_2D_dense<HE_TDL_2D_P1xP1>(kappa,radius,lc);
  toc();
  test = test || error>1;
  std::cout << "Relative error:\t"<<error<<" %"<<std::endl;

  return test;
}
