#include "tests_potential_dense.hpp"

using namespace bemtool;
int main(int argc, char *argv[]) {
  // Data
  Real kappa = 0.5;
  Real radius = 1.;
  Real lc = 0.0025;
	Real lc_output = 0.05;
  bool test =0;

  // Test
  tic();
  double error = 100*Test_potential_2D_dense<HE,2>(kappa, radius, lc, lc_output);
  toc();
  test = test || error>1;
	std::cout << "Relative error:\t"<<error<<" %"<<std::endl;

  return test;
}
