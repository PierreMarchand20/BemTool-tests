#include "tests_operator_slobodeckij.hpp"

int main(int argc, char *argv[]) {
  // Data
  Real radius = 0.5;
  Real lc = 0.005;

  return !(Test_operator_slobodeckij<1,P1_1D>(radius,lc)<0.01);
}
