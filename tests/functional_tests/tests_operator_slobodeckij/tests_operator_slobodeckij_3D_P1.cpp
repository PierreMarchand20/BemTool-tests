#include "tests_operator_slobodeckij.hpp"

int main(int argc, char *argv[]) {
  // Data
  Real radius = 0.5;
  Real lc = 0.09;

  return !(Test_operator_slobodeckij<2,P1_2D>(radius,lc)<0.01);
}
