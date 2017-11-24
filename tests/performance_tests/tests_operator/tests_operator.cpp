#include "tests_operator.hpp"

int main(int argc, char const *argv[]) {

  // Data
  Real kappa = 0.5;
  Real radius = 1.;
  Real lc = 0.0025;
  gmsh_circle("circle",radius,lc);


  Real error = 100*Test_operator_2D_dense<LA_SL_2D_P1xP1>(kappa,radius,lc);
  test = test || error<1;
  std::cout << "Relative error:\t"<<error<<" %"<<std::endl;

  
  //
  // // 2D Laplace - dense
  // Real err_LA_SL_2D_P1xP1 = 100*Test_operator_2D_dense<LA_SL_2D_P1xP1>(kappa,radius,lc);
  // Real err_LA_DL_2D_P1xP1 = 100*Test_operator_2D_dense<LA_DL_2D_P1xP1>(kappa,radius,lc);
  // Real err_LA_TDL_2D_P1xP1 = 100*Test_operator_2D_dense<LA_TDL_2D_P1xP1>(kappa,radius,lc);
  // Real err_LA_HS_2D_P1xP1 = 100*Test_operator_2D_dense<LA_HS_2D_P1xP1>(kappa,radius,lc);
  //
  // std::cout << "Relative error:\t"<<err_LA_SL_2D_P1xP1<<" %"<<std::endl;
  // std::cout << "Absolute error:\t"<<err_LA_DL_2D_P1xP1<<" %"<<std::endl;
  // std::cout << "Absolute error:\t"<<err_LA_TDL_2D_P1xP1<<" %"<<std::endl;
  // std::cout << "Relative error:\t"<<err_LA_HS_2D_P1xP1<<" %"<<std::endl;
  //
  // // 2D Helmholtz - dense
  // Real err_HE_SL_2D_P1xP1 = 100*Test_operator_2D_dense<HE_SL_2D_P1xP1>(kappa,radius,lc);
  // Real err_HE_DL_2D_P1xP1 = 100*Test_operator_2D_dense<HE_DL_2D_P1xP1>(kappa,radius,lc);
  // Real err_HE_TDL_2D_P1xP1 = 100*Test_operator_2D_dense<HE_TDL_2D_P1xP1>(kappa,radius,lc);
  // Real err_HE_HS_2D_P1xP1 = 100*Test_operator_2D_dense<HE_HS_2D_P1xP1>(kappa,radius,lc);
  //
  // std::cout << "Relative error:\t"<<err_HE_SL_2D_P1xP1<<" %"<<std::endl;
  // std::cout << "Relative error:\t"<<err_HE_DL_2D_P1xP1<<" %"<<std::endl;
  // std::cout << "Relative error:\t"<<err_HE_TDL_2D_P1xP1<<" %"<<std::endl;
  // std::cout << "Relative error:\t"<<err_HE_HS_2D_P1xP1<<" %"<<std::endl;

  return test;
}
