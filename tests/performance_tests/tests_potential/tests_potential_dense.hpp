#ifndef TESTS_POTENTIAL_DENSE_HPP
#define TESTS_POTENTIAL_DENSE_HPP

#include <bemtool/tools.hpp>
#include <bemtool/miscellaneous/eigen_wrap.hpp>

#include <bemtool-tests/miscellaneous/gmsh_calls.hpp>
#include <bemtool-tests/miscellaneous/refsolution.hpp>

using namespace bemtool;


template <EquationEnum OperatorType, int Dimension>
double Test_potential_2D_dense(Real kappa, Real radius, Real lc, Real lc_output) {
  std::cout << radius << " "<<lc << " " << lc_output <<std::endl;
  // Data
  Real kappa2 = kappa*kappa;
  int p = 1;

  // Mesh
  gmsh_circle("circle",radius,lc);
  Geometry node("circle.msh");
  Mesh1D mesh; mesh.Load(node,0);
  Orienting(mesh);
  int nb_elt = NbElt(mesh);
  gmsh_clean("circle");

  // Mesh
  gmsh_disc("disc",radius*0.9,lc_output);
  Geometry node_output("disc.msh");
  int nb_dof_output = NbNode(node_output);
  Mesh2D mesh_output; mesh_output.Load(node_output,0);
  Orienting(mesh_output);
  gmsh_clean("disc");

  // Dof
  Dof<P1_1D> dof(mesh);
  int nb_dof = NbDof(dof);

  // Operator
  Potential<PotKernel<OperatorType,SL_POT,2,P1_1D>> SL(mesh,kappa);
  Potential<PotKernel<OperatorType,DL_POT,2,P1_1D>> DL(mesh,kappa);

  EigenDense  SL_mat(nb_dof_output,nb_dof),DL_mat(nb_dof_output,nb_dof); Clear(SL_mat);Clear(DL_mat);
  progress bar("Assemblage\t",nb_dof_output);
  for(int j=0; j<nb_dof_output; j++){
    bar++;
    R3 x = node_output[j];
    for(int k=0; k<nb_elt; k++){
        SL_mat(j,dof[k]) += SL(x,k);
        DL_mat(j,dof[k]) += DL(x,k);
    }
  }
  bar.end();

  // Trace
  std::vector<Cplx> trace_dirichlet(nb_dof,0), trace_neumann(nb_dof,0);
  for (int i=0;i<nb_dof;i++){
    R3 x;
    x[0]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][0];
    x[1]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][1];
    x[2]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][2];

    double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    double theta = std::atan2 (x[1],x[0]);
    trace_dirichlet[i] = RefSolution<OperatorType,Dimension>::Compute(r,theta,p,radius,kappa);
    trace_neumann[i] = RefSolution<OperatorType,Dimension>::ComputeDerivative(r,theta,p,radius,kappa);

  }


  // Reference
  std::vector<Cplx> sol_ref(nb_dof_output);
  for (int i=0;i<nb_dof_output;i++){
    R3 x;
    x[0]=node_output[i][0];
    x[1]=node_output[i][1];
    x[2]=node_output[i][2];
    Real r     = std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    Real theta = std::atan2 (x[1],x[0]);

    sol_ref[i]= RefSolution<OperatorType,Dimension>::Compute(r,theta,p,radius,kappa);


  }


  // Radiated field
  std::vector<Cplx> sol(nb_dof_output,0);
  std::vector<Real> sol_real(nb_dof_output,0),sol_ref_real(nb_dof_output,0);

  add_mv_prod(sol,SL_mat,trace_neumann);
  add_mv_prod(sol,DL_mat,trace_dirichlet);

  double norm =0;
  double error =0;

  for (int i =0;i<nb_dof_output;i++){
    std::cout << sol[i]<<" "<<sol_ref[i]<<std::endl;
    norm += pow(abs(sol_ref[i]),2);
    error += pow(abs(sol_ref[i]-sol[i]),2);
    sol_real[i]=std::real(sol[i]);
    sol_ref_real[i]=std::real(sol_ref[i]);
  }

  norm = sqrt(norm);
  error = sqrt(error)/norm;

  // Save
  WritePointValGmsh(mesh_output,"test.msh",sol_real);
  WritePointValGmsh(mesh_output,"ref.msh",sol_ref_real);



  // Error
  return error;
}


#endif
