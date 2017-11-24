#ifndef TESTS_POTENTIAL_HMAT_HPP
#define TESTS_POTENTIAL_HMAT_HPP

#include <bemtool/tools.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>

#include <bemtool-ext/miscellaneous/gmsh_calls.hpp>
#include <bemtool-ext/miscellaneous/refsolution.hpp>
using namespace bemtool;


template <EquationEnum EquationType, int Dimension>
double Test_potential_2D_hmat(Real kappa, Real radius, Real lc, Real lc_output) {
  // Data
  Real kappa2 = kappa*kappa;
  int p = 1;

  // Get the rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Mesh
  if (rank==0){
    gmsh_circle("circle",radius,lc);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  Geometry node("circle.msh");
  Mesh1D mesh; mesh.Load(node,0);
  Orienting(mesh);
  int nb_elt = NbElt(mesh);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank==0){
    gmsh_clean("circle");
  }

  // Mesh
  if (rank==0){
    gmsh_disc("disc",radius*0.9,lc_output);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  Geometry node_output("disc.msh");
  int nb_dof_output = NbNode(node_output);
  Mesh2D mesh_output; mesh_output.Load(node_output,0);
  Orienting(mesh_output);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank==0){
    gmsh_clean("disc");
  }
  std::vector<htool::R3> x_output(nb_dof_output);
  for (int i=0;i<nb_dof_output;i++){
    x_output[i][0]=node_output[i][0];;
    x_output[i][1]=node_output[i][1];;
    x_output[i][2]=node_output[i][2];;
  }


  // Dof
  Dof<P1_1D> dof(mesh);
  int nb_dof = NbDof(dof);
  std::vector<htool::R3> x(nb_dof);
  for (int i=0;i<nb_dof;i++){
    x[i][0]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][0];
    x[i][1]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][1];
    x[i][2]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][2];
  }

  // Operator
  Potential<PotKernel<EquationType,SL_POT,2,P1_1D>> SL(mesh,kappa);
  Potential<PotKernel<EquationType,DL_POT,2,P1_1D>> DL(mesh,kappa);

  // POT_Generator
  POT_Generator<PotKernel<EquationType,SL_POT,2,P1_1D>,P1_1D> generator_SL(SL,dof,node_output);
  POT_Generator<PotKernel<EquationType,DL_POT,2,P1_1D>,P1_1D> generator_DL(DL,dof,node_output);

  // Cluster trees
  std::shared_ptr<htool::Cluster_tree> t=std::make_shared<htool::Cluster_tree>(x_output);
	std::shared_ptr<htool::Cluster_tree> s=std::make_shared<htool::Cluster_tree>(x);

  htool::HMatrix<htool::partialACA,Cplx> SL_hmat(generator_SL,t,x_output,s,x);
  htool::HMatrix<htool::partialACA,Cplx> DL_hmat(generator_DL,t,x_output,s,x);
  SL_hmat.print_infos();
  DL_hmat.print_infos();

  // Trace
  std::vector<Cplx> trace_dirichlet(nb_dof,0), trace_neumann(nb_dof,0);
  for (int i=0;i<nb_dof;i++){
    R3 x;
    x[0]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][0];
    x[1]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][1];
    x[2]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][2];

    double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    double theta = std::atan2 (x[1],x[0]);
    trace_dirichlet[i] = RefSolution<EquationType,Dimension>::Compute(r,theta,p,radius,kappa);
    std::cout << trace_dirichlet[i] <<std::endl;
    trace_neumann[i] = RefSolution<EquationType,Dimension>::ComputeDerivative(r,theta,p,radius,kappa);

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

    sol_ref[i]= RefSolution<EquationType,Dimension>::Compute(r,theta,p,radius,kappa);


  }


  // Radiated field
  std::vector<Cplx> sol_SL=SL_hmat*trace_neumann;
  std::vector<Cplx> sol_DL=DL_hmat*trace_dirichlet;
  std::vector<Real> sol_real(nb_dof_output,0),sol_ref_real(nb_dof_output,0);

  // add_mv_prod(sol,SL_mat,trace_neumann);
  // add_mv_prod(sol,DL_mat,trace_dirichlet);

  double norm =0;
  double error =0;

  for (int i =0;i<nb_dof_output;i++){
    // std::cout << sol_SL[i]+sol_DL[i]<<" "<<sol_ref[i]<<std::endl;
    norm += pow(abs(sol_ref[i]),2);
    error += pow(abs(sol_ref[i]-sol_SL[i]-sol_DL[i]),2);
    sol_real[i]=std::real(sol_SL[i]+sol_DL[i])/(pi);
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
