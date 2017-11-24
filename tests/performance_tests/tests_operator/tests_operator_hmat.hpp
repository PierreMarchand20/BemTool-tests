#ifndef TESTS_OPERATOR_HPP
#define TESTS_OPERATOR_HPP

#include <bemtool/tools.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>
#include <bemtool-ext/miscellaneous/gmsh_calls.hpp>

using namespace bemtool;

template <typename OperatorType>
Real Test_operator_2D_hmat(Real kappa, Real radius, Real lc) {
  // Data
  Real kappa2 = kappa*kappa;
  int n = 1;
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
  BIOp<OperatorType> V(mesh,mesh,kappa);

  // Generator
  BIO_Generator<OperatorType,P1_1D> generator(V,dof);

  // HMatrix
  htool::HMatrix<htool::partialACA,Cplx> A(generator,x);
  A.print_infos();

  // Eigenvector
  std::vector<Cplx> En(nb_dof);
  for(int j=0; j<nb_elt; j++){
    const N2&         jdof = dof[j];
    const array<2,R3> xdof = dof(j);
    for(int k=0; k<2; k++){
      En[jdof[k]] = exp(iu*n*std::atan2 (xdof[k][1],xdof[k][0]));
      // double r =sqrt(xdof[k][0]*xdof[k][0]+xdof[k][1]*xdof[k][1]+xdof[k][2]*xdof[k][2]);
      // En[jdof[k]] = pow( (xdof[k][0]+iu*xdof[k][1])/r, n);
    }
  }

  // Eigenvalue
  std::vector<Cplx> temp = A*En;
  Cplx sum =0.;
  for(int j=0; j<nb_dof; j++){
      sum+= conj(En[j])*temp[j];
  }

  // Error
  Cplx refsol = RefEigenvalue<OperatorType>::Compute(n,radius,kappa);
  if (std::abs(refsol)<1e-16){
    return sqrt(abs(sum - refsol));
  }
  else {
    return sqrt(abs(sum - refsol)/abs(refsol));
  }
}


template <typename OperatorType>
Real Test_operator_3D_hmat(Real kappa, Real radius, Real lc) {
  // Data
  Real kappa2 = kappa*kappa;
  int n = 1;
	// Get the rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Mesh
  if (rank==0){
    gmsh_sphere("sphere",radius,lc);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  Geometry node("sphere.msh");
  Mesh2D mesh; mesh.Load(node,0);
  Orienting(mesh);
  int nb_elt = NbElt(mesh);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank==0){
    gmsh_clean("sphere");
  }

  // Dof
  Dof<P1_2D> dof(mesh);
  int nb_dof = NbDof(dof);
  std::vector<htool::R3> x(nb_dof);
  for (int i=0;i<nb_dof;i++){
    x[i][0]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][0];
    x[i][1]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][1];
    x[i][2]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][2];
  }

  // Operator
  BIOp<OperatorType> V(mesh,mesh,kappa);

  // Generator
  BIO_Generator<OperatorType,P1_2D> generator(V,dof);

  // HMatrix
  htool::HMatrix<htool::partialACA,Cplx> A(generator,x);
  A.print_infos();

  // Eigenvector
  std::vector<Cplx> En(nb_dof);
  for(int j=0; j<nb_elt; j++){
    const N2&         jdof = dof[j];
    const array<2,R3> xdof = dof(j);
    for(int k=0; k<2; k++){
      En[jdof[k]] = exp(iu*n*std::atan2 (xdof[k][1],xdof[k][0]));
      // double r =sqrt(xdof[k][0]*xdof[k][0]+xdof[k][1]*xdof[k][1]+xdof[k][2]*xdof[k][2]);
      // En[jdof[k]] = pow( (xdof[k][0]+iu*xdof[k][1])/r, n);
    }
  }

  // Eigenvalue
  std::vector<Cplx> temp = A*En;
  Cplx sum =0.;
  for(int j=0; j<nb_dof; j++){
      sum+= conj(En[j])*temp[j];
  }

  // Error
  Cplx refsol = RefEigenvalue<OperatorType>::Compute(n,radius,kappa);
  if (std::abs(refsol)<1e-16){
    return sqrt(abs(sum - refsol));
  }
  else {
    return sqrt(abs(sum - refsol)/abs(refsol));
  }
}

#endif
