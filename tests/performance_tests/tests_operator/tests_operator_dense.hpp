#ifndef TESTS_OPERATOR_HPP
#define TESTS_OPERATOR_HPP

#include <bemtool/tools.hpp>
#include <bemtool/miscellaneous/eigen_wrap.hpp>
#include <bemtool-ext/miscellaneous/gmsh_calls.hpp>

using namespace bemtool;

template <typename OperatorType>
Real Test_operator_2D_dense(Real kappa, Real radius, Real lc) {
  // Data
  Real kappa2 = kappa*kappa;
  int n = 1 ;

  // Mesh
  gmsh_circle("circle",radius,lc);
  Geometry node("circle.msh");
  Mesh1D mesh; mesh.Load(node,0);
  Orienting(mesh);
  int nb_elt = NbElt(mesh);
  gmsh_clean("circle");

  // Dof
  Dof<P1_1D> dof(mesh);
  int nb_dof = NbDof(dof);

  // Operator
  BIOp<OperatorType> V(mesh,mesh,kappa);
  EigenDense  A(nb_dof,nb_dof); Clear(A);
  progress bar("Assemblage\t",nb_elt);
  for(int j=0; j<nb_elt; j++){
    bar++;
    for(int k=0; k<nb_elt; k++){
        A(dof[j],dof[k]) += V(j,k);
    }
  }
  bar.end();

  // Eigenvector
  std::vector<Cplx> En(nb_dof);
  for(int j=0; j<nb_elt; j++){
    const N2&         jdof = dof[j];
    const array<2,R3> xdof = dof(j);
    for(int k=0; k<2; k++){
      En[jdof[k]] = exp(iu*n*std::atan2 (xdof[k][1],xdof[k][0]));
      // En[jdof[k]] = pow( xdof[k][0]+iu*xdof[k][1], n);
    }
  }

  // Eigenvalue
  Cplx sum =0.;
  for(int j=0; j<nb_dof; j++){
    for(int k=0; k<nb_dof; k++){
      sum+= A(j,k)*conj(En[j])*En[k];
    }
  }

  // Error
  Cplx refsol = RefEigenvalue<OperatorType>::Compute(n,1.,kappa);
  if (std::abs(refsol)<1e-16){
    return sqrt(abs(sum - refsol));
  }
  else {
    return sqrt(abs(sum - refsol)/abs(refsol));
  }
}

#endif
