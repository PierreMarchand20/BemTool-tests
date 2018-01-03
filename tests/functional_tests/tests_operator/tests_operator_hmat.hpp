#ifndef TESTS_OPERATOR_HPP
#define TESTS_OPERATOR_HPP

#include <bemtool/tools.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>
#include <bemtool-tests/tools.hpp>
#include <stdio.h>  /* defines FILENAME_MAX */
#include <unistd.h>
#define GetCurrentDir getcwd
#include<iostream>

std::string GetCurrentWorkingDir( void ) {
  char buff[FILENAME_MAX];
  GetCurrentDir( buff, FILENAME_MAX );
  std::string current_working_dir(buff);
  return current_working_dir;
}
using namespace bemtool;

template <int dim, typename OperatorType, typename Discretization>
Real Test_operator_hmat(Real kappa, Real radius, Real lc) {
    std::cout << GetCurrentWorkingDir() << std::endl;
    // Data
    Real kappa2 = kappa*kappa;
    int n;
    if (dim==1){
        n = 1;
    }
    else if (dim==2){
        n =0;
    }

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string meshname;
    if (dim==1){
        meshname="../../../../data/tests/functional_tests/tests_operator/circle_"+NbrToStr(lc);
    }
    else if (dim==2){
        meshname="sphere";
    }

    // Mesh
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << (meshname+".msh").c_str() << std::endl;
    Geometry node((meshname+".msh").c_str());
    Mesh<dim> mesh; mesh.Load(node,0);
    Orienting(mesh);
    int nb_elt = NbElt(mesh);
    MPI_Barrier(MPI_COMM_WORLD);

    // Dof
    Dof<Discretization> dof(mesh);
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
    BIO_Generator<OperatorType,Discretization> generator(V,dof);

    // HMatrix
    htool::HMatrix<htool::partialACA,Cplx> A(generator,x);
    A.print_infos();

    // Eigenvector
    std::vector<Cplx> En(nb_dof);
    std::vector<double> En_real(nb_dof);
    for(int j=0; j<nb_elt; j++){
        const array<dim+1,int>&  jdof = dof[j];
        const array<dim+1,R3> xdof = dof(j);
        for(int k=0; k<dim+1; k++){
            if (dim==1){
                En[jdof[k]] = exp(iu*n*std::atan2 (xdof[k][1],xdof[k][0]));
            }
            else if (dim==2){
                double rho = std::sqrt(xdof[k][0]*xdof[k][0]+xdof[k][1]*xdof[k][1]+xdof[k][2]*xdof[k][2]);
                double phi = std::atan2(xdof[k][1],xdof[k][0]);
                double theta = std::acos(xdof[k][2]/rho);
                En[jdof[k]] = boost::math::spherical_harmonic(n,n,theta,phi);
                En_real[jdof[k]] = std::real(boost::math::spherical_harmonic(n,n,theta,phi));
            }

            // double r =sqrt(xdof[k][0]*xdof[k][0]+xdof[k][1]*xdof[k][1]+xdof[k][2]*xdof[k][2]);
            // En[jdof[k]] = pow( (xdof[k][0]+iu*xdof[k][1])/r, n);
        }
    }
WritePointValGmsh(dof,"test.msh",En_real);
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


// template <typename OperatorType>
// Real Test_operator_3D_hmat(Real kappa, Real radius, Real lc) {
//   // Data
//   Real kappa2 = kappa*kappa;
//   int n = 1;
// 	// Get the rank of the process
// 	int rank;
// 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//   // Mesh
//   if (rank==0){
//     gmsh_sphere("sphere",radius,lc);
//   }
//   MPI_Barrier(MPI_COMM_WORLD);
//   Geometry node("sphere.msh");
//   Mesh2D mesh; mesh.Load(node,0);
//   Orienting(mesh);
//   int nb_elt = NbElt(mesh);
//   MPI_Barrier(MPI_COMM_WORLD);
//   if (rank==0){
//     gmsh_clean("sphere");
//   }
//
//   // Dof
//   Dof<P1_2D> dof(mesh);
//   int nb_dof = NbDof(dof);
//   std::vector<htool::R3> x(nb_dof);
//   for (int i=0;i<nb_dof;i++){
//     x[i][0]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][0];
//     x[i][1]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][1];
//     x[i][2]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][2];
//   }
//
//   // Operator
//   BIOp<OperatorType> V(mesh,mesh,kappa);
//
//   // Generator
//   BIO_Generator<OperatorType,P1_2D> generator(V,dof);
//
//   // HMatrix
//   htool::HMatrix<htool::partialACA,Cplx> A(generator,x);
//   A.print_infos();
//
//   // Eigenvector
//   std::vector<Cplx> En(nb_dof);
//   for(int j=0; j<nb_elt; j++){
//     const N2&         jdof = dof[j];
//     const array<2,R3> xdof = dof(j);
//     for(int k=0; k<2; k++){
//       En[jdof[k]] = exp(iu*n*std::atan2 (xdof[k][1],xdof[k][0]));
//       // double r =sqrt(xdof[k][0]*xdof[k][0]+xdof[k][1]*xdof[k][1]+xdof[k][2]*xdof[k][2]);
//       // En[jdof[k]] = pow( (xdof[k][0]+iu*xdof[k][1])/r, n);
//     }
//   }
//
//   // Eigenvalue
//   std::vector<Cplx> temp = A*En;
//   Cplx sum =0.;
//   for(int j=0; j<nb_dof; j++){
//       sum+= conj(En[j])*temp[j];
//   }
//
//   // Error
//   Cplx refsol = RefEigenvalue<OperatorType>::Compute(n,radius,kappa);
//   if (std::abs(refsol)<1e-16){
//     return sqrt(abs(sum - refsol));
//   }
//   else {
//     return sqrt(abs(sum - refsol)/abs(refsol));
//   }
// }

#endif
