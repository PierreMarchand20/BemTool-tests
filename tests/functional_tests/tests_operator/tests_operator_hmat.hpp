#ifndef TESTS_OPERATOR_HPP
#define TESTS_OPERATOR_HPP

#include <bemtool/tools.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>
#include <bemtool-tests/tools.hpp>

using namespace bemtool;

template <int dim, typename OperatorType, typename Discretization>
Real Test_operator_hmat(Real kappa, Real radius, Real lc) {
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
        if (rank==0){
            gmsh_circle("circle",radius,lc);
        }
        meshname = "circle";
        // meshname="../../../../data/tests/functional_tests/tests_operator/circle_"+NbrToStr(lc);
    }
    else if (dim==2){
        if (rank==0){
            gmsh_sphere("sphere",radius,lc);
        }
        meshname = "sphere";
        // meshname="../../../../data/tests/functional_tests/tests_operator/sphere_"+NbrToStr(lc);
    }

    // Mesh
    MPI_Barrier(MPI_COMM_WORLD);
    Geometry node((meshname+".msh").c_str());
    Mesh<dim> mesh; mesh.Load(node);
    Orienting(mesh);
    int nb_elt = NbElt(mesh);
    MPI_Barrier(MPI_COMM_WORLD);

    // Dof
    Dof<Discretization> dof(mesh);
    int nb_dof = NbDof(dof);
    std::vector<double> x(3*nb_dof);
    for (int i=0;i<nb_dof;i++){
        x[3*i+0]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][0];
        x[3*i+1]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][1];
        x[3*i+2]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][2];
    }

    // Generator
    BIO_Generator<OperatorType,Discretization> generator(dof,kappa);

    // Clustering
    if (rank == 0)
        std::cout << "Creating cluster tree" << std::endl;
    std::shared_ptr<htool::Cluster<htool::PCAGeometricClustering>> t = std::make_shared<htool::Cluster<htool::PCAGeometricClustering>>();
    std::vector<int> tab(nb_dof);
    t->build(nb_dof, x.data(),2);

    // HMatrix
    htool::HMatrix<Cplx> A(t,t);
    A.set_eta(-1);
    A.build(generator,x.data());
    A.print_infos();

    // Eigenvector
    std::vector<Cplx> En(nb_dof);
    std::vector<double> En_real(nb_dof);
    for(int j=0; j<nb_elt; j++){
        const array<Dof<Discretization>::Trait::nb_dof_loc,int>&  jdof = dof[j];
        const array<Dof<Discretization>::Trait::nb_dof_loc,R3> xdof = dof(j);
        for(int k=0; k<Dof<Discretization>::Trait::nb_dof_loc; k++){
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

    // Eigenvalue
    std::vector<Cplx> temp = A*En;
    Cplx sum =0.;
    for(int j=0; j<nb_dof; j++){
        sum+= conj(En[j])*temp[j];
    }

    // Error
    Cplx refsol = RefEigenvalue<OperatorType>::Compute(n,radius,kappa);
    std::cout << sum <<  " "<< refsol << std::endl;
    if (std::abs(refsol)<1e-16){
        return sqrt(abs(sum - refsol));
    }
    else {
        return sqrt(abs(sum - refsol)/abs(refsol));
    }
}

#endif
