#include <iostream>
#include <complex>
#include <vector>
#include <bemtool/tools.hpp>
#include <bemtool/miscellaneous/mat_struct.hpp>
// #include <bemtool/miscellaneous/htool_wrap.hpp>
#include <bemtool-tests/tools.hpp>


using namespace bemtool;

template <int dim, typename Discretization>
Real Test_operator_slobodeckij(Real radius, Real lc){

    // // Get the rank of the process
    // int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Data
    std::string meshname;
    if (dim==1){
        // if (rank==0){
            gmsh_circle("circle",radius,lc);
        // }
        meshname = "circle";
        // meshname="../../../../data/tests/functional_tests/tests_operator/circle_"+NbrToStr(lc);
    }
    else if (dim==2){
        // if (rank==0){
            gmsh_sphere("sphere",radius,lc);
        // }
        meshname = "sphere";
        // meshname="../../../../data/tests/functional_tests/tests_operator/sphere_"+NbrToStr(lc);
    }

	// Mesh
    // MPI_Barrier(MPI_COMM_WORLD);
    Geometry node((meshname+".msh").c_str());
    Mesh<dim> mesh; mesh.Load(node);
    Orienting(mesh);
    int nb_elt = NbElt(mesh);
    // MPI_Barrier(MPI_COMM_WORLD);

    // Dof
    Dof<Discretization> dof(mesh);
    int nb_dof = NbDof(dof);
    // std::vector<htool::R3> x(nb_dof);
    // for (int i=0;i<nb_dof;i++){
    //     x[i][0]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][0];
    //     x[i][1]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][1];
    //     x[i][2]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][2];
    // }


    // BIOp_SLO
    BIOp_SLO<dim,Discretization,Discretization> biop(mesh,mesh);

    // // Generator
    // BIO_Generator_slo<dim,Discretization> generator(dof);

    // // HMatrix
    // std::cout << "HMATRIX"<<std::endl;
    // htool::HMatrix<htool::partialACA,Cplx> A(generator,x);
    // A.print_infos();

    // Building matrix
    // SubBIOp_slo<dim,Discretization,Discretization> subbiop(dof,dof);
    DenseMatrix<Cplx>  A(nb_dof,nb_dof);
    // std::vector<int> dofs(nb_dof);
    // std::iota (std::begin(dofs), std::end(dofs), 0);
    // subbiop.compute_neumann_block(dofs,dofs,A);
std::cout <<"ok"<<std::endl;
    progress bar("Assemblage\t",nb_elt);
    for(int j=0; j<nb_elt; j++){
        bar++;
        for(int k=0; k<nb_elt; k++){

            BlockMat matloc = biop(j,k);
            const array<dim+1,int>&  jdof = dof[j];
            const array<dim+1,int>&  kdof = dof[k];
            std::vector<int> slo_num = biop.get_slo_num();
            std::vector<int> I(NbRow(matloc),0);
            for (int l=0;l<size(jdof);l++){
                I[slo_num[l]]=jdof[l];
                I[slo_num[l+size(jdof)]]=kdof[l];
            } 

            for (int l=0; l<I.size();l++){
                for (int m=0; m<I.size();m++){
                   A(I[l],I[m])+=matloc(l,m); 
                }
            }
        }
    }
    bar.end();


    // Eigenvector
    std::vector<Cplx> En_1(nb_dof),En_2(nb_dof);
    std::vector<double> En_real(nb_dof);

    std::pair<int,int> n_1(3,1);
    std::pair<int,int> n_2(2,2);
    int m_1 = 1;
    int m_2 = 3;

    for(int j=0; j<nb_elt; j++){
        const array<dim+1,int>&  jdof = dof[j];
        const array<dim+1,R3> xdof = dof(j);
        for(int k=0; k<dim+1; k++){
            if (dim==1){
                En_1[jdof[k]] = exp(iu*m_1*std::atan2 (xdof[k][1],xdof[k][0]));
            }
            else if (dim==2){
                double rho = std::sqrt(xdof[k][0]*xdof[k][0]+xdof[k][1]*xdof[k][1]+xdof[k][2]*xdof[k][2]);
                double phi = std::atan2(xdof[k][1],xdof[k][0]);
                double theta = std::acos(xdof[k][2]/rho);
                En_1[jdof[k]] = boost::math::spherical_harmonic(n_1.first,n_1.second,theta,phi);
                // En_real[jdof[k]] = std::real(boost::math::spherical_harmonic(n_1.first,n_1.second,theta,phi));
            }

        }
    }

    for(int j=0; j<nb_elt; j++){
        const array<dim+1,int>&  jdof = dof[j];
        const array<dim+1,R3> xdof = dof(j);
        for(int k=0; k<dim+1; k++){
            if (dim==1){
                En_2[jdof[k]] = exp(iu*m_2*std::atan2 (xdof[k][1],xdof[k][0]));
            }
            else if (dim==2){
                double rho = std::sqrt(xdof[k][0]*xdof[k][0]+xdof[k][1]*xdof[k][1]+xdof[k][2]*xdof[k][2]);
                double phi = std::atan2(xdof[k][1],xdof[k][0]);
                double theta = std::acos(xdof[k][2]/rho);
                En_2[jdof[k]] = boost::math::spherical_harmonic(n_2.first,n_2.second,theta,phi);
                // En_real[jdof[k]] = std::real(boost::math::spherical_harmonic(n_2.first,n_2.second,theta,phi));
            }

        }
    }

    // Compute orthogonality
    std::vector<Cplx> temp(nb_dof);
    mv_prod(temp,A,En_1);
    // std::vector<Cplx> temp = A*En_1;
    Cplx result=0;
    for (int i=0;i<nb_dof;i++){
        result+= temp[i]*std::conj(En_2[i]);
    }
    return std::abs(result);


}
