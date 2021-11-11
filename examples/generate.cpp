#include <bemtool/tools.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>
#include <bemtool-tests/tools.hpp>

using namespace bemtool;

int main(int argc, char *argv[]) {
    // Initialize the MPI environment
    MPI_Init(&argc,&argv);

    // Get the rank of the process
    int rank,sizeWorld;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);

    // Command line
    // Check the number of parameters
    if (argc < 2) {
        // Tell the user how to run the program
        std::cerr << "Usage: " << argv[0] << " meshname" << std::endl;
        /* "Usage messages" are a conventional way of telling the user
         * how to run a program if they enter the command incorrectly.
         */
        return 1;
    }
    std::string meshname = argv[1];
    double kappa = 10;
    int overlap = 4;
    R3 dir; dir[0]=1./std::sqrt(2);dir[1]=1./std::sqrt(2);dir[2]=0;

    // HTOOL
    htool::SetNdofPerElt(1);
    htool::SetEpsilon(0.1);
    htool::SetEta(-1);

    // HPDDM
    HPDDM::Option& opt = *HPDDM::Option::get();
    opt.parse(argc, argv);
    if(rank != 0)
        opt.remove("verbosity");

    // Mesh
    if (rank==0)
        std::cout << "Loading mesh" << std::endl;
    Geometry node(meshname);
    Mesh1D mesh; mesh.Load(node,0);
    Orienting(mesh);
    mesh = unbounded;
    int nb_elt = NbElt(mesh);
    std::vector<R3> normal = NormalTo(mesh);

    // Dof
    if (rank==0)
        std::cout << "Create dof" << std::endl;
    Dof<P1_1D> dof(mesh);
    int nb_dof = NbDof(dof);
    std::vector<htool::R3> x(nb_dof);
    for (int i=0;i<nb_dof;i++){
      x[i][0]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][0];
      x[i][1]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][1];
      x[i][2]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][2];
    }
    WritePointValGmsh(dof,"mesh.msh",std::vector<int>(nb_dof,1));

    // Generator
    if (rank==0)
       std::cout << "Creating generators"<<std::endl;
    BIO_Generator<HE_HS_2D_P1xP1,P1_1D> generator_W(dof,kappa);

    // Clustering
    if (rank==0)
	   std::cout << "Creating cluster tree" << std::endl;
    std::shared_ptr<htool::GeometricClustering> t=std::make_shared<htool::GeometricClustering>();
    std::vector<int>tab(nb_dof);
    std::iota(tab.begin(),tab.end(),int(0));
    t->build(x,std::vector<double>(nb_dof,0),tab,std::vector<double>(nb_dof,1));


    // HMatrix
    if (rank==0)
	   std::cout << "Building Hmatrix" << std::endl;
    htool::HMatrix<Cplx,htool::partialACA,htool::GeometricClustering> W(generator_W,t,x);
    W.print_infos();

    // Right-hand side
    if (rank==0)
        std::cout << "Building rhs" << std::endl;
    std::vector<Cplx> rhs(nb_dof,0);
    std::vector<double> rhs_real(nb_dof,0),rhs_abs(nb_dof,0);
    for (int i=0;i<nb_elt;i++){
        const N2&         jdof = dof[i];
        const array<2,R3> xdof = dof(i);
        R2x2 M_local = MassP1(mesh[i]);
        C2 Uinc;
        Uinc[0]= iu*kappa*(dir,normal[i])*exp( iu*kappa*(xdof[0],dir) );
        Uinc[1]= iu*kappa*(dir,normal[i])*exp( iu*kappa*(xdof[1],dir) );

        for(int k=0;k<2;k++){
            rhs[jdof[k]] -= (M_local(k,0)*Uinc[0]+M_local(k,1)*Uinc[1]);
        }
    }

    if (rank==0){
        htool::vector_to_bytes(rhs,"rhs.bin");
    }

    // Overlap
    if (rank==0)
        std::cout << "Building partitions" << std::endl;
    std::vector<int> cluster_to_ovr_subdomain;
    std::vector<int> ovr_subdomain_to_global;
    std::vector<int> neighbors;
    std::vector<std::vector<int> > intersections;

    Partition(W.get_MasterOffset_t(), W.get_permt(),dof,cluster_to_ovr_subdomain,ovr_subdomain_to_global,neighbors,intersections);


    W.get_cluster_tree_t().save_cluster("cluster_"+NbrToStr(sizeWorld));
    htool::GeometricClustering test;

    htool::vector_to_bytes(cluster_to_ovr_subdomain,"cluster_to_ovr_subdomain_"+NbrToStr(sizeWorld)+"_"+NbrToStr(rank)+".bin");
    htool::vector_to_bytes(ovr_subdomain_to_global,"ovr_subdomain_to_global_"+NbrToStr(sizeWorld)+"_"+NbrToStr(rank)+".bin");
    htool::vector_to_bytes(neighbors,"neighbors_"+NbrToStr(sizeWorld)+"_"+NbrToStr(rank)+".bin");

    for (int i = 0 ;i<intersections.size();i++){
        htool::vector_to_bytes(intersections[i],"intersections_"+NbrToStr(sizeWorld)+"_"+NbrToStr(rank)+"_"+NbrToStr(i)+".bin");
    }

    // Solve
    std::vector<Cplx> sol(nb_dof,0);
    std::vector<double> sol_abs(nb_dof),sol_real(nb_dof);
    htool::DDM<Cplx,htool::partialACA,htool::GeometricClustering> ddm(generator_W,W,ovr_subdomain_to_global,cluster_to_ovr_subdomain,neighbors,intersections);
    opt.parse("-hpddm_schwarz_method none");
    ddm.facto_one_level();
    ddm.solve(rhs.data(),sol.data());
    ddm.print_infos();

    
    if (rank==0){
        htool::vector_to_bytes(sol,"sol.bin");
    }

    // Matrix
    htool::SubMatrix<Cplx> A=generator_W.get_submatrix(tab,tab);
    A.matrix_to_bytes("matrix.bin");


    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}
