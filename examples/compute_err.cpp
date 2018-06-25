#include <bemtool/tools.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>
#include <bemtool-tests/tools.hpp>

using namespace bemtool;

int main(int argc, char *argv[]) {
    // Initialize the MPI environment
    MPI_Init(&argc,&argv);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
    double kappa = 0;


    // HTOOL variable
    htool::SetNdofPerElt(1);
    htool::SetEpsilon(0.001);
    htool::SetEta(100);
    htool::SetMinClusterSize(10);

    // Mesh
    Geometry node(meshname);
    Mesh2D mesh; mesh.Load(node,0);
    Orienting(mesh);
    mesh = unbounded;
    int nb_elt = NbElt(mesh);

    // Dof
    Dof<P1_2D> dof(mesh);
    int nb_dof = NbDof(dof);
    std::vector<htool::R3> x(nb_dof);
    for (int i=0;i<nb_dof;i++){
      x[i][0]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][0];
      x[i][1]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][1];
      x[i][2]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][2];
    }

    // Generator
    BIO_Generator<LA_SL_3D_P1xP1,P1_2D> generator(dof,kappa);

  	// Cluster trees
    std::shared_ptr<htool::Cluster_tree> t=std::make_shared<htool::Cluster_tree>(x);
    t->print_infos();

    // HMatrix
    htool::HMatrix<htool::partialACA,Cplx> A(generator,t,x);
    A.print_infos();
    // Compute error
    if (rank==0) {
        std::cout << "Error : "<<Frobenius_absolute_error(A,generator)<< std::endl;
    }


    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}
