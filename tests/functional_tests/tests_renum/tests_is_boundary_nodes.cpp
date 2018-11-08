#include <bemtool/tools.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>
#include <bemtool-tests/tools.hpp>

using namespace bemtool;
int main(int argc, char *argv[]) {

    // Initialize the MPI environment
    MPI_Init(&argc,&argv);

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

    // Data
    std::string meshname = argv[1];

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

    //
    std::vector<double> output(nb_dof,0);
    std::vector<int> boundary=is_boundary_nodes(dof);

    WritePointValGmsh(dof,"test.msh",boundary);

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}
