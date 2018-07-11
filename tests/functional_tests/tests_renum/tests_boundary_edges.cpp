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
    std::vector<std::pair<int,int>> boundary_edges=is_boundary(dof);
    std::cout <<boundary_edges.size()<<std::endl;
    std::vector<double> output(nb_dof,0);

    std::unordered_set<int> s;
    for (int i=0;i<boundary_edges.size();i++){
        s.insert(boundary_edges[i].first);
        s.insert(boundary_edges[i].second);
    }
    std::vector<int> boundary(s.size());
    boundary.assign( s.begin(), s.end() );

    std::cout <<boundary.size()<<std::endl;

    for (int i=0;i<boundary.size();i++){
        std::cout << boundary[i]<< std::endl;
        output[boundary[i]]=1;
    }
    WritePointValGmsh(dof,"test.msh",output);

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}
