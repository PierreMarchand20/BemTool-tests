#ifndef TESTS_WEIGHTS
#define TESTS_WEIGHTS

#include <bemtool/tools.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>
#include <bemtool-tests/tools.hpp>

using namespace bemtool;

int main(int argc, char *argv[])
{
    ////////////////////////////////////////////  MPI
    // Initialize the MPI environment
    MPI_Init(&argc,&argv);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Get the number of the process
    int sizeWorld;
    MPI_Comm_size(MPI_COMM_WORLD,&sizeWorld);

    ////////////////////////////////////////////  Inputs
    HPDDM::Option& opt = *HPDDM::Option::get();
    opt.parse(argc, argv, rank == 0,{
        std::forward_as_tuple("meshinput=<input_file>", "input mesh file", HPDDM::Option::Arg::argument),
    });
    
    if(rank != 0)
        opt.remove("verbosity");

    std::string meshname,outputfile;
    if(opt.prefix("meshinput").size())
        meshname=opt.prefix("meshinput");
    else{
        std::cout << "Missing meshinput"<<std::endl;
        return 1;
    }

    // Mesh
    if (rank==0)
        std::cout << "Loading mesh" << std::endl;
    Geometry node(meshname);
    Mesh<2> mesh; mesh.Load(node,0);
    Orienting(mesh);
    int nb_elt = NbElt(mesh);

    // Dof
    if (rank==0)
       std::cout << "Creating dof" << std::endl;
    Dof<P1_2D> dof(mesh);
    int nb_dof = NbDof(dof);
    std::vector<htool::R3> x(nb_dof);
    for (int i=0;i<nb_dof;i++){
      x[i][0]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][0];
      x[i][1]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][1];
      x[i][2]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][2];
    }

    // Weights
    if (rank==0)
        std::cout<< "Compute weights"<<std::endl;
    std::vector<double> g(nb_dof,0);
    for (int i=0;i<nb_dof;i++){
        std::vector<N2> elts= dof.ToElt(i);
        for (int j=0;j<elts.size();j++){
            g[i]+=Vol((mesh)[j]);
            
        }
    }

    // Cluster trees
    if (rank==0)
       std::cout << "Creating cluster tree" << std::endl;
    std::shared_ptr<htool::Cluster_tree> t=std::make_shared<htool::Cluster_tree>(x);
    std::vector<int> labels=t->get_labels(4);
    t->print_infos();

    // Output
    WritePointValGmsh(dof,"output.msh",labels);

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}


#endif