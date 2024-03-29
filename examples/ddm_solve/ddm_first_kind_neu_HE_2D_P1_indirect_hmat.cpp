
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

    ////////////////////////////////////////////  Inputs
    HPDDM::Option& opt = *HPDDM::Option::get();
    opt.parse(argc, argv, rank == 0,{
        std::forward_as_tuple("meshinput=<input_file>", "input mesh file", HPDDM::Option::Arg::argument),
        std::forward_as_tuple("meshoutput=<output_file>", "output mesh file", HPDDM::Option::Arg::argument),
        std::forward_as_tuple("outputpath=<output>", "output path", HPDDM::Option::Arg::argument),
        std::forward_as_tuple("kappa=<1>", "kappa.", HPDDM::Option::Arg::numeric),
        std::forward_as_tuple("overlap=<1>", "overlap.", HPDDM::Option::Arg::integer),
        std::forward_as_tuple("save=<0>", "integer", HPDDM::Option::Arg::integer),
        std::forward_as_tuple("epsilon=<0.01>", "epsilon.", HPDDM::Option::Arg::numeric),
        std::forward_as_tuple("eta=<10>", "eta.", HPDDM::Option::Arg::numeric),
        std::forward_as_tuple("mintargetdepth=<1>", "mintargetdepth.", HPDDM::Option::Arg::numeric),
        std::forward_as_tuple("minsourcedepth=<1>", "minsourcedepth.", HPDDM::Option::Arg::numeric),
    });
    if(rank != 0)
        opt.remove("verbosity");

    std::string meshname,meshname_output,outputpath;
    // Meshname
    if(opt.prefix("meshinput").size())
        meshname=opt.prefix("meshinput");
    else{
        std::cout << "Missing meshinput"<<std::endl;
        return 1;
    }
    // Meshname output
    if(opt.prefix("meshoutput").size())
        meshname_output=opt.prefix("meshoutput");
    else{
        std::cout << "Missing meshoutput"<<std::endl;
        return 1;
    }
    // Output path
    if(opt.prefix("outputpath").size()){
        outputpath=opt.prefix("outputpath");
    }
    else{
        std::cout << "Missing outputpath"<<std::endl;
        return 1;
    }

    // Save
    int save =opt.app()["save"];

    // kappa
    double kappa = opt.app()["kappa"];

    // other datas
    R3 dir; dir[0]=1./std::sqrt(2);dir[1]=1./std::sqrt(2);dir[2]=0;

    //////////////////////////////////////////// HTOOL variable
    htool::SetNdofPerElt(1);
    htool::SetEpsilon(opt.app()["epsilon"]);
    htool::SetEta(opt.app()["eta"]);
    htool::SetMinTargetDepth(opt.app()["mintargetdepth"]);
    htool::SetMinSourceDepth(opt.app()["minsourcedepth"]);

    ////////////////////////////////////////////     Problem setup
    // Mesh
    if (rank==0)
        std::cout << "Loading mesh" << std::endl;
    Geometry node(meshname);
    Mesh1D mesh; mesh.Load(node);
    Orienting(mesh);
    mesh = unbounded;
    int nb_elt = NbElt(mesh);
    std::vector<R3> normal = NormalTo(mesh);

    // Mesh
    Geometry node_output(meshname_output);
    int nb_dof_output = NbNode(node_output);
    Mesh2D mesh_output; mesh_output.Load(node_output);
    std::vector<htool::R3> x_output(nb_dof_output);
    std::vector<Cplx> uinc(nb_dof_output);
    std::vector<double> uinc_real(nb_dof_output),uinc_abs(nb_dof_output);
    for (int i=0;i<nb_dof_output;i++){
        x_output[i][0]=node_output[i][0];
        x_output[i][1]=node_output[i][1];
        x_output[i][2]=node_output[i][2];
        double temp = dir[0]*x_output[i][0]+dir[1]*x_output[i][1]+dir[2]*x_output[i][2];
        uinc[i] = exp(iu*kappa*temp);
        uinc_real[i] = std::real(uinc[i]);
        uinc_abs[i] = std::abs(uinc[i]);
    }
    if (rank==0 && save>0){
        WritePointValGmsh(mesh_output,(outputpath+"uinc_real.msh").c_str(),uinc_real);
        WritePointValGmsh(mesh_output,(outputpath+"uinc_abs.msh").c_str(),uinc_abs);
    }

    // Dof
    if (rank==0)
       std::cout << "Creating dof" << std::endl;
    Dof<P1_1D> dof(mesh);
    int nb_dof = NbDof(dof);
    std::vector<htool::R3> x(nb_dof);
    for (int i=0;i<nb_dof;i++){
      x[i][0]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][0];
      x[i][1]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][1];
      x[i][2]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][2];
    }

  	// Operator
  	Potential<PotKernel<HE,DL_POT,2,P1_1D>> POT_DL(mesh,kappa);

    // Generator
    if (rank==0)
       std::cout << "Creating generators"<<std::endl;
    BIO_Generator<HE_HS_2D_P1xP1,P1_1D> generator_W(dof,kappa);
    POT_Generator<PotKernel<HE,DL_POT,2,P1_1D>,P1_1D> generator_DL(POT_DL,dof,node_output);

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

    // HMatrix
    if (rank==0)
	   std::cout << "Building Hmatrix" << std::endl;
    htool::HMatrix<Cplx,htool::partialACA,htool::GeometricClustering> W(generator_W,x);
    htool::HMatrix<Cplx,htool::partialACA,htool::GeometricClustering> DL(generator_DL,x_output,x);

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

    // Overlap
    if (rank==0)
        std::cout << "Building partitions" << std::endl;
    std::vector<int> cluster_to_ovr_subdomain;
    std::vector<int> ovr_subdomain_to_global;
    std::vector<int> neighbors;
    std::vector<std::vector<int> > intersections;

    Partition(W.get_MasterOffset_t(), W.get_permt(),dof,cluster_to_ovr_subdomain,ovr_subdomain_to_global,neighbors,intersections);

    // Visu overlap
    if (save>0){
        std::vector<double> part_overlap(nb_dof,0);
      	for (int i=0;i<ovr_subdomain_to_global.size();i++){
      		part_overlap[ovr_subdomain_to_global[i]]=1;
      	}
      	for (int i =0;i<cluster_to_ovr_subdomain.size();i++){
      		part_overlap[ovr_subdomain_to_global[cluster_to_ovr_subdomain[i]]]+=1;
      	}
        WritePointValGmsh(dof,(outputpath+"part_ovlerap_"+NbrToStr(rank)+".msh").c_str(),part_overlap);
    }

    // Solve
    std::vector<Cplx> sol(nb_dof,0);
    std::vector<double> sol_abs(nb_dof),sol_real(nb_dof);
    htool::DDM<Cplx,htool::partialACA,htool::GeometricClustering> ddm(generator_W,W,ovr_subdomain_to_global,cluster_to_ovr_subdomain,neighbors,intersections);
    opt.parse("-hpddm_schwarz_method asm");
    ddm.facto_one_level();
    ddm.solve(rhs.data(),sol.data());


    // Radiated field
    if (rank==0)
        std::cout << "Radiated field" << std::endl;
    std::vector<Cplx> sol_DL=DL*sol;
    std::vector<double> rad_real(nb_dof_output,0),rad_phase(nb_dof_output,0),rad_abs(nb_dof_output,0);

    for (int i =0;i<nb_dof_output;i++){
        rad_phase[i]=atan2(std::imag(sol_DL[i]+uinc[i]),std::real(sol_DL[i]+uinc[i]));
        rad_real[i]=std::real(sol_DL[i]+uinc[i]);
        rad_abs[i]=std::abs(sol_DL[i]+uinc[i]);
    }

    // Save
    W.print_infos();
    DL.print_infos();
    ddm.print_infos();
    W.save_infos((outputpath+"infos_W.txt").c_str());
    DL.save_infos((outputpath+"infos_DL.txt").c_str());
    ddm.save_infos((outputpath+"solve.txt").c_str());

    if (rank==0 && save>0){
        WritePointValGmsh(mesh_output,(outputpath+"rad_phase.msh").c_str(),rad_phase);
        WritePointValGmsh(mesh_output,(outputpath+"rad_real.msh").c_str(),rad_real);
        WritePointValGmsh(mesh_output,(outputpath+"rad_abs.msh").c_str(),rad_abs);
    }

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}
