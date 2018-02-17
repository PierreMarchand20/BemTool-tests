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
    if (argc < 6) {
        // Tell the user how to run the program
        std::cerr << "Usage: " << argv[0] << " meshname meshname_output outputpath frequency save" << std::endl;
        /* "Usage messages" are a conventional way of telling the user
         * how to run a program if they enter the command incorrectly.
         */
        return 1;
    }

    // Data
    R3 dir; dir[0]=1;dir[1]=0;dir[2]=0;
    std::string meshname = argv[1];
    std::string meshname_output = argv[2];
    std::string outputpath = argv[3];
    int save = StrToNbr<int>(argv[5]); // 0 no save, 1 sol, 2 video
    double frequency = StrToNbr<double>(argv[4])*1e9;
    double c0= 299792548;
    double kappa = 2*pi*frequency/c0;


    // HTOOL variable
    htool::SetNdofPerElt(1);
    htool::SetEpsilon(1e-3);
    htool::SetEta(1);
    htool::SetMinClusterSize(10);

    // HPDDM verbosity
    HPDDM::Option& opt = *HPDDM::Option::get();
    opt.parse(argc, argv, rank == 0);
    if(rank != 0)
    	opt.remove("verbosity");

    // Mesh
    Geometry node(meshname);
    Mesh2D mesh; mesh.Load(node,0);
    Orienting(mesh);
    mesh = unbounded;
    int nb_elt = NbElt(mesh);

    // Mesh
    Geometry node_output(meshname_output);
    int nb_dof_output = NbNode(node_output);
    Mesh2D mesh_output; mesh_output.Load(node_output,0);
    // Orienting(mesh_output);
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

    // Dof
    Dof<P1_2D> dof(mesh);
    int nb_dof = NbDof(dof);
    std::vector<htool::R3> x(nb_dof);
    for (int i=0;i<nb_dof;i++){
      x[i][0]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][0];
      x[i][1]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][1];
      x[i][2]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][2];
    }

    // Operator
    Potential<PotKernel<HE,SL_POT,3,P1_2D>> POT_SL(mesh,kappa);
    Potential<PotKernel<HE,DL_POT,3,P1_2D>> POT_DL(mesh,kappa);

    // Generator
    SubBIO_Generator<HE_SL_3D_P1xP1,P1_2D> generator_V(dof,kappa);
    SubBIO_Generator<HE_DL_3D_P1xP1,P1_2D> generator_K(dof,kappa);
    POT_Generator<PotKernel<HE,SL_POT,3,P1_2D>,P1_2D> generator_SL(POT_SL,dof,node_output);
    POT_Generator<PotKernel<HE,DL_POT,3,P1_2D>,P1_2D> generator_DL(POT_DL,dof,node_output);

    // Cluster trees
    std::shared_ptr<htool::Cluster_tree> t_output=std::make_shared<htool::Cluster_tree>(x_output);
    std::shared_ptr<htool::Cluster_tree> t=std::make_shared<htool::Cluster_tree>(x);

    // HMatrix
    htool::HMatrix<htool::partialACA,Cplx> V(generator_V,t,x);
    htool::HMatrix<htool::partialACA,Cplx> K(generator_K,t,x);
    htool::HMatrix<htool::partialACA,Cplx> SL(generator_SL,t_output,x_output,t,x);
    htool::HMatrix<htool::partialACA,Cplx> DL(generator_DL,t_output,x_output,t,x);

    // Right-hand side
    std::vector<Cplx> temp(nb_dof,0),gd(nb_dof,0),rhs(nb_dof,0);
    std::vector<double> rhs_real(nb_dof,0),rhs_abs(nb_dof,0),gd_real(nb_dof,0),gd_abs(nb_dof,0);
    for (int i=0;i<nb_elt;i++){
        const N3&         jdof = dof[i];
        const array<3,R3> xdof = dof(i);
        R3x3 M_local = MassP1(mesh[i]);
        C3 Uinc;
        Uinc[0]= exp( iu*kappa*(xdof[0],dir) );
        Uinc[1]= exp( iu*kappa*(xdof[1],dir) );
        Uinc[2]= exp( iu*kappa*(xdof[2],dir) );

        for(int k=0;k<3;k++){
            gd[jdof[k]] -= (M_local(k,0)*Uinc[0]+M_local(k,1)*Uinc[1]+M_local(k,2)*Uinc[2]);
            temp[jdof[k]]= Uinc[k];
        }
    }
    rhs=K*temp;
    std::transform (rhs.begin(), rhs.end(), gd.begin(), rhs.begin(), [](Cplx u,Cplx v){return u+0.5*v;});

    // Overlap
    std::vector<int> cluster_to_ovr_subdomain;
    std::vector<int> ovr_subdomain_to_global;
    std::vector<int> neighbors;
    std::vector<std::vector<int> > intersections;

    Partition(V.get_MasterOffset_t(), V.get_permt(),dof,cluster_to_ovr_subdomain,ovr_subdomain_to_global,neighbors,intersections);

    // Visu overlap
    if (save>0){
        std::vector<double> part_overlap(nb_dof,0);
        for (int i=0;i<ovr_subdomain_to_global.size();i++){
            part_overlap[ovr_subdomain_to_global[i]]=1;
        }
        for (int i =0;i<cluster_to_ovr_subdomain.size();i++){
            part_overlap[ovr_subdomain_to_global[cluster_to_ovr_subdomain[i]]]+=1;
        }
        WritePointValGmsh(dof,("part_ovlerap_"+NbrToStr(rank)+".msh").c_str(),part_overlap);
    }



    // Solve
    std::vector<Cplx> sol(nb_dof,0);
    std::vector<double> sol_abs(nb_dof),sol_real(nb_dof);
    htool::DDM<htool::partialACA,Cplx> ddm(generator_V,V,ovr_subdomain_to_global,cluster_to_ovr_subdomain,neighbors,intersections);

    ddm.solve(rhs.data(),sol.data());

    for (int i=0;i<nb_dof;i++){
        sol_abs[i]=std::abs(sol[i]);
        sol_real[i]=std::real(sol[i]);
        rhs_abs[i]=std::abs(rhs[i]);
        rhs_real[i]=std::real(rhs[i]);
        gd_abs[i]=std::abs(gd[i]);
        gd_real[i]=std::real(gd[i]);
    }


    // Save
    if (rank==0 && save>0){
        WritePointValGmsh(dof,(outputpath+"sol_abs.msh").c_str(),sol_abs);
        WritePointValGmsh(dof,(outputpath+"sol_real.msh").c_str(),sol_real);
        WritePointValGmsh(dof,(outputpath+"rhs_real.msh").c_str(),rhs_real);
        WritePointValGmsh(dof,(outputpath+"rhs_abs.msh").c_str(),rhs_abs);

        WritePointValGmsh(dof,(outputpath+"gd_real.msh").c_str(),gd_real);
        WritePointValGmsh(dof,(outputpath+"gd_abs.msh").c_str(),gd_abs);
    }

    // Radiated field
    std::vector<Cplx> sol_SL=SL*sol;
    std::vector<Cplx> sol_DL=DL*temp;
    std::vector<double> rad_real(nb_dof_output,0),rad_abs(nb_dof_output,0);

    for (int i =0;i<nb_dof_output;i++){
        rad_real[i]=std::real(sol_SL[i]-sol_DL[i]+uinc[i]);
        rad_abs[i]=std::abs(sol_SL[i]-sol_DL[i]+uinc[i]);
    }

    // Save
    V.print_infos();
    K.print_infos();
    SL.print_infos();
    DL.print_infos();
    V.save_infos((outputpath+"infos_V.txt").c_str());
    K.save_infos((outputpath+"infos_K.txt").c_str());
    SL.save_infos((outputpath+"infos_SL.txt").c_str());
    DL.save_infos((outputpath+"infos_DL.txt").c_str());
    ddm.save_infos((outputpath+"infos_V.txt").c_str(),std::ios_base::app);
    if (rank==0 && save>0){
        WritePointValGmsh(mesh_output,(outputpath+"rad_real.msh").c_str(),rad_real);
        WritePointValGmsh(mesh_output,(outputpath+"rad_abs.msh").c_str(),rad_abs);
    }



    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}
