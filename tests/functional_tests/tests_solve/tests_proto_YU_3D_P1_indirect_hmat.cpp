
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
    double kappa = StrToNbr<double>(argv[4]);

    // HTOOL variable
    htool::SetNdofPerElt(1);
    htool::SetEpsilon(0.01);
    htool::SetEta(-1);
    // htool::SetMinClusterSize(10);

    // HPDDM verbosity
    HPDDM::Option& opt = *HPDDM::Option::get();
    opt.parse(argc, argv, rank == 0);
    if(rank != 0)
        opt.remove("verbosity");

    // Mesh
    Geometry node(meshname);
    Mesh2D mesh; mesh.Load(node,0);
    Orienting(mesh);
    int nb_elt = NbElt(mesh);

    // Mesh
    Geometry node_output(meshname_output);
    int nb_dof_output = NbNode(node_output);
    Mesh2D mesh_output; mesh_output.Load(node_output,0);
    std::vector<htool::R3> x_output(nb_dof_output);
    std::vector<Cplx> uinc(nb_dof_output);
    std::vector<double> uinc_real(nb_dof_output),uinc_abs(nb_dof_output);
    for (int i=0;i<nb_dof_output;i++){
        x_output[i][0]=node_output[i][0];
        x_output[i][1]=node_output[i][1];
        x_output[i][2]=node_output[i][2];
    }
    if (rank==0 && save>0){
        WriteMeshParaview(mesh_output,(outputpath+"output_paraview_mesh.geo").c_str());
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
  	Potential<PotKernel<YU,DL_POT,3,P1_2D>> POT_DL(mesh,kappa);
    SubBIOp<BIOp<YU_HS_3D_P1xP1>> subBIO(dof,dof,kappa);
    SubBIOp<BIOp<YU_SL_3D_P1xP1>> subBIO_N(dof,dof,kappa);

    // Generator
    BIO_Generator<YU_HS_3D_P1xP1,P1_2D> generator_W(dof,kappa);
    BIO_Generator<YU_SL_3D_P1xP1,P1_2D> generator_V(dof,kappa);
    POT_Generator<PotKernel<YU,DL_POT,3,P1_2D>,P1_2D> generator_DL(POT_DL,dof,node_output);

  	// Cluster trees
  	std::shared_ptr<htool::Cluster_tree> t_output=std::make_shared<htool::Cluster_tree>(x_output);
    std::shared_ptr<htool::Cluster_tree> t=std::make_shared<htool::Cluster_tree>(x);

    // HMatrix
    htool::HMatrix<Cplx,htool::partialACA,htool::GeometricClustering> W(generator_W,t,x);
    htool::HMatrix<Cplx,htool::partialACA,htool::GeometricClustering> DL(generator_DL,t_output,x_output,t,x);

    // Right-hand side
    std::vector<Cplx> rhs(nb_dof,0);
    std::vector<double> rhs_real(nb_dof,0),rhs_abs(nb_dof,0);
    for (int i=0;i<nb_elt;i++){
        const N3&         jdof = dof[i];
        const array<3,R3> xdof = dof(i);
        R3x3 M_local = MassP1(mesh[i]);
        C3 Uinc;
        Uinc[0]= 100*(-xdof[0][0]-1.5)*(-xdof[0][0]-1.5>0);//exp( iu*kappa*(xdof[0],dir) );
        Uinc[1]= 100*(-xdof[1][0]-1.5)*(-xdof[1][0]-1.5>0);;//exp( iu*kappa*(xdof[1],dir) );
        Uinc[2]= 100*(-xdof[2][0]-1.5)*(-xdof[2][0]-1.5>0);

        for(int k=0;k<2;k++){
            rhs[jdof[k]] += (M_local(k,0)*Uinc[0]+M_local(k,1)*Uinc[1]+M_local(k,2)*Uinc[2]);
        }
    }

    // Overlap
    std::vector<int> cluster_to_ovr_subdomain;
    std::vector<int> ovr_subdomain_to_global;
    std::vector<int> neighbors;
    std::vector<std::vector<int> > intersections;

    Partition(t->get_masteroffset(),  t->get_perm(),dof,cluster_to_ovr_subdomain,ovr_subdomain_to_global,neighbors,intersections);


    ///////////// Eigenvalue problem
    // Renumbering
    int n_local = ovr_subdomain_to_global.size();
    int n_inside = cluster_to_ovr_subdomain.size();
    std::vector<int> renum(n_local,-1);
    std::vector<int> renum_to_global(n_local);

    for (int i=0;i<cluster_to_ovr_subdomain.size();i++){
      renum[cluster_to_ovr_subdomain[i]]=i;
      renum_to_global[i]=ovr_subdomain_to_global[cluster_to_ovr_subdomain[i]];
    }
    int count =cluster_to_ovr_subdomain.size();
    for (int i=0;i<n_local;i++){
      if (renum[i]==-1){
        renum[i]=count;
        renum_to_global[count++]=ovr_subdomain_to_global[i];
      }
    }

    // Mass matrix
    htool::Matrix<Cplx> Mi(n_local,n_local);
    htool::Matrix<Cplx> Mbi(n_local,n_local);
    std::map<int,NDofLoc<2>>  Ix;
    for(int k=0; k<renum_to_global.size(); k++){
        const std::vector<N2>& jj = dof.ToElt(renum_to_global[k]);
        for(int l=0; l<jj.size(); l++){
            const N2& j = jj[l]; Ix[j[0]][j[1]] = k;
        }
    }

    for(typename std::map<int,NDofLoc<2>>::iterator itx = Ix.begin(); itx!=Ix.end(); itx++){
        const int&   jx = itx->first;
        const NDofLoc<2>& nx = itx->second;

        mat<2,2,Real> M_local = MassP1(mesh[jx]);
        for(int kx=0; kx<2; kx++){
            if(nx[kx]!=-1){
                for(int ky=0; ky<2; ky++){
                    if(nx[ky]!=-1){
                        Mi(nx[kx],nx[ky]) += M_local(kx,ky);
                    }
                }
            }
        }
    }
    // for (int i=0;i<n_local;i++){
    //     for (int j=0;j<n_local;j++){
    //         Mbi(j,i)=Mi(i,j);
    //     }
    // }

    // if (rank==0){
    //     Mi.matlab_save("Mi.txt");
    //     Mbi.matlab_save("Mbi.txt");
    // }

    // LU facto
    // int lda=n_local;
    // int info;
    // std::vector<int> _ipiv(n_local);
    // HPDDM::Lapack<Cplx>::getrf(&n_local,&n_local,Mi.data(),&lda,_ipiv.data(),&info);
    // HPDDM::Lapack<Cplx>::getrf(&n_local,&n_local,Mbi.data(),&lda,_ipiv.data(),&info);
    // if (rank==0){
    //     Mi.matlab_save("LUi.txt");
    //     Mbi.matlab_save("LUbi.txt");
    // }

    // // Partition of unity
    // std::vector<double> D(n_local);
    // htool::Matrix<Cplx> Di(n_local,n_local);
    // fill(D.begin(),D.begin()+n_inside,1);
    // fill(D.begin()+n_inside,D.end(),0);
    //
    // // Building Ai
    // htool::Matrix<Cplx> Ai(n_local,n_local);
    // subBIO.compute_block(renum_to_global,renum_to_global,Ai);
    //
    // for (int i=0;i<n_local;i++){
    //     for (int j=0;j<n_local;j++){
    //         Ai(i,j)=Ai(i,j)*D[i]*D[j];
    //     }
    // }
    // // if (rank==0)
    // //     Ai.matlab_save("Ai.txt");

    // Building Bi
    htool::Matrix<Cplx> Bi(n_local,n_local);
    subBIO_N.compute_neumann_block(renum_to_global,renum_to_global,Bi);
    // if (rank==0)
    //     Bi.matlab_save("Bi.txt");

    // // M^-1
    // const char l='N';
    // lda=n_local;
    // int ldb=n_local;
    // int nrhs = 1;
    // // std::vector<Cplx> ones (n_local,1);
    // // HPDDM::Lapack<Cplx>::getrs(&l,&n_local,&nrhs,Mi.data(),&lda,_ipiv.data(),ones.data(),&ldb,&info);
    // HPDDM::Lapack<Cplx>::getrs(&l,&n_local,&n_local,Mi.data(),&lda,_ipiv.data(),Ai.data(),&ldb,&info);
    // HPDDM::Lapack<Cplx>::getrs(&l,&n_local,&n_local,Mi.data(),&lda,_ipiv.data(),Bi.data(),&ldb,&info);
    // HPDDM::Lapack<Cplx>::getrs(&l,&n_local,&nrhs,Mi.data(),&lda,_ipiv.data(),ones.data(),&ldb,&info);
    // if (rank==0){
    //     Bi.matlab_save("Ii.txt");
    //     Ai.matlab_save("Ji.txt");
    //     htool::matlab_save(ones,"test.txt");
    // }

    // // Final matrix
    // htool::Matrix<Cplx> Ci(n_local,n_local);
    // Bi.mvprod(Ai.data(),Ci.data(),n_local);
    // // if (rank==0)
    // //     Ci.matlab_save("Ci.txt");

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
    // htool::DDM<Cplx,htool::partialACA,htool::GeometricClustering> ddm(generator_V,V,ovr_subdomain_to_global,cluster_to_ovr_subdomain,neighbors,intersections);
    // ddm.solve(rhs.data(),sol.data());

    htool::Proto_DDM<htool::partialACA,Cplx> ddm(generator_W,W,ovr_subdomain_to_global,cluster_to_ovr_subdomain,neighbors,intersections);
    htool::Proto_HPDDM<htool::partialACA,Cplx> hpddm_op(W,ddm);
    hpddm_op.build_coarse_space(Mi,Bi);
    hpddm_op.facto_one_level();
    hpddm_op.solve(rhs.data(),sol.data());


    // Radiated field
    std::vector<Cplx> sol_DL=DL*sol;


    // Save
    W.print_infos();
    DL.print_infos();
    hpddm_op.print_infos();
    W.save_infos((outputpath+"infos_V.txt").c_str());
    DL.save_infos((outputpath+"infos_SL.txt").c_str());
    // ddm.save_infos((outputpath+"infos_V.txt").c_str(),std::ios_base::app);

    if (rank==0 && save>0){
        if (save==1){
            std::vector<Real> output_real(nb_dof_output),output_abs(nb_dof_output),output_trace_real(nb_dof),output_trace_abs(nb_dof);
            for (int j=0;j<nb_dof;j++){
                output_trace_real[j]=std::real(sol[j]);
                output_trace_abs[j]=std::abs(sol[j]);
            }
            for (int j=0;j<nb_dof_output;j++){
                output_real[j]=std::real(sol_DL[j]);
                output_abs[j]=std::abs(sol_DL[j]);
            }
            WritePointValParaview(mesh_output,(outputpath+"output_paraview_rad_real.scl").c_str(),output_real);
            WriteCaseParaview((outputpath+"output_paraview_rad_real.case").c_str(),"output_paraview_mesh.geo","rad_real","output_paraview_rad_real.scl");
            WritePointValGmsh(mesh_output,(outputpath+"output_gmsh_rad_real.msh").c_str(),output_real);
            WritePointValGmsh(mesh_output,(outputpath+"output_gmsh_rad_abs.msh").c_str(),output_abs);
            WritePointValGmsh(dof,(outputpath+"output_gmsh_trace_real.msh").c_str(),output_trace_real);
            WritePointValGmsh(dof,(outputpath+"output_gmsh_trace_abs.msh").c_str(),output_trace_abs);
        }

    }


    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}
