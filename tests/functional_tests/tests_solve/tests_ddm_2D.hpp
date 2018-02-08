#ifndef TESTS_SOLVE_HPP
#define TESTS_SOLVE_HPP

#include <bemtool/tools.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>
#include <bemtool-tests/tools.hpp>

using namespace bemtool;

int Test_ddm_HE_2D_P1_direct_dir_hmat(Real kappa, Real radius, Real lc, int save=0) {
  // Get the rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Data
  R3 dir; dir[0]=-sqrt(2)/2.;dir[1]=-sqrt(2)/2.;dir[2]=0;
	std::string outputpath = "../output/tests/performance_tests/tests_solve/ddm_HE_2D_P1_direct_dir_hmat/";

  // Mesh
  if (rank==0){
		gmsh_circle("circle",radius,lc);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  Geometry node("circle.msh");
  Mesh1D mesh; mesh.Load(node,0);
	Orienting(mesh);
	mesh = unbounded;
  int nb_elt = NbElt(mesh);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank==0){
    gmsh_clean("circle");
  }

  // Mesh
	Geometry node_output("../data/output.msh");
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
	if (rank==0 && save>0){
    WritePointValGmsh(mesh_output,(outputpath+"uinc_real.msh").c_str(),uinc_real);
    WritePointValGmsh(mesh_output,(outputpath+"uinc_abs.msh").c_str(),uinc_abs);
  }

  // Dof
  Dof<P1_1D> dof(mesh);
  int nb_dof = NbDof(dof);
  std::vector<htool::R3> x(nb_dof);
  for (int i=0;i<nb_dof;i++){
    x[i][0]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][0];
    x[i][1]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][1];
    x[i][2]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][2];
  }

	// Operator
	Potential<PotKernel<HE,SL_POT,2,P1_1D>> POT_SL(mesh,kappa);
  Potential<PotKernel<HE,DL_POT,2,P1_1D>> POT_DL(mesh,kappa);

  // Generator
  SubBIO_Generator<HE_SL_2D_P1xP1,P1_1D> generator_V(dof,kappa);
	SubBIO_Generator<HE_DL_2D_P1xP1,P1_1D> generator_K(dof,kappa);
	POT_Generator<PotKernel<HE,SL_POT,2,P1_1D>,P1_1D> generator_SL(POT_SL,dof,node_output);
  POT_Generator<PotKernel<HE,DL_POT,2,P1_1D>,P1_1D> generator_DL(POT_DL,dof,node_output);

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
    const N2&         jdof = dof[i];
    const array<2,R3> xdof = dof(i);
    R2x2 M_local = MassP1(mesh[i]);
    C2 Uinc;
    Uinc[0]= exp( iu*kappa*(xdof[0],dir) );
    Uinc[1]= exp( iu*kappa*(xdof[1],dir) );


    for(int k=0;k<2;k++){
      gd[jdof[k]] -= (M_local(k,0)*Uinc[0]+M_local(k,1)*Uinc[1]);
      temp[jdof[k]]= Uinc[k];
    }
  }
  rhs=K*temp;
  std::transform (rhs.begin(), rhs.end(), gd.begin(), rhs.begin(), [](Cplx u,Cplx v){return 0.5*v+u;});

  // Overlap
  std::vector<int> cluster_to_ovr_subdomain;
	std::vector<int> ovr_subdomain_to_global;
	std::vector<int> neighbors;
	std::vector<std::vector<int> > intersections;

	Partition(V.get_MasterOffset_t(), V.get_permt(),dof,cluster_to_ovr_subdomain,ovr_subdomain_to_global,neighbors,intersections);

  // Visu overlap
  std::vector<double> part_overlap(nb_dof,0);
	for (int i=0;i<ovr_subdomain_to_global.size();i++){
		part_overlap[ovr_subdomain_to_global[i]]=1;
	}
	for (int i =0;i<cluster_to_ovr_subdomain.size();i++){
		part_overlap[ovr_subdomain_to_global[cluster_to_ovr_subdomain[i]]]+=1;
	}

  // WritePointValGmsh(dof,("part_ovlerap_"+NbrToStr(rank)+".msh").c_str(),part_overlap);




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

	// Radiated field
	std::vector<Cplx> sol_SL=SL*sol;
	std::vector<Cplx> sol_DL=DL*temp;
	std::vector<double> rad_real(nb_dof_output,0),rad_phase(nb_dof_output,0),rad_abs(nb_dof_output,0);

	for (int i =0;i<nb_dof_output;i++){
		rad_phase[i]=atan2(std::imag(sol_SL[i]-sol_DL[i]+uinc[i]),std::real(sol_SL[i]-sol_DL[i]+uinc[i]));
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
  if (rank==0 && save >0){
		WritePointValGmsh(mesh_output,(outputpath+"rad_phase.msh").c_str(),rad_phase);
    WritePointValGmsh(mesh_output,(outputpath+"rad_real.msh").c_str(),rad_real);
    WritePointValGmsh(mesh_output,(outputpath+"rad_abs.msh").c_str(),rad_abs);
  }

  return 0;
}

int Test_ddm_HE_2D_P1_indirect_dir_hmat(Real kappa, Real radius, Real lc, int save=0) {
  // Get the rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Data
  R3 dir; dir[0]=-sqrt(2)/2.;dir[1]=-sqrt(2)/2.;dir[2]=0;
	std::string outputpath = "../output/tests/performance_tests/tests_solve/ddm_HE_2D_P1_indirect_dir_hmat/";

  // Mesh
  if (rank==0){
		gmsh_circle("circle",radius,lc);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  Geometry node("circle.msh");
  Mesh1D mesh; mesh.Load(node,0);
  Orienting(mesh);
	mesh = unbounded;
  int nb_elt = NbElt(mesh);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank==0 && save>0){
    WriteMeshParaview(mesh,"mesh_paraview.geo");
    gmsh_clean("circle");
  }

  // Mesh
	Geometry node_output("../data/output.msh");
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
	if (rank==0 && save>0){
        WriteMeshParaview(mesh_output,"mesh_output_paraview.geo");
        WritePointValGmsh(mesh_output,(outputpath+"uinc_real.msh").c_str(),uinc_real);
        WritePointValGmsh(mesh_output,(outputpath+"uinc_abs.msh").c_str(),uinc_abs);
    }

  // Dof
  Dof<P1_1D> dof(mesh);
  int nb_dof = NbDof(dof);
  std::vector<htool::R3> x(nb_dof);
  for (int i=0;i<nb_dof;i++){
    x[i][0]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][0];
    x[i][1]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][1];
    x[i][2]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][2];
  }

	// Operator
	Potential<PotKernel<HE,SL_POT,2,P1_1D>> POT_SL(mesh,kappa);

    // Generator
    SubBIO_Generator<HE_SL_2D_P1xP1,P1_1D> generator_V(dof,kappa);
	POT_Generator<PotKernel<HE,SL_POT,2,P1_1D>,P1_1D> generator_SL(POT_SL,dof,node_output);

	// Cluster trees
	std::shared_ptr<htool::Cluster_tree> t_output=std::make_shared<htool::Cluster_tree>(x_output);
  std::shared_ptr<htool::Cluster_tree> t=std::make_shared<htool::Cluster_tree>(x);

  // HMatrix
  htool::HMatrix<htool::partialACA,Cplx> V(generator_V,t,x);
	htool::HMatrix<htool::partialACA,Cplx> SL(generator_SL,t_output,x_output,t,x);
	// Right-hand side
	std::vector<Cplx> rhs(nb_dof,0);
  std::vector<double> rhs_real(nb_dof,0),rhs_abs(nb_dof,0);
  for (int i=0;i<nb_elt;i++){
    const N2&         jdof = dof[i];
    const array<2,R3> xdof = dof(i);
    R2x2 M_local = MassP1(mesh[i]);
    C2 Uinc;
    Uinc[0]= exp( iu*kappa*(xdof[0],dir) );
    Uinc[1]= exp( iu*kappa*(xdof[1],dir) );

    for(int k=0;k<2;k++){
      rhs[jdof[k]] -= (M_local(k,0)*Uinc[0]+M_local(k,1)*Uinc[1]);
    }
  }

  // Overlap
  std::vector<int> cluster_to_ovr_subdomain;
	std::vector<int> ovr_subdomain_to_global;
	std::vector<int> neighbors;
	std::vector<std::vector<int> > intersections;

	Partition(V.get_MasterOffset_t(), V.get_permt(),dof,cluster_to_ovr_subdomain,ovr_subdomain_to_global,neighbors,intersections);

  // Visu overlap
  std::vector<double> part_overlap(nb_dof,0);
	for (int i=0;i<ovr_subdomain_to_global.size();i++){
		part_overlap[ovr_subdomain_to_global[i]]=1;
	}
	for (int i =0;i<cluster_to_ovr_subdomain.size();i++){
		part_overlap[ovr_subdomain_to_global[cluster_to_ovr_subdomain[i]]]+=1;
	}

  // WritePointValGmsh(dof,("part_ovlerap_"+NbrToStr(rank)+".msh").c_str(),part_overlap);




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
  }

	// Radiated field
	std::vector<Cplx> sol_SL=SL*sol;
  std::vector<double> rad_real(nb_dof_output,0),rad_phase(nb_dof_output,0),rad_abs(nb_dof_output,0);

  for (int i =0;i<nb_dof_output;i++){
		rad_phase[i]=atan2(std::imag(sol_SL[i]+uinc[i]),std::real(sol_SL[i]+uinc[i]));
    rad_real[i]=std::real(sol_SL[i]+uinc[i]);
    rad_abs[i]=std::abs(sol_SL[i]+uinc[i]);
  }

	// Save
	V.print_infos();
  SL.print_infos();
  ddm.print_infos();
  V.save_infos((outputpath+"infos_V.txt").c_str());
  SL.save_infos((outputpath+"infos_SL.txt").c_str());
	ddm.save_infos((outputpath+"infos_V.txt").c_str(),std::ios_base::app);

  if (rank==0 && save>0){
      if (save==1){
          std::vector<Real> output_real(nb_dof_output),output_abs(nb_dof_output);
          for (int j=0;j<nb_dof_output;j++){
              output_real[j]=std::real((sol_SL[j]+uinc[j]);
              output_abs[j]=std::abs((sol_SL[j]+uinc[j]);
          }
          WritePointValParaview(mesh_output,(outputpath+"output_paraview_rad_real.scl").c_str(),output_real);
          WriteCaseParaview((outputpath+"output_paraview_rad_real.case").c_str(),"output_paraview_mesh.geo","rad_real","output_paraview_rad_real.scl");
          WritePointValGmsh(mesh_output,(outputpath+"rad_phase.msh").c_str(),rad_phase);
            WritePointValGmsh(mesh_output,(outputpath+"rad_real.msh").c_str(),rad_real);
            WritePointValGmsh(mesh_output,(outputpath+"rad_abs.msh").c_str(),rad_abs);
      }
      else if (save==2){
          int nb = 50;
          std::vector<double> times(nb);
          for (int i=0;i<nb;i++){
              std::vector<Real> output_real(nb_dof_output),output_abs(nb_dof_output);
              for (int j=0;j<nb_dof_output;j++){
                  output_real[j]=std::real((sol_SL[j]+uinc[j])*exp(-pi*iu*i*(2./nb)));
                  output_abs[j]=std::abs((sol_SL[j]+uinc[j])*exp(-pi*iu*i*(2./nb)));
              }
              times[i]=i;
              WritePointValParaview(mesh_output,("rad_real_paraview.scl"+NbrToStr(i)).c_str(),output_real);
      }
      WriteCaseParaview("mesh_output_paraview.case","mesh_output_paraview.geo","rad_real","rad_real_paraview.scl",times);
      WritePointValGmsh(mesh_output,(outputpath+"rad_phase.msh").c_str(),rad_phase);
        WritePointValGmsh(mesh_output,(outputpath+"rad_real.msh").c_str(),rad_real);
        WritePointValGmsh(mesh_output,(outputpath+"rad_abs.msh").c_str(),rad_abs);
    }

  }

  return 0;
}

#endif
