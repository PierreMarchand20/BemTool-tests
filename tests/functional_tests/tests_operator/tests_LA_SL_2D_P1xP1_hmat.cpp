#include "tests_operator_hmat.hpp"

int main(int argc, char *argv[]) {
  // Initialize the MPI environment
	MPI_Init(&argc,&argv);

	// Get the rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // HTOOL variable
  htool::SetNdofPerElt(1);
  htool::SetEpsilon(1e-6);
  htool::SetEta(1);
  htool::SetMinClusterSize(10);

  // Data
  Real kappa = 0.5;
  Real radius = 0.5;
  Real lc = 0.002;

  // Test
  int test = 0;
  htool::tic();
  Real error = 100*Test_operator_hmat<1,LA_SL_2D_P1xP1,P1_1D>(kappa,radius,lc);
  htool::toc();
	test = test || error>1;
  if (rank==0){
		std::cout << "Relative error:\t"<<error<<" %"<<std::endl;
	}

  // Finalize the MPI environment.
	MPI_Finalize();

  return test;
}
