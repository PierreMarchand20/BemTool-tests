#include "tests_potential_hmat.hpp"

using namespace bemtool;
int main(int argc, char *argv[]) {
  // Initialize the MPI environment
  MPI_Init(&argc,&argv);

  // Get the rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Data
  Real kappa = 0.5;
  Real radius = 1.;
  Real lc = 0.0025;
	Real lc_output = 0.05;
  bool test = 0;

  // HTOOL variable
  htool::SetNdofPerElt(1);
  htool::SetEpsilon(1e-6);
  htool::SetEta(1);
  htool::SetMinClusterSize(10);

  // Test
  htool::tic();
  double error = 100*Test_potential_2D_hmat<LA,2>(kappa, radius, lc, lc_output);
  htool::toc();
  test = test || error>1;
  if (rank){
		std::cout << "Relative error:\t"<<error<<" %"<<std::endl;
	}
  // Finalize the MPI environment.
  MPI_Finalize();

  return test;
}
