
#include "tests_solve_2D.hpp"

using namespace bemtool;
int main(int argc, char *argv[]) {
  // Initialize the MPI environment
  MPI_Init(&argc,&argv);

  // Get the rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Data
  Real kappa = 10;
  Real radius = 1.;
  Real lc = 0.01;
  bool test =0;

  // HTOOL variable
  htool::SetNdofPerElt(1);
  htool::SetEpsilon(1e-6);
  htool::SetEta(1);
  htool::SetMinClusterSize(10);

  // HPDDM verbosity
	HPDDM::Option& opt = *HPDDM::Option::get();
	opt.parse(argc, argv, rank == 0);
	if(rank != 0)
		opt.remove("verbosity");

  // Solve
  htool::tic();
  Test_solve_HE_2D_P1_dir_hmat(kappa, radius, lc);
  htool::toc();

  // Finalize the MPI environment.
  MPI_Finalize();
  return 0;
}
