#include <bemtool/tools.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>
#include <bemtool-tests/tools.hpp>
#include "bemtool-tests/miscellaneous/view.hpp"
#include <htool/visu/view.hpp>

using namespace htool;
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
        /*# Init #*/
        int rankWorld, sizeWorld;
        MPI_Comm_size(MPI_COMM_WORLD, &sizeWorld);
        MPI_Comm_rank(MPI_COMM_WORLD, &rankWorld);

        Scene s;

    	s.init();

    	bemtool::attach_ui(s);

    	statics& gv = Scene::gv;

    	s.run();

        MPI_Finalize();
        return 0;
    return 0;
}
