#include <bemtool/tools.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>
#include <bemtool-tests/tools.hpp>

using namespace bemtool;
int main(int argc, char *argv[]) {




// Mesh
gmsh_circle("circle",1,0.03);

Geometry node("circle.msh");
Mesh1D mesh; mesh.Load(node,0);
  Orienting(mesh);
  mesh = unbounded;
int nb_elt = NbElt(mesh);

double kappa = 1;
 BIOp<LA_DL_2D_P1xP1> V(mesh,mesh,kappa);
 Dof<P1_1D> dof(mesh);
 int nb_dof = NbDof(dof);

 std::vector<Cplx> test(nb_dof);

#pragma omp parallel for schedule(guided)
for (int i=0;i<nb_dof;i++){
    test[0]=V(dof.ToElt(i),dof.ToElt(i));
    std::cout<<test[i]<<std::endl;
}


}
