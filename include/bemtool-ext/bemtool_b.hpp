#ifndef TOOLS_HPP
#define TOOLS_HPP

//=====================//

#include "bemtool/equations.hpp"

#include "bemtool/calculus/calculus.hpp"

#include "bemtool/mesh/element.hpp"
#include "bemtool/mesh/adjacency.hpp"
#include "bemtool/mesh/normal.hpp"
#include "bemtool/mesh/mesh.hpp"

#include "bemtool/fem/fem.hpp"
#include "bemtool/fem/femP1.hpp"
#include "bemtool/fem/dof.hpp"
#include "bemtool/fem/shapefct.hpp"
#include "bemtool/fem/interpolation.hpp"

#include "bemtool/operator/operator.hpp"
#include "bemtool/operator/helmholtz_op.hpp"
#include "bemtool/operator/laplace_op.hpp"
#include "bemtool/operator/yukawa_op.hpp"
#include "bemtool/operator/maxwell_op.hpp"

#include "bemtool/potential/potential.hpp"
#include "bemtool/potential/helmholtz_pot.hpp"
#include "bemtool/potential/laplace_pot.hpp"
#include "bemtool/potential/yukawa_pot.hpp"

#include "bemtool/quadrature/dunavant.hpp"
#include "bemtool/quadrature/quad.hpp"
#include "bemtool/quadrature/quad_bem.hpp"
#include "bemtool/quadrature/quad_pot.hpp"

// #include "miscellaneous/eigen_wrap.hpp"
#include "bemtool/miscellaneous/output.hpp"
#include "bemtool/miscellaneous/misc.hpp"
#include "bemtool/miscellaneous/coordinates.hpp"
#include "bemtool/miscellaneous/specialfct.hpp"
#include "bemtool/miscellaneous/refeigenvalue.hpp"
#include "bemtool/miscellaneous/refsol.hpp"

#endif
