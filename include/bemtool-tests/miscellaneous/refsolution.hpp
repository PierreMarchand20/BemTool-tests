#ifndef BEMTOOLTESTS_MISC_REFSOLUTION_HPP
#define BEMTOOLTESTS_MISC_REFSOLUTION_HPP

#include <bemtool/equations.hpp>
#include <bemtool/operator/operator.hpp>
#include <bemtool/miscellaneous/specialfct.hpp>


namespace bemtool{

template <int, int>
struct AnalyticalSolution;

template<EquationEnum equation, int dim>
struct RefSolution{

  static inline Cplx Compute(Real r, Real theta, const int& n, const Real& radius=1., const Real& k=1.){
    return AnalyticalSolution<equation,dim>::Compute(r,theta,n,radius,k);
  }

  static inline Cplx ComputeDerivative(Real r, Real theta, const int& n, const Real& radius=1., const Real& k=1.){
    return AnalyticalSolution<equation,dim>::ComputeDerivative(r,theta,n,radius,k);
  }

  static inline Cplx Compute(Real r, Real theta, const N2& nm, const Real& radius=1., const Real& k=1.){
    return AnalyticalSolution<equation,dim>::Compute(r,theta,nm,radius,k);
  }
  static inline Cplx ComputeDerivative(Real r, Real theta, const N2& nm, const Real& radius=1., const Real& k=1.){
    return AnalyticalSolution<equation,dim>::ComputeDerivative(r,theta,nm,radius,k);
  }
};

/*==========
  LAPLACE 2D
  ==========*/
template <> struct AnalyticalSolution<LA,2>{
  static inline Cplx
  Compute(Real r, Real theta, const int& n, const Real& radius=1., const Real& k=1.){
    return pow(r,n)*exp(iu*abs(n)*theta)/pow(radius,n);
  }

  static inline Cplx
  ComputeDerivative(Real r, Real theta, const int& n, const Real& radius=1., const Real& k=1.){
    return n*pow(r,n-1)*exp(iu*abs(n)*theta)/pow(radius,n);
  }
};

/*============
  HELMHOLTZ 2D
  ============*/

template <> struct AnalyticalSolution<HE,2>{
  static inline Cplx
  Compute(Real r, Real theta, const int& n, const Real& radius=1., const Real& k=1.){
    return BesselJ(abs(n),k*r)*exp(iu*abs(n)*theta)/BesselJ(abs(n),k*radius);
  }

  static inline Cplx
  ComputeDerivative(Real r, Real theta, const int& n, const Real& radius=1., const Real& k=1.){
    return (k * DBesselJ_Dx(abs(n),k*radius)*exp(iu*abs(n)*theta))/BesselJ(abs(n),k*radius);
  }
};

/*============
  Modified Helmholtz 2D
  ============*/

template <> struct AnalyticalSolution<YU,2>{
  static inline Cplx
  Compute(Real r, Real theta, const int& n, const Real& radius=1., const Real& k=1.){
    Real r2=r*r;
    return Modified_BesselI(abs(n),k*r)*exp(iu*abs(n)*theta)/Modified_BesselI(abs(n),k*radius);
  }
  static inline Cplx
  ComputeDerivative(Real r, Real theta, const int& n, const Real& radius=1., const Real& k=1.){
    return (k * DModified_BesselI_Dx(abs(n),k*radius)*exp(iu*abs(n)*theta))/Modified_BesselI(abs(n),k*radius);
  }
};

}
#endif
