#ifndef __QuenchCorr_h_NONEQMSM_wanghan__
#define __QuenchCorr_h_NONEQMSM_wanghan__

#include "Defines.h"
#include "Distribution.h"
#include "Perturbation.h"
#include "Integrator.h"

class QuenchCorr
{
  Distribution_1d corr;
  double sigma;
  double dt;
  const Perturbation * pert;
  double calG (const std::vector<Dofs > & dofs,
	       const std::vector<Dofs > & dw);
  FILE * fp;
public:
  QuenchCorr (const double & x0,
	      const double & x1,
	      const unsigned & nx,
	      const double & v0,
	      const double & v1,
	      const unsigned & nv,
	      const EulerMaruyama & quench_inte0,	// quench integrator at equilibrium
	      const Perturbation & pert);		// perturbation
  ~QuenchCorr ();
  void deposit (const Dofs & endpoint,
		const std::vector<Dofs > & dofs,
		const std::vector<Dofs > & dw);
  void average ();
  const Distribution_1d & getDist () const {return corr;}
}
    ;



#endif
