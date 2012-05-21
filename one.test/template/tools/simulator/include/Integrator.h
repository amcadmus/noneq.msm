#ifndef __Integrator_NONEQMSM_h_wanghan__
#define __Integrator_NONEQMSM_h_wanghan__

#include "Defines.h"
#include "Force.h"

class Integrator
{
public:
  virtual void step (Dofs & dofs,
		     const double time) const = 0;
}
    ;

class EulerMaruyama : public Integrator
{
  double dt;
  double sqrtdt;
  double gamma;
  double sigma;
  const Perturbation * pert;
  const Force * force;
public:
  EulerMaruyama ();
  EulerMaruyama (const double & gamma,
		 const double & kT,
		 const double & dt,
		 const Perturbation * p = NULL,
		 const Force * f = NULL);
public:
  void reinit (const double & gamma,
	       const double & kT,
	       const double & dt,
	       const Perturbation * p = NULL,
	       const Force * f = NULL);
  virtual void step (Dofs & dofs,
		     const double time) const;
}
    ;


#endif
