#ifndef __Integrator_NONEQMSM_h_wanghan__
#define __Integrator_NONEQMSM_h_wanghan__

#include "Defines.h"
#include "Force.h"
#include "Perturbation.h"

class Integrator
{
public:
  virtual double getDt    () const =0;
  virtual double getSigma () const =0;
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
  mutable Dofs storedw;
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
  virtual double getDt    () const {return dt;}
  virtual double getSigma () const {return sigma;}
  virtual void step (Dofs & dofs,
		     const double time) const;
  const Dofs & getDw () const {return storedw;}
}
    ;


#endif
