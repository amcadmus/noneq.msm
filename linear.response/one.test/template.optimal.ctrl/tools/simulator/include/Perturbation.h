#ifndef __Perturbation_h_NONEQMSM_wanghan__
#define __Perturbation_h_NONEQMSM_wanghan__

#include "Defines.h"
#include "Force.h"
#include "Poly.h"
#include "Interpolation.h"
#include <vector>
#include <fstream>
#include <string>

using namespace std;

class Perturbation 
{
public:
  virtual double Fe (const double & time) const = 0;
  virtual int numMode () const;
  virtual void FeMode     (const double & time,
			   vector<double > & modes) const;
  virtual void FeModePref (vector<double > & pref ) const;
public:
  virtual void operator () (const Dofs & dofs,
  			    Dofs & pvalue) const = 0;
  virtual void operator () (const Dofs & dofs,
			    const double & time,
			    Dofs & pvalue) const = 0;
}
    ;

class DissipativeFlux 
{
  const Perturbation * pert;
  const Force * force;
  mutable Dofs pvalue;
  mutable vector<double > fvalue;
public:
  DissipativeFlux (const Perturbation * p,
		   const Force * f);
  double operator () (const Dofs & dofs) const;
}
    ;

// V(r) = -r
class PertConstTilt : public Perturbation
{
  double strength;
  double warmTime;
public:
  PertConstTilt (const double & s,
		 const double warmTime = 0.);
public:
  virtual double Fe (const double & time) const;
  virtual int numMode () const;
  virtual void FeMode (const double & time,
		       vector<double > & modes) const;
  virtual void FeModePref (vector<double > & pref ) const;
public:
  virtual void operator () (const Dofs & dofs,
  			    Dofs & pvalue) const;
  virtual void operator () (const Dofs & dofs,
			    const double & time,
			    Dofs & pvalue) const;
}
    ;


class PertConstTiltTable : public Perturbation
{
  PiecewisePoly pwl;
public:
  PertConstTiltTable (const vector<double > & xx,
		      const vector<double > & vv);
public:
  void reinit (const vector<double > & xx,
	       const vector<double > & vv);
  virtual double Fe (const double & time) const;
  virtual int numMode () const;
  virtual void FeMode (const double & time,
		       vector<double > & modes) const;
  virtual void FeModePref (vector<double > & pref ) const;
public:
  virtual void operator () (const Dofs & dofs,
  			    Dofs & pvalue) const;
  virtual void operator () (const Dofs & dofs,
			    const double & time,
			    Dofs & pvalue) const;
}
    ;

// V(r) = 1/(2 \pi \sigma^2) * \exp (-(r - \mu)^2 / (2 \sigma^2))
class PertGaussian : public Perturbation 
{
  double strength;
  double mu;
  double sigma;
  double warmTime;
private:
  double factor0;
  double factor1;  
public:
  PertGaussian (const double s = 1.,
		const double m = 0.,
		const double sig = 1.,
		const double warmTime = 0.);
public:
  virtual double Fe (const double & time) const;
  virtual int numMode () const;
  virtual void FeMode (const double & time,
		       vector<double > & modes) const;
  virtual void FeModePref (vector<double > & pref ) const;
public:
  virtual void operator () (const Dofs & dofs,
  			    Dofs & pvalue) const;
  virtual void operator () (const Dofs & dofs,
			    const double & time,
			    Dofs & pvalue) const;
}
    ;

#endif
