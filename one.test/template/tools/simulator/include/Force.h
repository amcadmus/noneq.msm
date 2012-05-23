#ifndef __Force_h_NONEQMSM_wanghan__
#define __Force_h_NONEQMSM_wanghan__

#include "Defines.h"

class Perturbation 
{
public:
  virtual void operator () (const Dofs & dofs,
			    Dofs & pvalue) const = 0;
  virtual double Fe (const double & time) const = 0;
  virtual void operator () (const Dofs & dofs,
			    const double & time,
			    Dofs & pvalue) const = 0;
}
    ;

class Force 
{
public:
  virtual double potential (const Dofs & dofs,
			    const double & time) const = 0;
  virtual void operator () (const Dofs & dofs,
			    vector<double> & fvalue) const = 0;
  virtual void operator () (const Dofs & dofs,
			    const double & time,
			    vector<double> & fvalue) const = 0;
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


class PertConstTilt : public Perturbation
{
  double strength;
  double warmTime;
public:
  PertConstTilt (const double & s,
		 const double warmTime = 0.);
  virtual void operator () (const Dofs & dofs,
			    Dofs & pvalue) const;
  virtual double Fe (const double & time) const;
  virtual void operator () (const Dofs & dofs,
			    const double & time,
			    Dofs & pvalue) const;
}
    ;

class DoubleWell : public Force
{
  double kk;
  double aa;
public:
  DoubleWell (const double & k,
	      const double & a);
  virtual double potential (const Dofs & dofs,
			    const double & time) const;
  virtual void operator () (const Dofs & dofs,
			    vector<double> & fvalue) const ;
  virtual void operator () (const Dofs & dofs,
			    const double & time,
			    vector<double> & fvalue) const ;
}
    ;


#endif
