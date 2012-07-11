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

// V(r) = -r
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
  virtual void operator () (const Dofs & dofs,
  			    Dofs & pvalue) const;
  virtual double Fe (const double & time) const;
  virtual void operator () (const Dofs & dofs,
			    const double & time,
			    Dofs & pvalue) const;
}
    ;


// U(r) = 0.5 * kk * (r^2 - aa^2)^2
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

// U(r) = 0.5 * kk * r^2
class SingleWell : public Force 
{
  double kk;
public:
  public:
  SingleWell (const double k = 8.);
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
