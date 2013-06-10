#ifndef __Force_h_NONEQMSM_wanghan__
#define __Force_h_NONEQMSM_wanghan__

#include "Defines.h"

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
