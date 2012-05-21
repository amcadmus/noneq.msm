#include "Force.h"

PertConstTilt::
PertConstTilt (const double & s)
    : strength(s)
{
}

void PertConstTilt::
operator () (const Dofs & dofs,
	     const double & time,
	     Dofs & pvalue) const
{
  pvalue.vv[0] = strength;
}

DoubleWell::
DoubleWell (const double & k,
	    const double & a)
    : kk(k), aa(a)
{
}

void DoubleWell::
operator () (const Dofs & dofs,
	     const double & time,
	     vector<double> & fvalue) const
{
  double tmp = 0;
  for (unsigned dd = 0; dd < NUMDOFS; ++dd){
    tmp += dofs.xx[dd] * dofs.xx[dd];
  }
  tmp -= aa * aa;

  for (unsigned dd = 0; dd < NUMDOFS; ++dd){
    fvalue[dd] = -2. * kk * tmp * dofs.xx[dd];
  }
}

double DoubleWell::
potential (const Dofs & dofs,
	   const double & time) const
{
  double tmp = 0;
  for (unsigned dd = 0; dd < NUMDOFS; ++dd){
    tmp += dofs.xx[dd] * dofs.xx[dd];
  }
  tmp -= aa * aa;
  
  return 0.5 * kk * tmp * tmp;
}



