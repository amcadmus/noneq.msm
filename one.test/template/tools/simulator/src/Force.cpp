#include "Force.h"

PertConstTilt::
PertConstTilt (const double & s,
	       const double w)
    : strength(s), warmTime(w)
{
}

void PertConstTilt::
operator () (const Dofs & dofs,
	     const double & time,
	     Dofs & pvalue) const
{
  double tmp = strength;
  if (time < warmTime){
    tmp = tmp * time / warmTime;
  }
  pvalue.vv[0] = tmp;
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



