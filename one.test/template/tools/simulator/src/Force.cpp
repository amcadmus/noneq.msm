#include "Force.h"

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

void DoubleWell::
operator () (const Dofs & dofs,
	     vector<double> & fvalue) const
{
  operator () (dofs, 0, fvalue);
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

SingleWell::
SingleWell (const double k)
    :kk(k)
{
}

void SingleWell::
operator () (const Dofs & dofs,
	     const double & time,
	     vector<double> & fvalue) const
{
  for (unsigned dd = 0; dd < NUMDOFS; ++dd){
    fvalue[dd] = -kk * dofs.xx[dd];
  }
}

void SingleWell::
operator () (const Dofs & dofs,
	     vector<double> & fvalue) const
{
  operator () (dofs, 0, fvalue);
}

double SingleWell::
potential (const Dofs & dofs,
	   const double & time) const
{
  double tmp = 0;
  for (unsigned dd = 0; dd < NUMDOFS; ++dd){
    tmp += dofs.xx[dd] * dofs.xx[dd];
  }
  
  return 0.5 * kk * tmp;
}

