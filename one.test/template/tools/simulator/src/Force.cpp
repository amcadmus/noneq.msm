#include "Force.h"

DissipativeFlux::
DissipativeFlux (const Perturbation * p,
		 const Force * f)
    :pert (p), force(f), fvalue(NUMDOFS)
{
}

double DissipativeFlux::
operator () (const Dofs & dof) const
{
  (*pert)  (dof, pvalue);
  (*force) (dof, fvalue);
  double tmp = 0.;
  for (unsigned dd = 0; dd < NUMDOFS; ++dd){
    tmp += pvalue.xx[dd] * fvalue[dd] - pvalue.vv[dd] * dof.vv[dd];
  }
  return tmp;
}


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
  (*this) (dofs, pvalue);
  pvalue.vv[0] *= Fe(time);
}

void PertConstTilt::
operator () (const Dofs & dofs,
	     Dofs & pvalue) const
{
  pvalue.vv[0] = strength;
}

double PertConstTilt::
Fe (const double & time) const
{
  double tmp = 1.;
  if (time < warmTime){
    tmp = time / warmTime;
  }
  return tmp;
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

