#include "Perturbation.h"

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
  pvalue.vv[0] = 1.0;
}

double PertConstTilt::
Fe (const double & time) const
{
  double tmp = strength;
  if (warmTime > 0.0 && time < warmTime){
    tmp *= time / warmTime;
  }
  return tmp;
}


PertGaussian::
PertGaussian (const double s,
	      const double m,
	      const double sig,
	      const double w)
    : strength(s), mu(m), sigma(sig), warmTime(w)
{
  factor0 = 1./(sigma * sqrt(2. * M_PI));
  factor1 = 1./(2. * sigma * sigma);
}

double PertGaussian::
Fe (const double & time) const
{
  double tmp = strength;
  if (warmTime > 0.0 && time < warmTime){
    tmp *= time / warmTime;
  }
  return tmp;
}

void PertGaussian::
operator () (const Dofs & dofs,
	     Dofs & pvalue) const
{
  double sum = 0.;
  for (unsigned dd = 0; dd < NUMDOFS; ++dd){
    sum += dofs.xx[dd] * dofs.xx[dd];
  }
  sum = sqrt(sum);
  double pre = factor0 * exp ( - (sum - mu) * (sum - mu) * factor1) * 2. * (sum - mu) * factor1 / sum;
  for (unsigned dd = 0; dd < NUMDOFS; ++dd){
    pvalue.vv[dd] = dofs.xx[dd] * pre;
  }
}

void PertGaussian::
operator () (const Dofs & dofs,
	     const double & time,
	     Dofs & pvalue) const
{
  (*this) (dofs, pvalue);
  pvalue.vv[0] *= Fe(time);
  // for (unsigned dd = 0; dd < NUMDOFS; ++dd){
  //   pvalue.vv[dd] *= Fe(time);
  // }
}


