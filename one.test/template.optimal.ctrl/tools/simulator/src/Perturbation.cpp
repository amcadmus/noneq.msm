#include "Perturbation.h"
#include <algorithm>

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

int Perturbation::
numMode () const
{
  return 1;
}

void Perturbation::
FeMode (const double & time,
	vector<double > & modes) const
{
  modes.resize(1);
  modes[0] = Fe(time);
}

void Perturbation::
FeModePref (vector<double > & perf) const
{
  perf.resize(1);
  perf[0] = 1;
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

int PertConstTilt::
numMode () const
{
  return 1;
}

void PertConstTilt::
FeMode (const double & time,
	vector<double > & modes) const
{
  modes.resize(1);
  modes[0] = 1.;
  if (warmTime > 0.0 && time < warmTime){
    modes[0] *= time / warmTime;
  }
}

void PertConstTilt::
FeModePref (vector<double > & pref ) const
{
  pref.resize(1);
  pref[0] = strength;
}




PertConstTiltTable::
PertConstTiltTable (const double & s,
		    const vector<double > & xx,
		    const vector<double > & vv)
{
  reinit (xx, vv);
}

void PertConstTiltTable::
reinit (const vector<double > & xx,
	const vector<double > & vv)
{
  Interpolation::piecewiseLinear (xx, vv, pwl);
}


void PertConstTiltTable::
operator () (const Dofs & dofs,
	     const double & time,
	     Dofs & pvalue) const
{
  (*this) (dofs, pvalue);
  pvalue.vv[0] *= Fe(time);
}

void PertConstTiltTable::
operator () (const Dofs & dofs,
	     Dofs & pvalue) const
{
  pvalue.vv[0] = 1.0;
}

double PertConstTiltTable::
Fe (const double & time) const
{
  return pwl.value (time);
}

int PertConstTiltTable::
numMode () const
{
  return pwl.get_x().size();
}

void PertConstTiltTable::
FeMode (const double & time,
	vector<double > & modes) const
{
  modes.resize(numMode());
  std::fill (modes.begin(), modes.end(), 0.0);
  int index = int(time / (pwl.get_x()[1] - pwl.get_x()[0]));
  if (index == numMode() - 1) index --;
  if (index >= 0 && index < numMode() - 1){
    modes[index  ] = (pwl.get_x()[index+1] - time) / (pwl.get_x()[index+1] - pwl.get_x()[index]);
    modes[index+1] = (time - pwl.get_x()[index  ]) / (pwl.get_x()[index+1] - pwl.get_x()[index]);
  }
}

void PertConstTiltTable::
FeModePref (vector<double > & pref ) const
{
  pref.resize(numMode());
  pwl.value (pwl.get_x(), pref);
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

int PertGaussian::
numMode () const
{
  return 1;
}

void PertGaussian::
FeMode (const double & time,
	vector<double > & modes) const
{
  modes.resize(1);
  modes[0] = 1.;
  if (warmTime > 0.0 && time < warmTime){
    modes[0] *= time / warmTime;
  }
}

void PertGaussian::
FeModePref (vector<double > & pref ) const
{
  pref.resize(1);
  pref[0] = strength;
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


