#include "Integrator.h"
#include "RandomGenerator.h"

void EulerMaruyama::
step (Dofs & dofs,
      const double time) const
{
  Dofs pvalue;
  if (pert != NULL){
    (*pert) (dofs, time, pvalue);
  }
  else {
    for (unsigned dd = 0; dd < NUMDOFS; ++dd){
      pvalue.xx[dd] = 0.;
      pvalue.vv[dd] = 0.;
    }
  }
  vector<double > fvalue (NUMDOFS);
  if (force != NULL){
    (*force) (dofs, time, fvalue);
  }
  else {
    for (unsigned dd = 0; dd < NUMDOFS; ++dd){
      fvalue[dd] = 0.;
    }
  }

  // printf ("%f\n", fvalue[0]);
  Dofs oldDofs(dofs);
  for (unsigned dd = 0; dd < NUMDOFS; ++dd){
    double rand = RandomGenerator_MT19937::gaussian();
    dofs.xx[dd] = oldDofs.xx[dd] + dt * oldDofs.vv[dd] + dt * pvalue.xx[dd];
    dofs.vv[dd] = oldDofs.vv[dd] + dt * fvalue[dd] - dt * gamma * oldDofs.vv[dd] + dt * pvalue.vv[dd]
	+ sqrtdt * sigma * rand;
    storedw.vv[dd] = rand;
    if (! (dofs.xx[dd] >= -10 && dofs.xx[dd] <= 10. && dofs.vv[dd]>= -20. && dofs.vv[dd] <= 20.)  ){
      fprintf (stderr, "dd: %d   old: %f %f    new: %f %f\n",
    	       dd,
    	       oldDofs.xx[dd], oldDofs.vv[dd],
    	       dofs.xx[dd], dofs.vv[dd]);
      exit (1);
    }
  }
}


EulerMaruyama::
EulerMaruyama ()
    : dt(0.), sqrtdt(0.), gamma(0.), sigma(0.),
      pert(NULL), force(NULL)
{
}

EulerMaruyama::
EulerMaruyama (const double & gamma,
	       const double & kT,
	       const double & d,
	       const Perturbation * p,
	       const Force * f)
{
  reinit (gamma, kT, d, p, f);
}

void EulerMaruyama::
reinit (const double & g,
	const double & kT,
	const double & d,
	const Perturbation * p,
	const Force * f)
{
  dt = d;
  sqrtdt = sqrt(dt);
  gamma = g;
  sigma = sqrt(2. * kT * gamma);
  pert = p;
  force = f;
  for (unsigned dd = 0; dd < NUMDOFS; ++dd){
    storedw.xx[dd] = 0.;
    // storedw.vv[dd] = 0.;
  }
}


