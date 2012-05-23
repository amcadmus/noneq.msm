#include "TimeCorrelation.h"

void TimeCorrelation::
reinit (const double & step_,
	const double & time_,
	const double & x0,
	const double & x1,
	const unsigned & nx,
	const double & v0,
	const double & v1,
	const unsigned & nv,
	const Integrator * equiInte_,
	const unsigned quenchNumStep_,
	const Integrator * quenchInte_,
	const DissipativeFlux * flux_)
{
  equiInte = equiInte_;
  quenchInte = quenchInte_;
  flux = flux_;
  step = step_;
  time = time_;
  dt = equiInte->getDt();
  if (step / dt < 1){
    std::cerr << "step should be larger than dt" << std::endl;
    exit (1);
  }
  nFrame = (time + dt) / step + 1;
  nStep = nFrame - 1;
  quenchNumStep = quenchNumStep_;
  calCorrFeq = (step + .5 * dt) / dt;

  dists0.resize (nFrame);
  dists1.resize (nFrame);
  
  for (unsigned ii = 0; ii < nFrame; ++ii){
    dists0[ii].reinit (x0, x1, nx, v0, v1, nv);
    dists1[ii].reinit (x0, x1, nx, v0, v1, nv);
  }
}

static inline unsigned
decreasePosi (const vector<double > & ring,
	      const unsigned &posi)
{
  int p (posi);
  if ((--p) < 0){
    p += int (ring.size());
  }
  return unsigned(p);
}

static inline unsigned
increasePosi (const vector<double > & ring,
	      const unsigned &posi)
{
  int p (posi);
  if ((++p) >= int(ring.size())){
    p -= int (ring.size());
  }
  return unsigned(p);
}

void TimeCorrelation::
calculate (const Dofs & initDof,
	   const double & totalTime)
{
  unsigned countStep = 0;
  Dofs dof (initDof);
  vector<double > savedFlux (nFrame, 0.);
  unsigned savePosi = 0;
  unsigned countValidSaved = 0;
  
  for (double now = 0; now <= totalTime; now += dt){
    if (countStep == calCorrFeq){
      std::cout << "now t: " << now << " / " << totalTime << std::endl;
      countStep = 0;
      Dofs quenchDof (dof);
      for (unsigned tt = 0; tt < quenchNumStep; ++tt){
	quenchInte->step (quenchDof, 0.);
      }
      savedFlux[savePosi] = (*flux) (dof);
      if (countValidSaved < nFrame) countValidSaved ++;
      unsigned considerPosi = savePosi;
      savePosi = increasePosi (savedFlux, savePosi);
      for (unsigned ii = 0; ii < countValidSaved; ++ii){
	dists0[ii].deposite (dof, -savedFlux[considerPosi]);
	dists1[ii].deposite (quenchDof, savedFlux[considerPosi]);
	considerPosi = decreasePosi (savedFlux, considerPosi);
      }
    }
    equiInte->step (dof, 0.);
    countStep ++;
  }

  for (unsigned kk = 0; kk < nFrame; ++kk){
    dists0[kk].average();
    dists1[kk].average();
    for (unsigned ii = 0; ii < dists0[kk].nx; ++ii){
      for (unsigned jj = 0; jj < dists0[kk].nv; ++jj){
	dists0[kk].values[ii][jj] -= dists1[kk].values[ii][jj];
      }
    }
  }
}

void TimeCorrelation::
print () const
{
  for (unsigned kk = 0; kk < nFrame; ++kk){
    double printTime = step * kk;
    int timeI = int(printTime);
    int timeF = int(100 * (printTime - timeI));
    char name[2048];
    sprintf (name, "corr.x.%05d.%02d.out", timeI, timeF);
    dists0[kk].print_x (string(name));
    sprintf (name, "corr.vx.%05d.%02d.out", timeI, timeF);
    dists0[kk].print_xv (string(name));
  }
}

