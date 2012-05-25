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
calCorr (const Dofs & initDof,
	 const double & totalTime)
{
  unsigned countStep = 0;
  Dofs dof (initDof);
  vector<double > savedFlux (nFrame, 0.);
  unsigned savePosi = 0;
  unsigned countValidSaved = 0;
  int printCount = 0;
  
  for (double now = 0; now <= totalTime; now += dt){
    if (countStep == calCorrFeq){
      if (++printCount == 100){
	std::cout << "now t: " << now << " / " << totalTime << std::endl;
	printCount = 0;
      }
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
	dists0[ii].deposite (dof, savedFlux[considerPosi]);
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
calIndicator (const vector<vector<double > > & old,
	      const double & idTime,
	      const double & idStep,
	      const double & beta,
	      const Perturbation & pert,
	      vector<vector<vector<double > > > & timeNew)
{
  timeNew.clear();
  if (old.size() != dists0[0].nx || old[0].size() != dists0[0].nv){
    cerr << "inconsistent size of grid" << endl;
    exit (1);
  }
  for (double nowTime = 0; nowTime <= idTime + 0.5 * idStep; nowTime += idStep){
    vector<vector<double > > tmp (old);
    for (double inteTime = 0; inteTime <= nowTime - 0.5 * step; inteTime += step){
      double Fe0 = pert.Fe (nowTime - inteTime);
      double Fe1 = pert.Fe (nowTime - inteTime - step);
      int corr0Idx, corr1Idx;
      if (inteTime <= time + 0.5 * step){
	corr0Idx = (inteTime + 0.5 * step) / step;
      }
      else {
	corr0Idx = -1;
      }
      if (inteTime <= time - 0.5 * step){
	corr1Idx = (inteTime + 1.5 * step) / step;
      }
      else {
	corr1Idx = -1;
      }
	
      for (unsigned ii = 0; ii < dists0[0].nx; ++ii){
	for (unsigned jj = 0; jj < dists0[0].nv; ++jj){
	  double value0 = 0;
	  double value1 = 0;
	  if (corr0Idx >= 0) {
	    value0 = dists0[corr0Idx].values[ii][jj];
	  }
	  if (corr1Idx >= 0) {
	    value1 = dists0[corr1Idx].values[ii][jj];
	  }
	  
	  tmp[ii][jj] -= beta * 0.5 * step * (Fe0 * value0 + Fe1 * value1);
	}
      }
    }
    timeNew.push_back (tmp);
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

void TimeCorrelation::
save (const string & filename) const
{
  FILE * fp = fopen (filename.c_str(), "w");
  if (fp == NULL){
    cerr << "cannot open file " << filename << std::endl;
    exit(1) ;
  }

  fwrite (&(dists0[0].x0), sizeof(double), 1, fp);
  fwrite (&(dists0[0].x1), sizeof(double), 1, fp);
  fwrite (&(dists0[0].v0), sizeof(double), 1, fp);
  fwrite (&(dists0[0].v1), sizeof(double), 1, fp);
  fwrite (&(dists0[0].nx), sizeof(unsigned), 1, fp);
  fwrite (&(dists0[0].nv), sizeof(unsigned), 1, fp);

  fwrite (&step, sizeof(double), 1, fp);
  fwrite (&time, sizeof(double), 1, fp);
  fwrite (&nFrame, sizeof(unsigned), 1, fp);
  fwrite (&nStep,  sizeof(unsigned), 1, fp);
  // cout << step << endl;
  // cout << time << endl;
  // cout << nFrame << endl;
  // cout << nStep << endl;

  double * buff = (double *) malloc (sizeof(double) * dists0[0].nx * dists0[0].nv);

  for (unsigned kk = 0; kk < nFrame; ++kk){
    for (unsigned ii = 0; ii < dists0[0].nx; ++ii) {
      for (unsigned jj = 0; jj < dists0[0].nv; ++jj) {
	buff[ii*dists0[0].nv +jj] = dists0[kk].values[ii][jj];
      }
    }
    fwrite (buff, sizeof(double), dists0[0].nx*dists0[0].nv, fp);
  }

  free (buff);
  fclose (fp);
}

void TimeCorrelation::
load (const string & filename) 
{
  FILE * fp = fopen (filename.c_str(), "r");
  if (fp == NULL){
    cerr << "cannot open file " << filename << std::endl;
    exit(1) ;
  }

  double x0, x1, v0, v1;
  unsigned nx, nv;
  size_t rv;
  rv = fread (&x0, sizeof(double), 1, fp);
  rv = fread (&x1, sizeof(double), 1, fp);
  rv = fread (&v0, sizeof(double), 1, fp);
  rv = fread (&v1, sizeof(double), 1, fp);
  rv = fread (&nx, sizeof(unsigned), 1, fp);
  rv = fread (&nv, sizeof(unsigned), 1, fp);

  rv = fread (&step, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error reading corr file " << endl;
    exit(1);
  }
  rv = fread (&time, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error reading corr file " << endl;
    exit(1);
  }
  rv = fread (&nFrame, sizeof(unsigned), 1, fp);
  if (rv != 1){
    cerr << "error reading corr file " << endl;
    exit(1);
  }
  rv = fread (&nStep,  sizeof(unsigned), 1, fp);
  if (rv != 1){
    cerr << "error reading corr file " << endl;
    exit(1);
  }

  dists0.resize (nFrame);
  dists1.resize (nFrame);
  
  for (unsigned ii = 0; ii < nFrame; ++ii){
    dists0[ii].reinit (x0, x1, nx, v0, v1, nv);
    dists1[ii].reinit (x0, x1, nx, v0, v1, nv);
  }

  double * buff = (double *) malloc (sizeof(double) * nx * nv);

  for (unsigned kk = 0; kk < nFrame; ++kk){
    rv = fread (buff, sizeof(double), nx*nv, fp);
    if (rv != nx*nv){
      cerr << "error reading corr file " << endl;
      exit(1);
    }
    for (unsigned ii = 0; ii < nx; ++ii) {
      for (unsigned jj = 0; jj < nv; ++jj) {
	dists0[kk].values[ii][jj] = buff[ii*nv +jj];
      }
    }
  }

  free (buff);
  fclose (fp);
}

