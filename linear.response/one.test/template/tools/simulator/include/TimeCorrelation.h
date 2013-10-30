#ifndef __TimeCorrelation_h_NONEQMSM_wanghan__
#define __TimeCorrelation_h_NONEQMSM_wanghan__

#include "Defines.h"
#include "Force.h"
#include "Integrator.h"
#include "Distribution.h"

class TimeCorrelation 
{
  double step;
  double time;
  double dt;
  unsigned nFrame;
  unsigned nStep;
  unsigned calCorrFeq;
  std::vector<Distribution_1d > dists0;
  std::vector<Distribution_1d > dists1;
  const Integrator * equiInte;
  double quenchAtTime;
  unsigned quenchNumStep;
  const Integrator * quenchInte;
  const DissipativeFlux * flux;
public:
  void reinit (const double & step,
	       const double & time,
	       const double & x0,
	       const double & x1,
	       const unsigned & nx,
	       const double & v0,
	       const double & v1,
	       const unsigned & nv,
	       const Integrator * equiInte,
	       const Integrator * quenchInte,
	       const unsigned & quenchNumStep,
	       const double & quenchAtTime,
	       const DissipativeFlux * flux);
  void calCorr (const Dofs & initDof,
		const double & totalTime);
  void calIndicator (const vector<vector<double > > & old,
		     const double & idTime,
		     const double & idStep,
		     const double & beta,
		     const Perturbation & pert,
		     vector<vector<vector<double > > > & timeNew);
  void print () const ;
  void save (const string & filename) const;
  void load (const string & filename);
}
    ;


#endif
