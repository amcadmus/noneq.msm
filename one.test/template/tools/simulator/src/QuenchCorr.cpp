#include "QuenchCorr.h"

QuenchCorr::
QuenchCorr (const double & x0,
	    const double & x1,
	    const unsigned & nx,
	    const double & v0,
	    const double & v1,
	    const unsigned & nv,
	    const EulerMaruyama & quench_inte0,
	    const Perturbation & pert_)
{
  corr.reinit (x0, x1, nx, v0, v1, nv);
  sigma = quench_inte0.getSigma();
  dt = quench_inte0.getDt();
  pert = &pert_;
}

void QuenchCorr::
deposit (const Dofs & endpoint,
	 const std::vector<Dofs > & dofs,
	 const std::vector<Dofs > & dw)
{
  double Gx = calG(dofs, dw);
  corr.deposite (endpoint, Gx);
}

double QuenchCorr::
calG (const std::vector<Dofs > & dofs,
      const std::vector<Dofs > & dw)
{
  double sum = 0.;
  for (unsigned ii = 0; ii < dofs.size(); ++ii){
    Dofs pvalue;
    pert->operator()(dofs[ii], pvalue);
    double dot = 0.;
    for (unsigned dd = 0; dd < NUMDOFS; ++dd){
      dot += pvalue.vv[dd] * dw[ii].vv[dd];
    }
    sum += 1./sigma * dot;
  }
  return sum;
}

void QuenchCorr::
average ()
{
  corr.average();
}
