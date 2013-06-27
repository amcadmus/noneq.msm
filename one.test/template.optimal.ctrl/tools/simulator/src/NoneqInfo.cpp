#include "NoneqInfo.h"

NoneqInfo::
NoneqInfo ()
    : ppert (NULL)
{
}


void NoneqInfo::
reinit (const double & beta_,
	const double & x0_,
	const double & x1_,
	const double & dt_,
	const double & noneqTime_,
	const double & noneqCheckFeq_,
	const Perturbation & pert)
{
  beta = beta_;
  x0 = x0_;
  x1 = x1_;
  dt = dt_;
  noneqTime = noneqTime_;
  noneqCheckFeq = noneqCheckFeq_;
  ppert = &pert;

  numMode = ppert->numMode();  
  numCheck = int((noneqTime + 0.5 * noneqCheckFeq) / noneqCheckFeq) + 1;
  checkNumFeq = int((noneqCheckFeq + 0.5 * dt) / dt);
  
  order0.resize(numCheck);
  order0punish.resize(numCheck);

  for (int tt = 0; tt < numCheck; ++tt){
    order0[tt] = 0.;
    order0punish[tt] = 0.;
  }

  ntraj = 0;
  newTraj ();
  ntraj = 0;
}


void NoneqInfo::
newTraj ()
{
  ntraj ++;

  punish = 0;
  countNoneq = 0;
  countNoneqSeg = 0;
}

int NoneqInfo::
inSet (const Dofs & x) 
{
  if (x.xx[0] > (0.5 * (x0 + x1))){
    return 1;
  }
  else {
    return 0;
  }
}


void NoneqInfo::
depositMainTraj (const Dofs & oldx,
		 const Dofs & newx,
		 const double & sigma,
		 const Dofs & dw)
{
  if (countNoneq == 0){
    order0[0] += (-inSet(oldx));
    // countNoneq ++;
    // return;
  }
    
  double nowTime = countNoneq * dt;
  Dofs pvalueTime;
  ppert->operator () (oldx, nowTime, pvalueTime);
  
  double tmp2 = 0.;
  for (unsigned dd = 0; dd < NUMDOFS; ++dd){
    tmp2 += pvalueTime.vv[dd] * pvalueTime.vv[dd];
  }
  punish += beta * 0.5 * tmp2 * dt;

  countNoneq ++;

  if (countNoneq == (countNoneqSeg + 1) * checkNumFeq){
    countNoneqSeg ++;
    if (countNoneqSeg >= numCheck){
      std::cerr << "# the non eq traj is too long, problematic!" << std::endl;
      exit (1);
    }
    order0[countNoneqSeg] += (- inSet(newx));
    order0punish[countNoneqSeg] += (punish);
  }
}



void NoneqInfo::
average ()
{
  for (int tt = 0; tt < numCheck; ++tt){
    order0[tt]				/= double(ntraj);
    order0punish[tt]			/= double(ntraj);
  }
}

#include <mpi.h>
using namespace MPI;

void NoneqInfo::
collectLast () 
{
  double* tmporder0s = (double *) malloc (sizeof(double) * numCheck);
  double* tmporder0r = (double *) malloc (sizeof(double) * numCheck);
  for (int ii = 0; ii < numCheck; ++ii) {
    tmporder0s[ii] = order0[ii];
  }
  COMM_WORLD.Allreduce (tmporder0s, tmporder0r, numCheck, MPI_DOUBLE, SUM);
  int size = COMM_WORLD.Get_size();
  for (int ii = 0; ii < numCheck; ++ii) {
    order0[ii] = tmporder0r[ii] / double(size);
  }
  
  free (tmporder0r);
  free (tmporder0s);
}

