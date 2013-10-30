#include "NoneqResponseInfo.h"

NoneqResponseInfo::
NoneqResponseInfo ()
    : ppert (NULL)
{
}


void NoneqResponseInfo::
reinit (const double & beta_,
	const double & x0_,
	const double & x1_,
	const unsigned & nx_,
	const double & v0_,
	const double & v1_,
	const unsigned & nv_,
	const double & dt_,
	const double & noneqTime_,
	const double & noneqCheckFeq_,
	const Perturbation & pert)
{
  beta = beta_;
  x0 = x0_;
  x1 = x1_;
  nx = nx_;
  v0 = v0_;
  v1 = v1_;
  nv = nv_;
  dt = dt_;
  noneqTime = noneqTime_;
  noneqCheckFeq = noneqCheckFeq_;
  ppert = &pert;

  numMode = ppert->numMode();  
  numCheck = int((noneqTime + 0.5 * noneqCheckFeq) / noneqCheckFeq) + 1;
  checkNumFeq = int((noneqCheckFeq + 0.5 * dt) / dt);
  
  order0.resize(numCheck);
  order0error.resize(numCheck);
  order0punish.resize(numCheck);
  order1.resize(numCheck);
  order1error.resize(numCheck);
  order2.resize(numCheck);
  
  for (int jj = 0; jj < numCheck; ++jj){
    order1[jj].resize (numMode);
    order1error[jj].resize (numMode);
    order2[jj].resize (numMode);
    for (int kk = 0; kk < numMode; ++kk){
      order2[jj][kk].resize (numMode);
    }
  }

  for (int tt = 0; tt < numCheck; ++tt){
    order0[tt] = 0.;
    order0error[tt] = 0.;
    order0punish[tt] = 0.;
    for (int jj = 0; jj < numMode; ++jj){
      order1[tt][jj] = 0.;
      order1error[tt][jj] = 0.;
      for (int kk = 0; kk < numMode; ++kk){
	order2[tt][jj][kk] = 0.;
      }
    }
  }

  ntraj = 0;
  newTraj ();
  ntraj = 0;
}


void NoneqResponseInfo::
newTraj ()
{
  ntraj ++;

  punish = 0;
  Gj.resize (numMode);
  Hjk.resize (numMode);
  for (int jj = 0; jj < numMode; ++jj){
    Gj[jj] = 0.;
    Hjk[jj].resize (numMode);
    for (int kk = 0; kk < numMode; ++kk){
      Hjk[jj][kk] = 0.;
    }
  }
  countNoneq = 0;
  countNoneqSeg = 0;
}

int NoneqResponseInfo::
inSet (const Dofs & x) 
{
  if (x.xx[0] > (0.5 * (x0 + x1))){
    return 1;
  }
  else {
    return 0;
  }
}


void NoneqResponseInfo::
depositMainTraj (const Dofs & oldx,
		 const Dofs & newx,
		 const double & sigma,
		 const Dofs & dw)
{
  if (countNoneq == 0){
    order0[0] += (-inSet(oldx));
    order0error[0] += (-inSet(oldx)) * (-inSet(oldx));
    // countNoneq ++;
    // return;
  }
    
  double nowTime = countNoneq * dt;
  Dofs pvalue;
  ppert->operator () (oldx, pvalue);
  Dofs pvalueTime;
  ppert->operator () (oldx, nowTime, pvalueTime);
  std::vector<double > modes;
  ppert->FeMode (nowTime, modes);
  
  double tmp0 = 0.;
  double tmp1 = 0.;
  double tmp2 = 0.;
  for (unsigned dd = 0; dd < NUMDOFS; ++dd){
    tmp0 += pvalue.vv[dd] * dw.vv[dd];
    tmp1 += pvalue.vv[dd] * pvalue.vv[dd];
    tmp2 += pvalueTime.vv[dd] * pvalueTime.vv[dd];
  }
  for (int jj = 0; jj < numMode; ++jj){
    Gj[jj] += 1./sigma * modes[jj] * tmp0;
    // Gj[jj] += - 1./sigma * modes[jj] * tmp0;
    // for (int kk = 0; kk < numMode; ++kk){
    //   Hjk[jj][kk] += 1./(sigma*sigma) * modes[jj] * modes[kk] * tmp1 * dt;
    // }
  }
  punish += beta * 0.5 * tmp2 * dt;

  countNoneq ++;

  if (countNoneq == (countNoneqSeg + 1) * checkNumFeq){
    countNoneqSeg ++;
    if (countNoneqSeg >= numCheck){
      std::cerr << "# the non eq traj is too long, problematic!" << std::endl;
      exit (1);
    }
    // order0[countNoneqSeg] += (- inSet(newx) + punish);
    order0[countNoneqSeg] += (- inSet(newx));
    order0error[countNoneqSeg] += (- inSet(newx)) * (- inSet(newx));
    order0punish[countNoneqSeg] += (punish);
    for (int jj = 0; jj < numMode; ++jj){
      // order1[countNoneqSeg][jj] +=  (- inSet(newx) + punish) * Gj[jj];
      order1[countNoneqSeg][jj] +=  (- inSet(newx)) * Gj[jj];
      order1error[countNoneqSeg][jj] +=  (- inSet(newx)) * Gj[jj] * (- inSet(newx)) * Gj[jj];
      // for (int kk = 0; kk < numMode; ++kk){
      // 	order2[countNoneqSeg][jj][kk] = (- inSet(newx) + punish) * (Gj[jj] * Gj[kk] - Hjk[jj][kk]);
      // }
    }
  }
}



void NoneqResponseInfo::
average ()
{
  for (int tt = 0; tt < numCheck; ++tt){
    order0[tt]				/= double(ntraj);
    order0error[tt]			/= double(ntraj);
    order0punish[tt]			/= double(ntraj);
    for (int jj = 0; jj < numMode; ++jj){
      order1[tt][jj]			/= double(ntraj);
      order1error[tt][jj]		/= double(ntraj);
      // for (int kk = 0; kk < numMode; ++kk){
      // 	order2[tt][jj][kk]		/= double(ntraj);
      // }
    }
  }
}

#include <mpi.h>
using namespace MPI;

void NoneqResponseInfo::
collectLast () 
{
  double tmporder0r;
  double * tmporder1s = (double *) malloc (sizeof(double) * numMode);
  double * tmporder1r = (double *) malloc (sizeof(double) * numMode);

  COMM_WORLD.Allreduce (&(order0.back()), &tmporder0r, 1, MPI_DOUBLE, SUM);
  for (int ii = 0; ii < numMode; ++ii){
    tmporder1s[ii] = order1.back()[ii];
  }
  COMM_WORLD.Allreduce (tmporder1s, tmporder1r, numMode, MPI_DOUBLE, SUM);
  int size = COMM_WORLD.Get_size();
  for (int ii = 0; ii < numMode; ++ii){
    order1.back()[ii] = tmporder1r[ii] / double(size);
  }
  
  free (tmporder1s);
  free (tmporder1r);
}


void NoneqResponseInfo::
collect () 
{
  double * tmporder0s = (double *) malloc (sizeof(double) * numCheck);
  double * tmporder0r = (double *) malloc (sizeof(double) * numCheck);

  for (int ii = 0; ii < numCheck; ++ii){
    tmporder0s[ii] = order0[ii];
  }
  COMM_WORLD.Allreduce (tmporder0s, tmporder0r, numCheck, MPI_DOUBLE, SUM);
  int size = COMM_WORLD.Get_size();
  for (int ii = 0; ii < numCheck; ++ii){
    order0[ii] = tmporder0r[ii] / double(size);
  }
  for (int ii = 0; ii < numCheck; ++ii){
    tmporder0s[ii] = order0error[ii];
  }
  COMM_WORLD.Allreduce (tmporder0s, tmporder0r, numCheck, MPI_DOUBLE, SUM);
  size = COMM_WORLD.Get_size();
  for (int ii = 0; ii < numCheck; ++ii){
    order0error[ii] = tmporder0r[ii] / double(size);
  }
  
  free (tmporder0s);
  free (tmporder0r);

  double * tmporder1s = (double *) malloc (sizeof(double) * numCheck * numMode);
  double * tmporder1r = (double *) malloc (sizeof(double) * numCheck * numMode);

  for (int ii = 0; ii < numCheck; ++ii){
    for (int jj = 0; jj < numMode; ++jj){
      tmporder1s[ii*numMode + jj] = order1[ii][jj];
    }
  }
  COMM_WORLD.Allreduce (tmporder1s, tmporder1r, numCheck * numMode, MPI_DOUBLE, SUM);
  for (int ii = 0; ii < numCheck; ++ii){
    for (int jj = 0; jj < numMode; ++jj){
      order1[ii][jj] = tmporder1r[ii*numMode + jj] / double(size);
    }
  }
  for (int ii = 0; ii < numCheck; ++ii){
    for (int jj = 0; jj < numMode; ++jj){
      tmporder1s[ii*numMode + jj] = order1error[ii][jj];
    }
  }
  COMM_WORLD.Allreduce (tmporder1s, tmporder1r, numCheck * numMode, MPI_DOUBLE, SUM);
  for (int ii = 0; ii < numCheck; ++ii){
    for (int jj = 0; jj < numMode; ++jj){
      order1error[ii][jj] = tmporder1r[ii*numMode + jj] / double(size);
    }
  }
  
  free (tmporder1s);
  free (tmporder1r);


  double * tmporder2s = (double *) malloc (sizeof(double) * numCheck * numMode * numMode);
  double * tmporder2r = (double *) malloc (sizeof(double) * numCheck * numMode * numMode);

  for (int ii = 0; ii < numCheck; ++ii){
    for (int jj = 0; jj < numMode; ++jj){
      for (int kk = 0; kk < numMode; ++kk){ 
	tmporder2s[ii*numMode*numMode + jj*numMode + kk] = order2[ii][jj][kk];
      }
    }
  }
  COMM_WORLD.Allreduce (tmporder2s, tmporder2r, numCheck * numMode * numMode, MPI_DOUBLE, SUM);
  for (int ii = 0; ii < numCheck; ++ii){
    for (int jj = 0; jj < numMode; ++jj){
      for (int kk = 0; kk < numMode; ++kk){ 
	order2[ii][jj][kk] = tmporder2r[ii*numMode*numMode + jj*numMode + kk] / double(size);
      }
    }
  }
  
  free (tmporder2s);
  free (tmporder2r);

}


// void NoneqResponseInfo::
// calculate (const double & time,
// 	   const Perturbation & pert1,
// 	   Distribution_1d & dist,
// 	   Distribution_1d & quench_dist,
// 	   const int order)
// {
//   int index = int((time + noneqCheckFeq * 0.5) / noneqCheckFeq);
//   if (index > numCheck){
//     std::cerr << "calculate: time too late!" << std::endl;
//     exit (1);
//   }
//   double nowTime = index * noneqCheckFeq;

//   vector<double > pref0, pref1;
//   ppert->FeModePref (pref0);
//   pert1.FeModePref (pref1);

//   double Fe0, Fe1;
//   Fe0 = ppert->Fe (nowTime);
//   Fe1 = pert1.Fe (nowTime);

//   for (int ii = 0; ii < numMode; ++ii){
//     pref1[ii] -= pref0[ii];
//   }
//   Fe1 -= Fe0;

//   dist.reinit (x0, x1, nx, v0, v1, nv);
//   // order 0
//   if (order >= 0){
//     dist.add (1., order0[index]);
//   }
//   // order 1
//   if (order >= 1){
//     for (int jj = 0; jj < numMode; ++jj){
//       dist.add (pref1[jj], order1[index][jj]);
//     }
//   }
//   // order 2
//   if (order >= 2){
//     for (int jj = 0; jj < numMode; ++jj){
//       for (int kk = 0; kk < numMode; ++kk){
// 	dist.add (0.5 * pref1[jj] * pref1[kk], order2[index][jj][kk]);
//       }
//     }
//   }

//   quench_dist.reinit (x0, x1, nx, v0, v1, nv);
//   // order 0
//   if (order >= 0){
//     quench_dist.add (1., quench_order0[index]);
//   }
//   // order 1
//   if (order >= 1){
//     quench_dist.add (Fe1, quench_order1_term1[index]);
//     for (int jj = 0; jj < numMode; ++jj){
//       quench_dist.add (pref1[jj], quench_order1_term2[index][jj]);
//     }
//   }
//   // order 2
//   if (order >= 2){
//     quench_dist.add (0.5 * Fe1 * Fe1, quench_order2_term1[index]);
//     for (int jj = 0; jj < numMode; ++jj){
//       quench_dist.add (Fe1 * pref1[jj], quench_order2_term3[index][jj]);
//       for (int kk = 0; kk < numMode; ++kk){
// 	quench_dist.add (0.5 * pref1[jj] * pref1[kk], quench_order2_term2[index][jj][kk]);
//       }
//     }  
//   }
// }




