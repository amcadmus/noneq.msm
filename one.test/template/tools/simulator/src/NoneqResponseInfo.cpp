#include "NoneqResponseInfo.h"

NoneqResponseInfo::
NoneqResponseInfo ()
    : ppert (NULL)
{
}


void NoneqResponseInfo::
reinit (const double & x0_,
	const double & x1_,
	const unsigned & nx_,
	const double & v0_,
	const double & v1_,
	const unsigned & nv_,
	const double & dt_,
	const double & noneqTime_,
	const double & noneqCheckFeq_,
	const double & quenchTime_,
	const Perturbation & pert)
{
  x0 = x0_;
  x1 = x1_;
  nx = nx_;
  v0 = v0_;
  v1 = v1_;
  nv = nv_;
  dt = dt_;
  noneqTime = noneqTime_;
  noneqCheckFeq = noneqCheckFeq_;
  quenchTime = quenchTime_;
  ppert = &pert;

  numMode = ppert->numMode();  
  numCheck = int((noneqTime + 0.5 * noneqCheckFeq) / noneqCheckFeq) + 1;
  checkNumFeq = int((noneqCheckFeq + 0.5 * dt) / dt);
  quenchNumStep = int((quenchTime + 0.5 * dt) / dt);
  
  order0.resize(numCheck);
  order1.resize(numCheck);
  order2.resize(numCheck);
  
  quench_order0.resize(numCheck);
  quench_order1_term1.resize(numCheck);
  quench_order1_term2.resize(numCheck);
  quench_order2_term1.resize(numCheck);
  quench_order2_term2.resize(numCheck);
  quench_order2_term3.resize(numCheck);

  for (int jj = 0; jj < numCheck; ++jj){
    order1[jj].resize (numMode);
    order2[jj].resize (numMode);
    quench_order1_term2[jj].resize (numMode);
    quench_order2_term2[jj].resize (numMode);
    quench_order2_term3[jj].resize (numMode);
    for (int kk = 0; kk < numMode; ++kk){
      order2[jj][kk].resize (numMode);
      quench_order2_term2[jj][kk].resize (numMode);
    }
  }

  for (int tt = 0; tt < numCheck; ++tt){
    order0[tt]				.reinit (x0, x1, nx, v0, v1, nv);
    quench_order0[tt]			.reinit (x0, x1, nx, v0, v1, nv);
    quench_order1_term1[tt]		.reinit (x0, x1, nx, v0, v1, nv);
    quench_order2_term1[tt]		.reinit (x0, x1, nx, v0, v1, nv);
    for (int jj = 0; jj < numMode; ++jj){
      order1[tt][jj]			.reinit (x0, x1, nx, v0, v1, nv);
      quench_order1_term2[tt][jj]	.reinit (x0, x1, nx, v0, v1, nv);
      quench_order2_term3[tt][jj]	.reinit (x0, x1, nx, v0, v1, nv);
      for (int kk = 0; kk < numMode; ++kk){
	order2[tt][jj][kk]		.reinit (x0, x1, nx, v0, v1, nv);
	quench_order2_term2[tt][jj][kk] .reinit (x0, x1, nx, v0, v1, nv);
      }
    }
  }


  newTraj ();
}


void NoneqResponseInfo::
newTraj ()
{
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

void NoneqResponseInfo::
newQuenchTraj ()
{
  G0 = 0.;
  H00 = 0.;
  countQuench = 0;  
}


void NoneqResponseInfo::
depositMainTraj (const Dofs & oldx,
		 const Dofs & newx,
		 const double & sigma,
		 const Dofs & dw)
{
  // if (countNoneq - countNoneqSeg * checkNumFeq == 0){
  if (countNoneq == 0){
    order0[0].deposite (oldx);
    countNoneq ++;
    return;
  }
    
  double nowTime = countNoneq * dt;
  Dofs pvalue;
  ppert->operator () (oldx, pvalue);
  std::vector<double > modes;
  ppert->FeMode (nowTime, modes);
  
  double tmp0 = 0.;
  double tmp1 = 0.;
  for (unsigned dd = 0; dd < NUMDOFS; ++dd){
    tmp0 += pvalue.vv[dd] * dw.vv[dd];
    tmp1 += pvalue.vv[dd] * pvalue.vv[dd];
  }
  for (int jj = 0; jj < numMode; ++jj){
    Gj[jj] += 1./sigma * modes[jj] * tmp0;
    for (int kk = 0; kk < numMode; ++kk){
      Hjk[jj][kk] += 1./sigma * modes[jj] * modes[kk] * tmp1 * dt;
    }
  }

  countNoneq ++;
  // std::cout << "noneq count is " << countNoneq
  // 	    << "  countNoneqSeg is " << countNoneqSeg 
  // 	    << std::endl;
  if (countNoneq == (countNoneqSeg + 1) * checkNumFeq){
    countNoneqSeg ++;
    if (countNoneqSeg >= numCheck){
      std::cerr << "# the non eq traj is too long, problematic!" << std::endl;
      exit (1);
    }
    order0[countNoneqSeg].deposite (newx);
    for (int jj = 0; jj < numMode; ++jj){
      order1[countNoneqSeg][jj].deposite (newx, Gj[jj]);
      for (int kk = 0; kk < numMode; ++kk){
	order2[countNoneqSeg][jj][kk].deposite (newx, Hjk[jj][kk]);
      }
    }
  }
}

void NoneqResponseInfo::
depositQuenchTraj (const Dofs & oldx,
		   const Dofs & newx,
		   const double & sigma,
		   const Dofs & dw)
{
  // double nowTime = countQuench * dt;
  
  Dofs pvalue;
  ppert->operator () (oldx, pvalue);

  double tmp0 = 0.;
  double tmp1 = 0.;
  for (unsigned dd = 0; dd < NUMDOFS; ++dd){
    tmp0 += pvalue.vv[dd] * dw.vv[dd];
    tmp1 += pvalue.vv[dd] * pvalue.vv[dd];
  }
  G0 += tmp0 / sigma;
  H00 += tmp1 * dt / sigma;

  countQuench ++;
  // std::cout << "quench count is " << countQuench << std::endl;
  if (countQuench == quenchNumStep){
    quench_order0[countNoneqSeg].deposite (newx);
    if (countNoneqSeg != 0){
      quench_order1_term1[countNoneqSeg].deposite (newx, G0);
      quench_order2_term1[countNoneqSeg].deposite (newx, H00);
      for (int jj = 0; jj < numMode; ++jj){
	quench_order1_term2[countNoneqSeg][jj].deposite (newx, Gj[jj]);
	quench_order2_term3[countNoneqSeg][jj].deposite (newx, G0 * Gj[jj]);
	for (int kk = 0; kk < numMode; ++kk){
	quench_order2_term2[countNoneqSeg][jj][kk].deposite (newx, Gj[jj] * Gj[kk]);
	}
      }
    }
  }
}


void NoneqResponseInfo::
average ()
{
  for (int tt = 0; tt < numCheck; ++tt){
    order0[tt]				.average ();
    quench_order0[tt]			.average ();
    quench_order1_term1[tt]		.average ();
    quench_order2_term1[tt]		.average ();
    for (int jj = 0; jj < numMode; ++jj){
      order1[tt][jj]			.average ();
      quench_order1_term2[tt][jj]	.average ();
      quench_order2_term3[tt][jj]	.average ();
      for (int kk = 0; kk < numMode; ++kk){
	order2[tt][jj][kk]		.average ();
	quench_order2_term2[tt][jj][kk] .average ();
      }
    }
  }
  // for (unsigned tt = 0; tt < numCheck; ++tt){
  //   order0[tt]				.average ();
  //   quench_order0[tt].average ();
  //   quench_order1_term1[tt].average ();
  //   quench_order2_term1[tt].average ();
}


void NoneqResponseInfo::
calculate (const double & time,
	   const Perturbation & pert1,
	   Distribution_1d & dist,
	   Distribution_1d & quench_dist)
{
  int index = int((time + noneqCheckFeq * 0.5) / noneqCheckFeq);
  if (index > numCheck){
    std::cerr << "calculate: time too late!" << std::endl;
    exit (1);
  }
  double nowTime = index * noneqCheckFeq;

  vector<double > pref0, pref1;
  ppert->FeModePref (pref0);
  pert1.FeModePref (pref1);

  double Fe0, Fe1;
  Fe0 = ppert->Fe (nowTime);
  Fe1 = pert1.Fe (nowTime);

  for (int ii = 0; ii < numMode; ++ii){
    pref1[ii] -= pref0[ii];
  }
  Fe1 -= Fe0;

  // order 0
  dist.reinit (x0, x1, nx, v0, v1, nv);
  dist.add (1., order0[index]);
  // order 1
  for (int jj = 0; jj < numMode; ++jj){
    dist.add (pref1[jj], order1[index][jj]);
  }
  // order 2
  for (int jj = 0; jj < numMode; ++jj){
    for (int kk = 0; kk < numMode; ++kk){
      dist.add (-0.5 * pref1[jj] * pref1[kk], order2[index][jj][kk]);
    }
  }

  // order 0
  quench_dist.reinit (x0, x1, nx, v0, v1, nv);
  quench_dist.add (1., quench_order0[index]);
  // order 1
  quench_dist.add (Fe1, quench_order1_term1[index]);
  for (int jj = 0; jj < numMode; ++jj){
    quench_dist.add (pref1[jj], quench_order1_term2[index][jj]);
  }
  // order 2
  quench_dist.add (-0.5 * Fe1 * Fe1, quench_order2_term1[index]);
  for (int jj = 0; jj < numMode; ++jj){
    quench_dist.add (Fe1 * pref1[jj], quench_order2_term3[index][jj]);
    for (int kk = 0; kk < numMode; ++kk){
      quench_dist.add (-0.5 * pref1[jj] * pref1[kk], quench_order2_term2[index][jj][kk]);
    }
  }  
}




