#include "NoneqResponseInfo.h"

NoneqResponseInfo::
NoneqResponseInfo ()
    : ppert (NULL)
{
}


void NoneqResponseInfo::
reinit (const double & x0_,
	const double & x1_,
	const double & dt_,
	const double & noneqTime_,
	const double & noneqCheckFeq_,
	const Perturbation & pert)
{
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
  order1.resize(numCheck);
  order2.resize(numCheck);
  
  for (int jj = 0; jj < numCheck; ++jj){
    order1[jj].resize (numMode);
    order2[jj].resize (numMode);
    for (int kk = 0; kk < numMode; ++kk){
      order2[jj][kk].resize (numMode);
    }
  }

  for (int tt = 0; tt < numCheck; ++tt){
    order0[tt] = 0.;
    for (int jj = 0; jj < numMode; ++jj){
      order1[tt][jj] = 0.;
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

inline int NoneqResponseInfo::
inSet (const Dofs & x) 
{
  if (x.xx[0] > (0.5 * (x0 + x1))){
    return 1;
  }
  else {
    return 0;
  }
}

inline double NoneqResponseInfo::
trajObservable (const Traj & traj)
{
  if (! traj.full()) return 0.;
  int posi = traj.tailPosi();
  const Dofs & oldx (traj.getValue (posi+1));
  const Dofs & newx (traj.getValue (posi));

  if ( ( inSet (oldx) && !inSet(newx) ) ||
       ( inSet (newx) && !inSet(oldx) ) ){
    return 1.;
  }
  else {
    return 0.;
  }    
}

void NoneqResponseInfo::
depositMainTraj (const Traj & traj,
		 const double & sigma,
		 const Dofs & dw)
{  
  if (countNoneq == 0){
    // order0[0] += (-inSet(oldx));
    // countNoneq ++;
    // return;
  }
  int posi = traj.tailPosi();
  const Dofs & oldx (traj.getValue(posi-1));
    
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
      Hjk[jj][kk] += 1./(sigma*sigma) * modes[jj] * modes[kk] * tmp1 * dt;
    }
  }

  countNoneq ++;

  if (countNoneq == (countNoneqSeg + 1) * checkNumFeq){
    countNoneqSeg ++;
    if (countNoneqSeg >= numCheck){
      std::cerr << "# the non eq traj is too long, problematic!" << std::endl;
      exit (1);
    }
    order0[countNoneqSeg] += trajObservable (traj);
    for (int jj = 0; jj < numMode; ++jj){
      order1[countNoneqSeg][jj] +=  trajObservable(traj) * Gj[jj];
      for (int kk = 0; kk < numMode; ++kk){
      	order2[countNoneqSeg][jj][kk] += trajObservable(traj) * (Gj[jj] * Gj[kk] - Hjk[jj][kk]);
      }
    }
  }
}



void NoneqResponseInfo::
average ()
{
  for (int tt = 0; tt < numCheck; ++tt){
    order0[tt]				/= double(ntraj);
    for (int jj = 0; jj < numMode; ++jj){
      order1[tt][jj]			/= double(ntraj);
      for (int kk = 0; kk < numMode; ++kk){
      	order2[tt][jj][kk]		/= double(ntraj);
      }
    }
  }
}

#include <mpi.h>
using namespace MPI;

void NoneqResponseInfo::
collect () 
{
  double * tmporder0s = (double *) malloc (sizeof(double) * numCheck);
  double * tmporder0r = (double *) malloc (sizeof(double) * numCheck);

  // COMM_WORLD.Allreduce (&(order0.back()), &tmporder0r, 1, MPI_DOUBLE, SUM);
  for (int ii = 0; ii < numCheck; ++ii){
    tmporder0s[ii] = order0[ii];
  }
  COMM_WORLD.Allreduce (tmporder0s, tmporder0r, numCheck, MPI_DOUBLE, SUM);
  int size = COMM_WORLD.Get_size();
  for (int ii = 0; ii < numCheck; ++ii){
    order0[ii] = tmporder0r[ii] / double(size);
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

void NoneqResponseInfo::
calculate (const double & time,
	   const Perturbation & pert1,
	   double & rate,
	   const int order)
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
  if (order >= 0){
    rate = order0[index];
  }
  // order 1
  if (order >= 1){
    for (int jj = 0; jj < numMode; ++jj){
      rate += pref1[jj] * order1[index][jj];
    }
  }
  // order 2
  if (order >= 2){
    for (int jj = 0; jj < numMode; ++jj){
      for (int kk = 0; kk < numMode; ++kk){
  	rate += 0.5 * pref1[jj] * pref1[kk] * order2[index][jj][kk];
      }
    }
  }
}



void NoneqResponseInfo::
save (const string & filename) const
{
  FILE * fp = fopen (filename.c_str(), "w");

  size_t rv;
  rv = fwrite (&x0, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error writing corr file " << endl;
    exit(1);
  }
  rv = fwrite (&x1, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error writing corr file " << endl;
    exit(1);
  }

  rv = fwrite (&dt, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error writing corr file " << endl;
    exit(1);
  }
  rv = fwrite (&noneqTime, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error writing corr file " << endl;
    exit(1);
  }
  rv = fwrite (&noneqCheckFeq, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error writing corr file " << endl;
    exit(1);
  }
  rv = fwrite (&numMode, sizeof(int), 1, fp);
  if (rv != 1){
    cerr << "error writing corr file " << endl;
    exit(1);
  }

  double * buff0 = (double *) malloc (sizeof(double) * numCheck);
  double * buff1 = (double *) malloc (sizeof(double) * numCheck * numMode);
  double * buff2 = (double *) malloc (sizeof(double) * numCheck * numMode * numMode);

  for (int ii = 0; ii < numCheck; ++ii){
    buff0[ii] = order0[ii];
    for (int jj = 0; jj < numMode; ++jj){
      buff1[ii*numMode+jj] = order1[ii][jj];
      for (int kk = 0; kk < numMode; ++kk){
	buff2[ii * numMode * numMode + jj * numMode + kk] = order2[ii][jj][kk];
      }
    }
  }
  rv = fwrite (buff0, sizeof(double), numCheck, fp);
  if (int(rv) != numCheck){
    cerr << "error writing corr file " << endl;
    exit(1);
  }
  rv = fwrite (buff1, sizeof(double), numCheck * numMode, fp);
  if (int(rv) != numCheck * numMode){
    cerr << "error writing corr file " << endl;
    exit(1);
  }
  rv = fwrite (buff2, sizeof(double), numCheck * numMode * numMode, fp);
  if (int(rv) != numCheck * numMode * numMode){
    cerr << "error writing corr file " << endl;
    exit(1);
  }
  
  free (buff0);
  free (buff1);
  free (buff2);
  
  fclose (fp);
}


void NoneqResponseInfo::
load (const string & filename)
{
  FILE * fp = fopen (filename.c_str(), "r");

  size_t rv;
  rv = fread (&x0, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error reading corr file " << endl;
    exit(1);
  }
  rv = fread (&x1, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error reading corr file " << endl;
    exit(1);
  }

  rv = fread (&dt, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error reading corr file " << endl;
    exit(1);
  }
  rv = fread (&noneqTime, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error reading corr file " << endl;
    exit(1);
  }
  rv = fread (&noneqCheckFeq, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error reading corr file " << endl;
    exit(1);
  }
  numCheck = int((noneqTime + 0.5 * noneqCheckFeq) / noneqCheckFeq) + 1;
  checkNumFeq = int((noneqCheckFeq + 0.5 * dt) / dt);

  rv = fread (&numMode, sizeof(int), 1, fp);
  if (rv != 1){
    cerr << "error reading corr file " << endl;
    exit(1);
  }
  if (numMode != ppert->numMode()){
    cerr << "inconsistent perturbation with saved!" << endl;
    exit (1);
  }
  
  order0.resize(numCheck);
  order1.resize(numCheck);
  order2.resize(numCheck);

  for (int jj = 0; jj < numCheck; ++jj){
    order1[jj].resize (numMode);
    order2[jj].resize (numMode);
    for (int kk = 0; kk < numMode; ++kk){
      order2[jj][kk].resize (numMode);
    }
  }  
  
  double * buff0 = (double *) malloc (sizeof(double) * numCheck);
  double * buff1 = (double *) malloc (sizeof(double) * numCheck * numMode);
  double * buff2 = (double *) malloc (sizeof(double) * numCheck * numMode * numMode);

  rv = fread (buff0, sizeof(double), numCheck, fp);
  if (int(rv) != numCheck){
    cerr << "error writing corr file " << endl;
    exit(1);
  }
  rv = fread (buff1, sizeof(double), numCheck * numMode, fp);
  if (int(rv) != numCheck * numMode){
    cerr << "error writing corr file " << endl;
    exit(1);
  }
  rv = fread (buff2, sizeof(double), numCheck * numMode * numMode, fp);
  if (int(rv) != numCheck * numMode * numMode){
    cerr << "error writing corr file " << endl;
    exit(1);
  }

  for (int ii = 0; ii < numCheck; ++ii){
    order0[ii] = buff0[ii];
    for (int jj = 0; jj < numMode; ++jj){
      order1[ii][jj] = buff1[ii*numMode+jj];
      for (int kk = 0; kk < numMode; ++kk){
	order2[ii][jj][kk] = buff2[ii * numMode * numMode + jj * numMode + kk];
      }
    }
  }
  
  free (buff0);
  free (buff1);
  free (buff2);

  newTraj ();
  fclose (fp);
}
