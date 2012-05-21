#include "Distribution.h"

Distribution_1d::
Distribution_1d (const double & x0_,
		 const double & x1_,
		 const unsigned & nx_,
		 const double & v0_,
		 const double & v1_,
		 const unsigned & nv_)
{
  reinit (x0_, x1_, nx_, v0_, v1_, nv_);
}


void Distribution_1d::
reinit (const double & x0_,
	const double & x1_,
	const unsigned & nx_,
	const double & v0_,
	const double & v1_,
	const unsigned & nv_)
{
  x0 = x0_;
  x1 = x1_;
  v0 = v0_;
  v1 = v1_;
  nx = nx_;
  nv = nv_;

  hx = (x1 - x0) / double (nx);
  hv = (v1 - v0) / double (nv);
  valuepp = 1./(hx*hv);
  
  gridx.resize(nx);
  gridv.resize(nv);

  for (unsigned ii = 0; ii < gridx.size(); ++ii){
    gridx[ii] = x0 + (0.5 + ii) * hx;
  }
  for (unsigned ii = 0; ii < gridv.size(); ++ii){
    gridv[ii] = v0 + (0.5 + ii) * hv;
  }

  values.resize(nx);
  for (unsigned ii = 0; ii < nx; ++ii){
    values[ii].resize(nv);
  }

  clear ();
}


void Distribution_1d::
clear ()
{
  nframe = 0.;
  for (unsigned ii = 0; ii < nx; ++ii){
    for (unsigned jj = 0; jj < nv; ++jj){
      values[ii][jj] = 0.;
    }
  }
}


void Distribution_1d::
average ()
{
  if (nframe == 0.) return;
  for (unsigned ii = 0; ii < nx; ++ii){
    for (unsigned jj = 0; jj < nv; ++jj){
      values[ii][jj] /= double(nframe);
    }
  }  
}

#include <fstream>

void Distribution_1d::
print_xv (FILE * fp) const
{
  for (unsigned ii = 0; ii < nx; ++ii){
    for (unsigned jj = 0; jj < nv; ++jj){
      fprintf (fp, "%f %f %f\n", gridx[ii], gridv[jj], values[ii][jj]);
    }
    fprintf (fp, "\n");
  }
}

void Distribution_1d::
print_x (FILE * fp) const
{
  for (unsigned ii = 0; ii < nx; ++ii){
    double avg = 0.;
    for (unsigned jj = 0; jj < nv; ++jj){
      avg += values[ii][jj] * hv;
    }
    fprintf (fp, "%f %f\n", gridx[ii], avg);
  }
}


