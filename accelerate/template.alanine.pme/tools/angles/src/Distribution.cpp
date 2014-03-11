#include "Distribution.h"

MetastableSet::
MetastableSet (const double & phb,
	       const double & phe,
	       const double & psb,
	       const double & pse)
    : phi_b (phb), phi_e (phe), psi_b (psb), psi_e (pse)
{
}

static bool
metastableSet_inrange (const double & bb,
		       const double & ee,
		       const double & vv)
{
  if (bb < ee){
    return (vv >= bb && vv < ee);
  }
  else {
    return ( (vv >= bb && vv <= 180) || (vv >= -180 && vv < ee) );
  }
}


bool MetastableSet::
inSet (const double & phi,
       const double & psi) const
{
  return ( metastableSet_inrange(phi_b, phi_e, phi) &&
	   metastableSet_inrange(psi_b, psi_e, psi) );
}


Distribution_1d::
Distribution_1d (const double x0_,
		 const double x1_,
		 const unsigned nx_,
		 const double v0_,
		 const double v1_,
		 const unsigned nv_,
		 const unsigned nDataInBlock)
{
  reinit (x0_, x1_, nx_, v0_, v1_, nv_, nDataInBlock);
}


void Distribution_1d::
reinit (const double & x0_,
	const double & x1_,
	const unsigned & nx_,
	const double & v0_,
	const double & v1_,
	const unsigned & nv_,
	const unsigned nDataInBlock)
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
    for (unsigned jj = 0; jj < nv; ++jj){
      values[ii][jj].reinit (nDataInBlock);
    }
  }

  clear ();
}


void Distribution_1d::
clear ()
{
  nframe = 0.;
  for (unsigned ii = 0; ii < nx; ++ii){
    for (unsigned jj = 0; jj < nv; ++jj){
      values[ii][jj].clear ();
    }
  }
}

void Distribution_1d::
deposite (const double & psi,
	  const double & phi)
{
  deposite (psi, phi, 1.);
}

void Distribution_1d::
deposite (const double & psi,
	  const double & phi,
	  const double& scale)
{
  int ix = (psi - x0) / hx;
  int iv = (phi - v0) / hv;
  if ((ix < 0 || ix >= int(nx)) || (iv < 0 || iv >= int(nv))){
    // return;
    if (ix < 0) ix += int(nx);
    if (ix >= int(nx)) ix -= int(nx);
    if (iv < 0) iv += int(nv);
    if (iv >= int(nv)) iv -= int(nv);
  }
  nframe += 1.;
  for (int ii = 0; ii < int(nx); ++ii){
    for (int jj = 0; jj < int(nv); ++jj){
      if (ii == ix && jj == iv){
	values[ix][iv].deposite(valuepp * scale);
      }
      else {
	values[ii][jj].deposite(0.);
      }
    }
  }
}

void Distribution_1d::
average ()
{
  if (nframe == 0.) return;
  for (unsigned ii = 0; ii < nx; ++ii){
    for (unsigned jj = 0; jj < nv; ++jj){
      values[ii][jj].calculate ();
    }
  }  
}

#include <fstream>

void Distribution_1d::
print_xv (const string & filename) const
{
  FILE * fp = fopen (filename.c_str(), "w");
  if (fp == NULL){
    std::cerr << "Distribution_1d::print_xv: cannot open file " << std::endl;
    return;
  }
  print_xv (fp);
  fclose(fp);
}

void Distribution_1d::
print_x (const string & filename) const
{
  FILE * fp = fopen (filename.c_str(), "w");
  if (fp == NULL){
    std::cerr << "Distribution_1d::print_x: cannot open file " << std::endl;
    return;
  }
  print_x (fp);
  fclose(fp);
}

void Distribution_1d::
print_v (const string & filename) const
{
  FILE * fp = fopen (filename.c_str(), "w");
  if (fp == NULL){
    std::cerr << "Distribution_1d::print_x: cannot open file " << std::endl;
    return;
  }
  print_v (fp);
  fclose(fp);
}


void Distribution_1d::
print_xv (FILE * fp) const
{
  for (unsigned ii = 0; ii < nx; ++ii){
    for (unsigned jj = 0; jj < nv; ++jj){
      fprintf (fp, "%f %f %.16e %.16e\n", gridx[ii], gridv[jj], values[ii][jj].getAvg(), values[ii][jj].getAvgError());
    }
    fprintf (fp, "\n");
  }
}

void Distribution_1d::
print_x (FILE * fp) const
{
  for (unsigned ii = 0; ii < nx; ++ii){
    double avg = 0.;
    double avgerr = 0.;
    for (unsigned jj = 0; jj < nv; ++jj){
      avg += values[ii][jj].getAvg() * hv;
      avgerr += values[ii][jj].getAvgError() * values[ii][jj].getAvgError() * hv * hv;
    }
    fprintf (fp, "%f %.16e %.16e\n", gridx[ii], avg, sqrt(avgerr));
  }
}

void Distribution_1d::
print_v (FILE * fp) const
{
  for (unsigned ii = 0; ii < nx; ++ii){
    double avg = 0.;
    double avgerr = 0.;
    for (unsigned jj = 0; jj < nv; ++jj){
      avg += values[jj][ii].getAvg() * hx;
      avgerr += values[jj][ii].getAvgError() * values[jj][ii].getAvgError() * hx * hx;
    }
    fprintf (fp, "%f %.16e %.16e\n", gridv[ii], avg, sqrt(avgerr));
  }
}


// void Distribution_1d::
// substract (const Distribution_1d & d)
// {
//   if (d.nx != nx || d.nv != nv){
//     std::cerr << "unmatch distributions, do nothing" << std::endl;
//     return;
//   }

//   for (unsigned ii = 0; ii < nx; ++ii){
//     for (unsigned jj = 0; jj < nv; ++jj){
//       values[ii][jj] -= d.values[ii][jj];
//     }
//   }
// }

// void Distribution_1d::
// add (const double & scalor,
//      const Distribution_1d & d)
// {
//   if (d.nx != nx || d.nv != nv){
//     std::cerr << "unmatch distributions, do nothing" << std::endl;
//     return;
//   }

//   for (unsigned ii = 0; ii < nx; ++ii){
//     for (unsigned jj = 0; jj < nv; ++jj){
//       values[ii][jj] += scalor * d.values[ii][jj];
//     }
//   }
// }


// void Distribution_1d::
// save (FILE * fp) const
// {
//   size_t rv;
//   rv = fwrite (&x0, sizeof(double), 1, fp);
//   if (rv != 1){
//     cerr << "error writing corr file " << endl;
//     exit(1);
//   }
//   rv = fwrite (&x1, sizeof(double), 1, fp);
//   if (rv != 1){
//     cerr << "error writing corr file " << endl;
//     exit(1);
//   }
//   rv = fwrite (&v0, sizeof(double), 1, fp);
//   if (rv != 1){
//     cerr << "error writing corr file " << endl;
//     exit(1);
//   }
//   rv = fwrite (&v1, sizeof(double), 1, fp);
//   if (rv != 1){
//     cerr << "error writing corr file " << endl;
//     exit(1);
//   }
//   rv = fwrite (&nx, sizeof(unsigned), 1, fp);
//   if (rv != 1){
//     cerr << "error writing corr file " << endl;
//     exit(1);
//   }
//   rv = fwrite (&nv, sizeof(unsigned), 1, fp);
//   if (rv != 1){
//     cerr << "error writing corr file " << endl;
//     exit(1);
//   }

//   double * buff = (double *) malloc (sizeof(double) * nx * nv);
//   for (unsigned ii = 0; ii < nx; ++ii) {
//     for (unsigned jj = 0; jj < nv; ++jj) {
//       buff[ii*nv +jj] = values[ii][jj];
//     }
//   }
//   rv = fwrite (buff, sizeof(double), nx*nv, fp);
//   if (rv != nx*nv){
//     cerr << "error writing corr file " << endl;
//     exit(1);
//   }

//   free (buff);
// }


// bool Distribution_1d::
// load (FILE * fp)
// {
//   size_t rv;
//   rv = fread (&x0, sizeof(double), 1, fp);
//   if (rv != 1){
//     cerr << "error reading corr file or EOF is reached" << endl;
//     return false;
//   }
//   rv = fread (&x1, sizeof(double), 1, fp);
//   if (rv != 1){
//     cerr << "error reading corr file " << endl;
//     exit(1);
//   }
//   rv = fread (&v0, sizeof(double), 1, fp);
//   if (rv != 1){
//     cerr << "error reading corr file " << endl;
//     exit(1);
//   }
//   rv = fread (&v1, sizeof(double), 1, fp);
//   if (rv != 1){
//     cerr << "error reading corr file " << endl;
//     exit(1);
//   }
//   rv = fread (&nx, sizeof(unsigned), 1, fp);
//   if (rv != 1){
//     cerr << "error reading corr file " << endl;
//     exit(1);
//   }
//   rv = fread (&nv, sizeof(unsigned), 1, fp);
//   if (rv != 1){
//     cerr << "error reading corr file " << endl;
//     exit(1);
//   }
//   reinit (x0, x1, nx, v0, v1, nv);

//   double * buff = (double *) malloc (sizeof(double) * nx * nv);
//   rv = fread (buff, sizeof(double), nx*nv, fp);
//   if (rv != nx*nv){
//     cerr << "error reading corr file " << endl;
//     exit(1);
//   }
//   for (unsigned ii = 0; ii < nx; ++ii) {
//     for (unsigned jj = 0; jj < nv; ++jj) {
//       values[ii][jj] = buff[ii*nv +jj];
//     }
//   }
//   free (buff);

//   return true;
// }


