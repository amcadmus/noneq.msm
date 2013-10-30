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
		 const unsigned nv_)
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
  backup_values.resize(nx);
  for (unsigned ii = 0; ii < nx; ++ii){
    backup_values[ii].resize(nv);
  }

  backup_number = 1e6;
  backup_unbacked_count = 0.;
  clear ();
}


void Distribution_1d::
clear ()
{
  nframe = 0.;
  for (unsigned ii = 0; ii < nx; ++ii){
    for (unsigned jj = 0; jj < nv; ++jj){
      values[ii][jj] = 0.;
      backup_values[ii][jj] = 0.;
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
    return;
  }
  else {
    nframe += 1.;
    backup_unbacked_count += 1.;
    values[ix][iv] += valuepp * scale;
    if (backup_unbacked_count >= backup_number){
      // std::cout << "backuped" << std::endl;
      /* fprintf (stderr, "backup: %.16e %.16e\n", backup_unbacked_count, backup_number); */
      for (unsigned ii = 0; ii < nx; ++ii){
	for (unsigned jj = 0; jj < nv; ++jj){
	  /* if (ii == nx/2 && jj == nv/2){ */
	  /*   fprintf (stderr,"average: %.16e %.16e    %.16e\n" , backup_values[ii][jj], values[ii][jj], nframe); */
	  /* } */
	  backup_values[ii][jj] += values[ii][jj];
	  values[ii][jj] = 0.;
	}
      }
      backup_unbacked_count = 0.;
    }
  }
}

void Distribution_1d::
average ()
{
  if (nframe == 0.) return;
  for (unsigned ii = 0; ii < nx; ++ii){
    for (unsigned jj = 0; jj < nv; ++jj){
      // if (ii == nx/2 && jj == nv/2){
      // 	fprintf (stderr,"average: %.16e %.16e     %.16e\n" , backup_values[ii][jj], values[ii][jj], nframe);
      // }
      values[ii][jj] += backup_values[ii][jj];
      values[ii][jj] /= double(nframe);
    }
  }  
}

#include <fstream>

void Distribution_1d::
print_xv (const string & filename) const
{
  FILE * fp = fopen (filename.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << std::endl;
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
    std::cerr << "cannot open file " << std::endl;
    return;
  }
  print_x (fp);
  fclose(fp);
}


void Distribution_1d::
print_xv (FILE * fp) const
{
  for (unsigned ii = 0; ii < nx; ++ii){
    for (unsigned jj = 0; jj < nv; ++jj){
      fprintf (fp, "%f %f %.16e\n", gridx[ii], gridv[jj], values[ii][jj]);
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
    fprintf (fp, "%f %.16e\n", gridx[ii], avg);
  }
}


void Distribution_1d::
substract (const Distribution_1d & d)
{
  if (d.nx != nx || d.nv != nv){
    std::cerr << "unmatch distributions, do nothing" << std::endl;
    return;
  }

  for (unsigned ii = 0; ii < nx; ++ii){
    for (unsigned jj = 0; jj < nv; ++jj){
      values[ii][jj] -= d.values[ii][jj];
    }
  }
}

void Distribution_1d::
add (const double & scalor,
     const Distribution_1d & d)
{
  if (d.nx != nx || d.nv != nv){
    std::cerr << "unmatch distributions, do nothing" << std::endl;
    return;
  }

  for (unsigned ii = 0; ii < nx; ++ii){
    for (unsigned jj = 0; jj < nv; ++jj){
      values[ii][jj] += scalor * d.values[ii][jj];
    }
  }
}


void Distribution_1d::
save (FILE * fp) const
{
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
  rv = fwrite (&v0, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error writing corr file " << endl;
    exit(1);
  }
  rv = fwrite (&v1, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error writing corr file " << endl;
    exit(1);
  }
  rv = fwrite (&nx, sizeof(unsigned), 1, fp);
  if (rv != 1){
    cerr << "error writing corr file " << endl;
    exit(1);
  }
  rv = fwrite (&nv, sizeof(unsigned), 1, fp);
  if (rv != 1){
    cerr << "error writing corr file " << endl;
    exit(1);
  }

  double * buff = (double *) malloc (sizeof(double) * nx * nv);
  for (unsigned ii = 0; ii < nx; ++ii) {
    for (unsigned jj = 0; jj < nv; ++jj) {
      buff[ii*nv +jj] = values[ii][jj];
    }
  }
  rv = fwrite (buff, sizeof(double), nx*nv, fp);
  if (rv != nx*nv){
    cerr << "error writing corr file " << endl;
    exit(1);
  }

  free (buff);
}


bool Distribution_1d::
load (FILE * fp)
{
  size_t rv;
  rv = fread (&x0, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error reading corr file or EOF is reached" << endl;
    return false;
  }
  rv = fread (&x1, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error reading corr file " << endl;
    exit(1);
  }
  rv = fread (&v0, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error reading corr file " << endl;
    exit(1);
  }
  rv = fread (&v1, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error reading corr file " << endl;
    exit(1);
  }
  rv = fread (&nx, sizeof(unsigned), 1, fp);
  if (rv != 1){
    cerr << "error reading corr file " << endl;
    exit(1);
  }
  rv = fread (&nv, sizeof(unsigned), 1, fp);
  if (rv != 1){
    cerr << "error reading corr file " << endl;
    exit(1);
  }
  reinit (x0, x1, nx, v0, v1, nv);

  double * buff = (double *) malloc (sizeof(double) * nx * nv);
  rv = fread (buff, sizeof(double), nx*nv, fp);
  if (rv != nx*nv){
    cerr << "error reading corr file " << endl;
    exit(1);
  }
  for (unsigned ii = 0; ii < nx; ++ii) {
    for (unsigned jj = 0; jj < nv; ++jj) {
      values[ii][jj] = buff[ii*nv +jj];
    }
  }
  free (buff);

  return true;
}


