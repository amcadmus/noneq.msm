#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
#include <cmath>
#include <boost/program_options.hpp>
#define MaxLineLength 10240

#include "Interpolation.h"

namespace po = boost::program_options;
using namespace std;

// class SwitchFunction
// {
//   double start, end, inter;
// public:
//   SwitchFunction (const double & start_,
// 		  const double & end_);
//   double u  (const double & rr) const;
//   double up (const double & rr) const;
// }
//     ;

class SwitchFunction
{
  Poly myp;
  Poly mypp;
  double start, end, inter;
public:
  SwitchFunction (const double & start_,
		  const double & end_);
  double u  (const double & rr) const;
  double up (const double & rr) const;
}
    ;

class C1
{
public:
  double u  (const double & rr) const;
  double up (const double & rr) const;
}
    ;

class C1_Smooth
{
  C1 c1;
  SwitchFunction sw;
public:
  C1_Smooth (const double & smooth_start,
	     const double & smooth_end);
  double u  (const double & rr) const;
  double up (const double & rr) const;
}
    ;

class C6
{
public:
  double u  (const double & rr) const;
  double up (const double & rr) const;
}
    ;

class C6_Smooth
{
  C6 c6;
  SwitchFunction sw;
public:
  C6_Smooth (const double & smooth_start,
	     const double & smooth_end);
  double u  (const double & rr) const;
  double up (const double & rr) const;
}
    ;

class C12
{
public:
  double u  (const double & rr) const;
  double up (const double & rr) const;
}
    ;

class C12_Smooth
{
  C12 c12;
  SwitchFunction sw;
public:
  C12_Smooth (const double & smooth_start,
	      const double & smooth_end);
  double u  (const double & rr) const;
  double up (const double & rr) const;
}
    ;

SwitchFunction::
SwitchFunction (const double & start_,
		const double & end_)
    : start(start_), end(end_)
{
  inter = end - start;
  Interpolation::piece6OrderInterpol (start, end, 0, 1, 0, 0, 0, 0, myp);
  mypp = myp;
  mypp.derivative ();
}

double SwitchFunction::
u  (const double & rr) const 
{
  if (rr < start) {
    return 0.;
  }
  else if (rr > end){
    return 1.;
  }
  else {
    return myp.value (rr);
  }
}

double SwitchFunction::
up (const double & rr) const 
{
  if (rr < start) {
    return 0.;
  }
  else if (rr > end){
    return 0.;
  }
  else {
    return mypp.value (rr);
  }
}

// SwitchFunction::
// SwitchFunction (const double & start_,
// 		const double & end_)
//     : start(start_), end(end_)
// {
//   inter = end - start;
// }

// double SwitchFunction::
// u  (const double & rr) const 
// {
//   if (rr < start) {
//     return 0.;
//   }
//   else if (rr > end){
//     return 1.;
//   }
//   else {
//     return 0.5 + 0.5 * cos ( (rr - start) / inter * M_PI + M_PI );
//   }
// }

// double SwitchFunction::
// up (const double & rr) const 
// {
//   if (rr < start) {
//     return 0.;
//   }
//   else if (rr > end){
//     return 0.;
//   }
//   else {
//     return - 0.5 * sin ( (rr - start) / inter * M_PI + M_PI ) * M_PI / inter;
//   }
// }



double C1::
u  (const double & rr) const 
{
  return 1./(rr);
}

double C1::
up (const double & rr) const 
{
  return -1./(rr * rr);
}

C1_Smooth::
C1_Smooth (const double & smooth_start,
	   const double & smooth_end)
    : sw(smooth_start, smooth_end)
{
}

double C1_Smooth::
u (const double & rr) const
{
  return (1. - sw.u(rr)) * c1.u(rr);
}

double C1_Smooth::
up (const double & rr) const
{
  return (1. - sw.u(rr)) * c1.up(rr) + ( - sw.up(rr)) * c1.u(rr);
}


double C6::
u  (const double & rr) const 
{
  return 1./(rr * rr * rr * rr * rr * rr);
}

double C6::
up (const double & rr) const 
{
  return -6./(rr * rr * rr * rr * rr * rr * rr);
}

C6_Smooth::
C6_Smooth (const double & smooth_start,
	   const double & smooth_end)
    : sw(smooth_start, smooth_end)
{
}

double C6_Smooth::
u (const double & rr) const
{
  return (1. - sw.u(rr)) * c6.u(rr);
}

double C6_Smooth::
up (const double & rr) const
{
  return (1. - sw.u(rr)) * c6.up(rr) + ( - sw.up(rr)) * c6.u(rr);
}


double C12::
u  (const double & rr) const 
{
  double rr6 = rr * rr * rr * rr * rr * rr;
  return 1./(rr6 * rr6);
}

double C12::
up (const double & rr) const 
{
  double rr6 = rr * rr * rr * rr * rr * rr;
  return -12./(rr6 * rr6 * rr);
}

C12_Smooth::
C12_Smooth (const double & smooth_start,
	   const double & smooth_end)
    : sw(smooth_start, smooth_end)
{
}

double C12_Smooth::
u (const double & rr) const
{
  return (1. - sw.u(rr)) * c12.u(rr);
}

double C12_Smooth::
up (const double & rr) const
{
  return (1. - sw.u(rr)) * c12.up(rr) + ( - sw.up(rr)) * c12.u(rr);
}


int main(int argc, char * argv[])
{
  std::string ofile;
  double erf;
  double er;
  double rc, rsmooth;
  double ext = 1.2;
  double dr;
  double scaleEleStart, scaleEleEnd;
  double scaleVdwStart, scaleVdwEnd;
  double scale;
  double eleScale0, vdwScale0;
  double fmax = 1e30;
  
  po::options_description desc ("Allow options");
  desc.add_options()
    ("help,h", "print this message")
    ("er", po::value<double > (&er)->default_value (1.0), "reletive permitivity")
    ("erf", po::value<double > (&erf)->default_value (78.0), "dielectric constant")
    ("rc", po::value<double > (&rc)->default_value (1.3), "cut-off radius")
    ("r-smooth", po::value<double > (&rsmooth)->default_value (1.2), "cut-off radius")
    ("table-ext", po::value<double > (&ext)->default_value (1.2), "table extention")
    ("dr", po::value<double > (&dr)->default_value (0.002), "table step")
    ("scale", po::value<double > (&scale)->default_value (1.0), "scale")
    ("scale-ele-scale-0", po::value<double > (&eleScale0)->default_value (1.0), "ele scale0")
    ("scale-vdw-scale-0", po::value<double > (&vdwScale0)->default_value (1.0), "vdw scale0")
    ("scale-ele-start", po::value<double > (&scaleEleStart)->default_value (0.1), "scale start r")
    ("scale-ele-end", po::value<double > (&scaleEleEnd)->default_value (0.18), "scale end r")
    ("scale-vdw-start", po::value<double > (&scaleVdwStart)->default_value (0.3), "scale start r")
    ("scale-vdw-end", po::value<double > (&scaleVdwEnd)->default_value (0.6), "scale end r")
    ("output,o", po::value<std::string > (&ofile)->default_value ("table.xvg"), "the output table of potential");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }
  
  double rup = rc + ext;
  int nn = rup / dr + 1;
  
  FILE * fout = fopen (ofile.c_str(), "w");

  if (rsmooth < rc){
    C1_Smooth pc1 (rsmooth, rc);
    C6_Smooth pc6 (rsmooth, rc);
    C12_Smooth pc12 (rsmooth, rc);
    SwitchFunction psf1 (scaleEleStart, scaleEleEnd);
    SwitchFunction psf2 (scaleVdwStart, scaleVdwEnd);

    for (int ii = 0; ii < nn; ++ii){
      double rr = ii * dr;
      if (rr == 0){
	fprintf (fout, "%.10e  %.10e %.10e  %.10e %.10e  %.10e %.10e\n",
		 rr, 0., 0., 0., 0., 0., 0.);
      }
      else {
	double vc = (scale - eleScale0) * psf1.u(rr) * pc1.u(rr) + eleScale0 * pc1.u(rr);
	double fc = - ( (scale - eleScale0) * (psf1.up(rr) * pc1.u(rr) + psf1.u(rr) * pc1.up(rr)) + eleScale0 * pc1.up(rr));
	double vc6 = - ((scale - vdwScale0) * psf2.u(rr) * pc6.u(rr) + vdwScale0 * pc6.u(rr));
	double fc6 = ( (scale - vdwScale0) * (psf2.up(rr) * pc6.u(rr) + psf2.u(rr) * pc6.up(rr)) + vdwScale0 * pc6.up(rr));
	double vc12 = (scale - vdwScale0) * psf2.u(rr) * pc12.u(rr) + vdwScale0 * pc12.u(rr);
	double fc12 = - ( (scale - vdwScale0) * (psf2.up(rr) * pc12.u(rr) + psf2.u(rr) * pc12.up(rr)) + vdwScale0 * pc12.up(rr));
	if (fc > fmax) {
	  fc = 0;
	  vc = 0;
	}
	if (fc6 > fmax) {
	  fc6 = 0;
	  vc6 = 0;
	}
	if (fc12 > fmax) {
	  fc12 = 0;
	  vc12 = 0;
	}	
	fprintf (fout, "%.10e  %.10e %.10e  %.10e %.10e  %.10e %.10e\n",
		 rr, vc, fc, vc6, fc6, vc12, fc12
		 );
	// fprintf (fout, "%.10e  %.10e %.10e  %.10e %.10e  %.10e %.10e\n",
	// 	 rr,
	// 	 (scale - scale0) * psf1.u(rr) * pc1.u(rr) + scale0 * pc1.u(rr),
	// 	 - ( (scale - scale0) * (psf1.up(rr) * pc1.u(rr) + psf1.u(rr) * pc1.up(rr)) + scale0 * pc1.up(rr)),
	// 	 - ((scale - scale0) * psf2.u(rr) * pc6.u(rr) + scale0 * pc6.u(rr)),
	// 	 ( (scale - scale0) * (psf2.up(rr) * pc6.u(rr) + psf2.u(rr) * pc6.up(rr)) + scale0 * pc6.up(rr)),
	// 	 (scale - scale0) * psf2.u(rr) * pc12.u(rr) + scale0 * pc12.u(rr),
	// 	 - ( (scale - scale0) * (psf2.up(rr) * pc12.u(rr) + psf2.u(rr) * pc12.up(rr)) + scale0 * pc12.up(rr))
	// 	 );
      }
    }
  }
  else {
    C1 pc1;
    C6 pc6;
    C12 pc12;
    SwitchFunction psf1 (scaleEleStart, scaleEleEnd);
    SwitchFunction psf2 (scaleVdwStart, scaleVdwEnd);

    for (int ii = 0; ii < nn; ++ii){
      double rr = ii * dr;
      if (rr == 0){
	fprintf (fout, "%.10e  %.10e %.10e  %.10e %.10e  %.10e %.10e\n",
		 rr, 0., 0., 0., 0., 0., 0.);
      }
      else {
	double vc = (scale - eleScale0) * psf1.u(rr) * pc1.u(rr) + eleScale0 * pc1.u(rr);
	double fc = - ( (scale - eleScale0) * (psf1.up(rr) * pc1.u(rr) + psf1.u(rr) * pc1.up(rr)) + eleScale0 * pc1.up(rr));
	double vc6 = - ((scale - vdwScale0) * psf2.u(rr) * pc6.u(rr) + vdwScale0 * pc6.u(rr));
	double fc6 = ( (scale - vdwScale0) * (psf2.up(rr) * pc6.u(rr) + psf2.u(rr) * pc6.up(rr)) + vdwScale0 * pc6.up(rr));
	double vc12 = (scale - vdwScale0) * psf2.u(rr) * pc12.u(rr) + vdwScale0 * pc12.u(rr);
	double fc12 = - ( (scale - vdwScale0) * (psf2.up(rr) * pc12.u(rr) + psf2.u(rr) * pc12.up(rr)) + vdwScale0 * pc12.up(rr));
	if (fc > fmax) {
	  fc = 0;
	  vc = 0;
	}
	if (fc6 > fmax) {
	  fc6 = 0;
	  vc6 = 0;
	}
	if (fc12 > fmax) {
	  fc12 = 0;
	  vc12 = 0;
	}
	fprintf (fout, "%.10e  %.10e %.10e  %.10e %.10e  %.10e %.10e\n",
		 rr, vc, fc, vc6, fc6, vc12, fc12
		 );
	// fprintf (fout, "%.10e  %.10e %.10e  %.10e %.10e  %.10e %.10e\n",
	// 	 rr,
	// 	 (scale - scale0) * psf1.u(rr) * pc1.u(rr) + scale0 * pc1.u(rr),
	// 	 - ( (scale - scale0) * (psf1.up(rr) * pc1.u(rr) + psf1.u(rr) * pc1.up(rr)) + scale0 * pc1.up(rr)),
	// 	 - ((scale - scale0) * psf2.u(rr) * pc6.u(rr) + scale0 * pc6.u(rr)),
	// 	 ( (scale - scale0) * (psf2.up(rr) * pc6.u(rr) + psf2.u(rr) * pc6.up(rr)) + scale0 * pc6.up(rr)),
	// 	 (scale - scale0) * psf2.u(rr) * pc12.u(rr) + scale0 * pc12.u(rr),
	// 	 - ( (scale - scale0) * (psf2.up(rr) * pc12.u(rr) + psf2.u(rr) * pc12.up(rr)) + scale0 * pc12.up(rr))
	// 	 );
      }
    }
  }    
  
  fclose (fout);
  
  return 0;
}
