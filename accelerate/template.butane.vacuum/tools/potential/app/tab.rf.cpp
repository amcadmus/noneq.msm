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

namespace po = boost::program_options;
using namespace std;

class ReactionField 
{
  double er, erf, rc, rc3;
public:
  ReactionField (const double & er_,
		 const double & erf_,
		 const double & rc_);
  double u  (const double & rr) const;
  double up (const double & rr) const;
}
    ;

class SwitchFunction
{
  double start, end, inter;
public:
  SwitchFunction (const double & start_,
		  const double & end_);
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

class C12
{
public:
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
    return 0.5 + 0.5 * cos ( (rr - start) / inter * M_PI + M_PI );
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
    return - 0.5 * sin ( (rr - start) / inter * M_PI + M_PI ) * M_PI / inter;
  }
}


ReactionField::
ReactionField (const double & er_,
	       const double & erf_,
	       const double & rc_)
    : er(er_), erf(erf_), rc(rc_)
{
  rc3 = rc * rc * rc;
}

double ReactionField::
u  (const double & rr) const 
{
  if (rr < rc) {
    return
	1. / (er * rr) 
	* ( 1. + (erf - er) / (2. * erf + er) * rr * rr * rr / rc3 )
	- 1. / (er * rc)
	* ( 3. * erf / (2. * erf + er) );
  }
  else {
    return 0.;
  }
}

double ReactionField::
up (const double & rr) const 
{
  if (rr < rc) {
    return
	- 1. / (er * rr * rr) 
	* ( 1. + (erf - er) / (2. * erf + er) * rr * rr * rr / rc3 )
	+ 1. / (er * rr)
	* ( (erf - er) / (2. * erf + er) * 3. * rr * rr / rc3 );
  }
  else {
    return 0.;
  }
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


int main(int argc, char * argv[])
{
  std::string ofile;
  double erf;
  double er;
  double rc;
  double ext = 1.2;
  double dr;
  double scaleEleStart, scaleEleEnd;
  double scaleVdwStart, scaleVdwEnd;
  double scale;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("er", po::value<double > (&er)->default_value (1.0), "reletive permitivity")
      ("erf", po::value<double > (&erf)->default_value (78.0), "dielectric constant")
      ("rc", po::value<double > (&rc)->default_value (1.3), "cut-off radius")
      ("table-ext", po::value<double > (&ext)->default_value (1.2), "table extention")
      ("dr", po::value<double > (&dr)->default_value (0.002), "table step")
      ("scale", po::value<double > (&scale)->default_value (1.0), "scale")
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
  
  ReactionField prf (er, erf, rc);
  C6 pc6;
  C12 pc12;
  SwitchFunction psf1 (scaleEleStart, scaleEleEnd);
  SwitchFunction psf2 (scaleVdwStart, scaleVdwEnd);

  double rup = rc + ext;
  int nn = rup / dr + 1;
  
  FILE * fout = fopen (ofile.c_str(), "w");
  for (int ii = 0; ii < nn; ++ii){
    double rr = ii * dr;
    if (rr == 0){
      fprintf (fout, "%f  %f %f  %f %f  %f %f\n",
	       rr, 0., 0., 0., 0., 0., 0.);
    }
    else {
      fprintf (fout, "%f  %f %f  %f %f  %f %f\n",
	       rr,
	       (scale - 1) * psf1.u(rr) * prf.u(rr) + prf.u(rr),
	       - ( (scale - 1) * (psf1.up(rr) * prf.u(rr) + psf1.u(rr) * prf.up(rr)) + prf.up(rr)),
	       - ((scale - 1) * psf2.u(rr) * pc6.u(rr) + pc6.u(rr)),
	       ( (scale - 1) * (psf2.up(rr) * pc6.u(rr) + psf2.u(rr) * pc6.up(rr)) + pc6.up(rr)),
	       (scale - 1) * psf2.u(rr) * pc12.u(rr) + pc12.u(rr),
	       - ( (scale - 1) * (psf2.up(rr) * pc12.u(rr) + psf2.u(rr) * pc12.up(rr)) + pc12.up(rr))
	  );
    }
  }
  fclose (fout);
  
  return 0;
}
