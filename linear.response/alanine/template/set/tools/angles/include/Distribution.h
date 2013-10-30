#ifndef __Distribution_h_NONEQMSM_wanghan__
#define __Distribution_h_NONEQMSM_wanghan__

#include "Defines.h"
#include <vector>
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

using namespace std;

class Distribution_1d
{
public:
  double x0, x1, v0, v1;
  double hx, hv;
  double valuepp;
  unsigned nx, nv;
  vector<vector<double > > values;
  vector<vector<double > > backup_values;
  double backup_number;
  double backup_unbacked_count;
  vector<double > gridx;
  vector<double > gridv;
  double nframe;
public:
  // Distribution_1d () {};
  Distribution_1d (const double x0 = -180,
		   const double x1 = 180.,
		   const unsigned nx = 36,
		   const double v0 = -180,
		   const double v1 = 180.,
		   const unsigned nv = 36);
  void reinit (const double & x0,
	       const double & x1,
	       const unsigned & nx,
	       const double & v0,
	       const double & v1,
	       const unsigned & nv);
public:
  void clear ();
  void deposite (const double & psi,
		 const double & phi);
  void deposite (const double & psi,
		 const double & phi,
		 const double & scale);
  void average ();
  void calTraj ();
  void print_x  (const string & filename) const;
  void print_xv (const string & filename) const;
  void print_x  (FILE * fp) const;
  void print_xv (FILE * fp) const;
  void substract (const Distribution_1d & d);
  void add (const double & scalor,
	    const Distribution_1d & d);
  double getNframe () const {return nframe;}
public:
  void save (FILE * fp) const;
  bool load (FILE * fp);
}
    ;

class MetastableSet 
{
private:
  double phi_b, phi_e;
  double psi_b, psi_e;
public:
  MetastableSet (const double & phb,
		 const double & phe,
		 const double & psb,
		 const double & pse);
  bool inSet (const double & phi,
	      const double & psi) const;
}
    ;



#endif
