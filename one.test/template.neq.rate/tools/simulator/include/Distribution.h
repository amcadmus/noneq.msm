#ifndef __Distribution_h_NONEQMSM_wanghan__
#define __Distribution_h_NONEQMSM_wanghan__

#include "Defines.h"


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
  Distribution_1d () {};
  Distribution_1d (const double & x0,
		   const double & x1,
		   const unsigned & nx,
		   const double & v0,
		   const double & v1,
		   const unsigned & nv);
  void reinit (const double & x0,
	       const double & x1,
	       const unsigned & nx,
	       const double & v0,
	       const double & v1,
	       const unsigned & nv);
public:
  void clear ();
  void deposite (const Dofs & dof);
  void deposite (const Dofs & dof, const double & scale);
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
  void load (FILE * fp);
}
    ;




#endif
