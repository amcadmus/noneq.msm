#ifndef __Distribution_h_NONEQMSM_wanghan__
#define __Distribution_h_NONEQMSM_wanghan__

#include "Defines.h"


class Distribution_1d
{
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
  void average ();
  void calTraj ();
  void print_x  (const string & filename) const;
  void print_xv (const string & filename) const;
  void print_x  (FILE * fp) const;
  void print_xv (FILE * fp) const;
}
    ;

inline void Distribution_1d::
deposite (const Dofs & dof)
{
  int ix = (dof.xx[0] - x0) / hx;
  int iv = (dof.vv[0] - v0) / hv;
  if ((ix < 0 || ix >= int(nx)) || (iv < 0 || iv >= int(nv))){
    return;
  }
  else {
    nframe += 1.;
    backup_unbacked_count += 1.;
    values[ix][iv] += valuepp;
    if (backup_unbacked_count >= backup_number){
      // std::cout << "backuped" << std::endl;
      for (unsigned ii = 0; ii < nx; ++ii){
	for (unsigned jj = 0; jj < nv; ++jj){
	  backup_values[ii][jj] += values[ii][jj];
	  values[ii][jj] = 0.;
	}
      }
      backup_unbacked_count = 0.;
    }
  }
}



#endif
