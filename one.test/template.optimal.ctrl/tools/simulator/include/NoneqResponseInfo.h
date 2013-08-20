#ifndef __NoneqResponseInfo_h_NONEQMSM_wanghan__
#define __NoneqResponseInfo_h_NONEQMSM_wanghan__

#include <vector>
#include "Defines.h"
#include "Distribution.h"
#include "Perturbation.h"

using namespace std;

class NoneqResponseInfo 
{
  double x0, x1, v0, v1;
  unsigned nx, nv;

  double beta;
  
  double dt;
  double noneqTime;
  double noneqCheckFeq;
  int numCheck;
  int checkNumFeq;

  int countNoneq;
  int countNoneqSeg;

  const Perturbation * ppert;
  int numMode;
  
  vector<double >				order0;		// time
  vector<double >				order0error;	// time
  vector<vector<double > >			order1;		// time x vec
  vector<vector<double > >			order1error;	// time x vec
  vector<vector<vector<double > > >		order2;		// time x mat  
  vector<double >				order0punish;	// time
  
  vector<double >				Gj;		// vec 
  vector<vector<double > >			Hjk;		// mat

  double					punish;
  int						ntraj;
public:
  NoneqResponseInfo ();
  // ~NoneqResponseInfo ();
public:
  void reinit (const double & beta,
	       const double & x0,
	       const double & x1,
	       const unsigned & nx,
	       const double & v0,
	       const double & v1,
	       const unsigned & nv,
	       const double & dt,
	       const double & noneqTime,
	       const double & noneqCheckFeq,
	       const Perturbation & pert);
  int inSet (const Dofs & x) ;
  void newTraj ();
  void depositMainTraj (const Dofs & oldx,
			const Dofs & newx,
			const double & sigma,
			const Dofs & dw);
  void average ();
  void collectLast () ;
  void collect () ;
public:
  // void calculate (const double & time,
  // 		  const Perturbation & pert1,
  // 		  double & dist,
  // 		  double & quench_dist,
  // 		  const int order = 2);
  const vector<double > &	   get_order0 () const {return order0;}
  const vector<double > &	   get_order0error () const {return order0error;}
  const vector<double > &	   get_order0punish () const {return order0punish;}
  const vector<vector<double > > & get_order1 () const {return order1;}
  const vector<vector<double > > & get_order1error () const {return order1error;}
}
    ;


#endif
