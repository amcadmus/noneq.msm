#ifndef __NoneqResponseInfo_h_NONEQMSM_wanghan__
#define __NoneqResponseInfo_h_NONEQMSM_wanghan__

#include <vector>
#include <list>
#include "Defines.h"
#include "Distribution.h"
#include "Perturbation.h"
#include "Traj.h"

using namespace std;

class NoneqResponseInfo 
{
  double x0, x1;
  
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
  vector<vector<double > >			order1;		// time x vec
  vector<vector<vector<double > > >		order2;		// time x mat  
  
  vector<double >				Gj;		// vec 
  vector<vector<double > >			Hjk;		// mat

  int						ntraj;
public:
  NoneqResponseInfo ();
  // ~NoneqResponseInfo ();
public:
  void reinit (const double & x0,
	       const double & x1,
	       const double & dt,
	       const double & noneqTime,
	       const double & noneqCheckFeq,
	       const Perturbation & pert);
  int inSet (const Dofs & x) ;
  void newTraj ();
  double trajObservable (const Traj & traj);
  void depositMainTraj (const Traj & traj,
			const double & sigma,
			const Dofs & dw);
  void average ();
  void collect () ;
public:
  void calculate (const double & time,
		  const Perturbation & pert1,
		  double & rate,
		  const int order = 1);
  const vector<double > &	   get_order0 () const {return order0;}
  const vector<vector<double > > & get_order1 () const {return order1;}
public:
  void save (const string & filename) const;
  void load (const string & filename);
}
    ;


#endif
