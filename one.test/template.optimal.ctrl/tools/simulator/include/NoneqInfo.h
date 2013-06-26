#ifndef __NoneqInfo_h_NONEQMSM_wanghan__
#define __NoneqInfo_h_NONEQMSM_wanghan__

#include <vector>
#include "Defines.h"
#include "Distribution.h"
#include "Perturbation.h"

using namespace std;

class NoneqInfo 
{
  double x0, x1;

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
  vector<double >				order0punish;	// time
  
  double					punish;
  int						ntraj;
public:
  NoneqInfo ();
  // ~NoneqInfo ();
public:
  void reinit (const double & beta,
	       const double & x0,
	       const double & x1,
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
public:
  const vector<double > &	   get_order0 () const {return order0;}
  const vector<double > &	   get_order0punish () const {return order0punish;}
}
    ;


#endif
