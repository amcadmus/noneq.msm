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
  
  double dt;
  double noneqTime;
  double noneqCheckFeq;
  int numCheck;
  int checkNumFeq;
  double quenchTime;
  int quenchNumStep;

  int countNoneq;
  int countNoneqSeg;
  int countQuench;

  const Perturbation * ppert;
  int numMode;
  
  vector<Distribution_1d >			order0;		// time
  vector<vector<Distribution_1d > >		order1;		// time x vec
  vector<vector<vector<Distribution_1d > > >	order2;		// time x mat  
  
  vector<Distribution_1d >			quench_order0;		// time
  vector<Distribution_1d >			quench_order1_term1;	// time
  vector<vector<Distribution_1d > >		quench_order1_term2;	// time x vec
  vector<Distribution_1d >			quench_order2_term1;	// time
  vector<vector<vector<Distribution_1d > > >	quench_order2_term2;	// time x mat
  vector<vector<Distribution_1d > >		quench_order2_term3;	// time x vec
  
  vector<double >				Gj;		// vec 
  vector<vector<double > >			Hjk;		// mat 
  double					G0;		// 1
  double					H00;		// 1

public:
  NoneqResponseInfo ();
  // ~NoneqResponseInfo ();
public:
  void reinit (const double & x0,
	       const double & x1,
	       const unsigned & nx,
	       const double & v0,
	       const double & v1,
	       const unsigned & nv,
	       const double & dt,
	       const double & noneqTime,
	       const double & noneqCheckFeq,
	       const double & quenchTime,
	       const Perturbation & pert);
  void newTraj ();
  void newQuenchTraj ();
  void depositMainTraj (const Dofs & oldx,
			const Dofs & newx,
			const double & sigma,
			const Dofs & dw);
  void depositQuenchTraj (const Dofs & oldx,
			  const Dofs & newx,
			  const double & sigma,
			  const Dofs & dw);
  void average ();
public:
  void calculate (const double & time,
		  const Perturbation & pert1,
		  Distribution_1d & dist,
		  Distribution_1d & quench_dist,
		  const int order = 2);
public:
  void save (const string & filename) const;
  void load (const string & filename);
}
    ;


#endif
