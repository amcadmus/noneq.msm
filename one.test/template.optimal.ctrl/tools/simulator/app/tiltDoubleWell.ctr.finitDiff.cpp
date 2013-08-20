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

#include "Force.h"
#include "Perturbation.h"
#include "Integrator.h"
#include "Distribution.h"
#include "NoneqInfo.h"
#include "Stopwatch.h"

#include <boost/program_options.hpp>
#include "RandomGenerator.h"
#include <mpi.h>

namespace po = boost::program_options;
using namespace MPI;

void readFile (const string & filename,
	       std::vector<double > & tt,
	       std::vector<double > & value)
{
  tt.clear ();
  value.clear ();
  FILE * fp = fopen (filename.c_str(), "r");
  if (fp == NULL){
    cerr << "cannot open file " << filename << endl;
    exit (1);
  }

  double tmp1, tmp2;
  while (fscanf (fp, "%lf%lf", &tmp1, &tmp2) == 2){
    tt.push_back (tmp1);
    value.push_back (tmp2);
  }
  
  fclose (fp);
}


void cal_order1punish (const double & beta,
		       const vector<double > & tt,		// size num mode
		       const vector<double > & ttvalue,
		       vector<double > & order1punish)
{
  const std::vector<double > & vv (ttvalue);
  order1punish.resize (tt.size());
  order1punish[0] = 0.5 * beta * (2./3. * vv[0] + 1./3. * vv[1]) * (tt[1] - tt[0]);
  int nn = tt.size() - 1;
  order1punish[nn] = 0.5 * beta * (1./3. * vv[nn-1] + 2./3. * vv[nn]) * (tt[nn] - tt[nn-1]);
  
  for (unsigned ii = 1; ii < tt.size() - 1; ++ii){
    order1punish[ii] = 0.5 * beta * ((1./3. * vv[ii-1] + 2./3. * vv[ii]) * (tt[ii] - tt[ii-1]) +
				     (1./3. * vv[ii+1] + 2./3. * vv[ii]) * (tt[ii+1] - tt[ii]) );
  }
}

double runTraj (const vector<double > & tt,		// size num mode
		const vector<double > & ttvalue,
		vector<double > & order0,		// size noneq check
		vector<double > & order0error,		// size noneq check
		vector<double > & order0punish,		// size noneq check
		const DoubleWell & dw,
		const double & dt,
		const double & gamma,
		const double & kT,
		const double & nst,
		const double & x0,
		const double & x1,
		const double & beta,
		const double & noneqTime,
		const double & noneqCheckFeq,
		const double & branchFeq,
		const int & nstprint,
		const long unsigned int seed)
{
  int count = 0;
  int countBranch = 0;
  double time = 0;
  static int myiter = 0;
  int rank = COMM_WORLD.Get_rank();
  int branchNumFeq = (branchFeq+0.5*dt) / dt;

  PertConstTiltTable pert (tt, ttvalue);

  EulerMaruyama inte (gamma, kT, dt,
		      NULL,
		      dynamic_cast<const Force *> (&dw),
		      myiter + seed + rank*2);
  EulerMaruyama noneqInte (gamma, kT, dt,
			   dynamic_cast<Perturbation *> (&pert),
			   dynamic_cast<const Force *> (&dw),
			   myiter + seed + rank*2 + 1);
  myiter ++;
  
  double inteSigma = noneqInte.getSigma();

  Dofs xx;
  xx.xx[0] = 0.;
  xx.vv[0] = 0.;

  NoneqInfo noneqInfo;
  noneqInfo.reinit (beta, x0, x1, dt, noneqTime, noneqCheckFeq, pert);

  for (double ii = 0.; ii < nst+0.1; ii += 1.){
    inte.step (xx, 0.);
    count ++;
    countBranch ++;
    time += dt;
    if (int(nstprint) == count && rank == 0){
      count = 0;
      printf ("# traj time: %f %f %f    \r", time, xx.xx[0], xx.vv[0]);
      fflush (stdout);
    }
    if (countBranch == int(branchNumFeq)){
      countBranch = 0;
      Dofs branchXX (xx);
      Dofs branchXX_old (xx);
      noneqInfo.newTraj ();
      // branching
      for (double ttNoneq = 0.; ttNoneq < noneqTime-0.5*dt; ttNoneq += dt){
	// printf ("%f %d %d\n", ttNoneq, countNoneqCheck, noneqCheckNumFeq);
	branchXX_old = branchXX;
	noneqInte.step (branchXX, ttNoneq);
	Dofs dw = noneqInte.getDw ();
	noneqInfo.depositMainTraj (branchXX_old, branchXX, inteSigma, dw);
      }
    }
  }
  noneqInfo.average ();
  noneqInfo.collectLast ();
  
  order0 = noneqInfo.get_order0();
  order0error = noneqInfo.get_order0error();
  order0punish = noneqInfo.get_order0punish();
  if (rank == 0) printf ("\n# rank: %d order0: %e\n", rank, order0.back());

  return order0.back() + order0punish.back();
}



int main(int argc, char * argv[])
{
  MPI::Init (argc, argv);
  int rank = COMM_WORLD.Get_rank();  
  
  double gamma;
  double kT;
  double T;
  double dt;
  double pertSt0;
  double aa;
  double kk;
  double nst;
  unsigned nstprint = 1;
  
  double noneqTime;
  double noneqCheckFeq;
  double branchFeq;
  unsigned noneqCheckNumFeq;
  unsigned branchNumFeq;
  double x0, x1, v0, v1;
  unsigned nx, nv;
  int order;
  double beta;
  double timeResolution;
  unsigned nTimeFrame;
  string sfile, lfile;
  double gradientDescentStep;
  double finiteDiffStep;
  string ifile;
  
  unsigned long seed;

  // Stopwatch sw_eq, sw_noneq, sw_quench;
  Stopwatch sw_total;
  
  po::options_description desc ("Allow options");
  desc.add_options()
    ("help,h", "print this message")
    ("dt,d", po::value<double > (&dt)->default_value(0.0001), "time step [ps]")
    ("beta,b", po::value<double > (&beta)->default_value(1.0), "the punishment constant, unknown unit...")
    ("nst,n", po::value<double > (&nst)->default_value(10000), "number of time step")
    ("print-feq,p", po::value<unsigned > (&nstprint)->default_value(10000), "print frequency")
    ("gamma,g", po::value<double > (&gamma)->default_value(1.), "gamma [ps^-1]")
    ("temperature,t", po::value<double > (&T)->default_value(300.), "temperature [K]")
    ("double-well-k,k", po::value<double > (&kk)->default_value(8.0), "k parameter of the double well potential [kJ/(mol nm^2)]")
    ("double-well-a,a", po::value<double > (&aa)->default_value(1.0), "a parameter of the double well potential [nm]")
    ("branch-feq", po::value<double > (&branchFeq)->default_value(1.), "branch frequency [ps]")
    ("noneq-check-feq", po::value<double > (&noneqCheckFeq)->default_value(10.), "non-equilibrium branch check frequency [ps]")
    ("noneq-time", po::value<double > (&noneqTime)->default_value(200.), "non-equilibrium simulation time [ps]")
    ("time-resolution", po::value<double > (&timeResolution)->default_value(10.), "time resolution of the ctrl [ps]")      
    ("pert-strength0",po::value<double > (&pertSt0)->default_value(1.0), "perturbation strength 0 [kJ/(mol nm)]")
    ("order", po::value<int > (&order)->default_value(2), "order of response")
    // ("save-corr", po::value<string > (&sfile), "save correlation")
    // ("load-corr", po::value<string > (&lfile), "load saved correlation")
    ("x-low", po::value<double > (&x0)->default_value (-2.0), "the lower bound of x range considered")
    ("x-up",  po::value<double > (&x1)->default_value ( 2.0), "the upper bound of x range considered")
    ("v-low", po::value<double > (&v0)->default_value (-8.0), "the lower bound of v range considered")
    ("v-up",  po::value<double > (&v1)->default_value ( 8.0), "the upper bound of v range considered")
    ("x-grid", po::value<unsigned > (&nx)->default_value (50), "the number of grid point of x")
    ("v-grid", po::value<unsigned > (&nv)->default_value (50), "the number of grid point of v")
    ("gradient-descent-step", po::value<double > (&gradientDescentStep)->default_value (0.1), "step size of GD")
    ("finite-diff-step", po::value<double > (&finiteDiffStep)->default_value (0.5), "step size of finite diff to calculate the gradient")
    ("init,i", po::value<string > (&ifile), "initial guess of GD")
    ("seed",po::value<unsigned long > (&seed)->default_value(1), "random seed");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  RandomGenerator_MT19937::init_genrand (seed+rank);
  kT = T * 1.38 * 6.02 * 1e-3;
  branchNumFeq = (branchFeq+0.5*dt) / dt;
  noneqCheckNumFeq = (noneqCheckFeq+0.5*dt) / dt;

  if (rank == 0){
    std::cout << "###################################################" << std::endl;
    std::cout << "# T: " << T << " [K]" << std::endl;
    std::cout << "# kT: " << kT << " [kJ/mol]" << std::endl;
    std::cout << "# k: " << kk << " [kJ/(mol nm^4)]" << std::endl;
    std::cout << "# a: " << aa << " [nm]" << std::endl;
    std::cout << "# barrier: " << 0.5 * kk * aa*aa*aa*aa << " [kJ/mol]" << std::endl;
    std::cout << "# gamma: " << gamma << " [ps^-1]" << std::endl;
    std::cout << "# pert st0 (init)   : " << pertSt0 << " [kJ/(mol nm)]" << std::endl;
    std::cout << "# branch Feq: " << branchFeq << " [ps]" << std::endl;
    std::cout << "# branch every: " << branchNumFeq << " steps" << std::endl;
    std::cout << "# noneq check Feq: " << noneqCheckFeq << " [ps]" << std::endl;
    std::cout << "# noneq check every: " << noneqCheckNumFeq << " steps" << std::endl;
    std::cout << "# noneq time: " << noneqTime << " [ps]" << std::endl;
    std::cout << "# xrange: [ " << x0 << " , " << x1 << " ] " << std::endl;
    std::cout << "# vrange: [ " << v0 << " , " << v1 << " ] " << std::endl;
    std::cout << "# nx: " << nx << std::endl;
    std::cout << "# nv: " << nv << std::endl;
    std::cout << "###################################################" << std::endl;  
  }
  

  // initial sets
  DoubleWell dw (kk, aa);

  nTimeFrame = unsigned((noneqTime + 0.5 * timeResolution) / timeResolution) + 1;
  vector<double > tt (nTimeFrame);
  vector<double > ttvalue (nTimeFrame);
  if (vm.count("init") ){
    readFile (ifile, tt, ttvalue);
    if (tt.size() != nTimeFrame){
      cerr << "inconsistent time frame with init input file\n";
      return 1;
    }
  }
  else {
    for (unsigned ii = 0; ii < nTimeFrame; ++ii){
      tt[ii] = timeResolution * ii;
      ttvalue[ii] = pertSt0 * tt[ii] / noneqTime;
    }
  }
  
  sw_total.start();

  for (unsigned iter = 0; iter < 100 ; ++iter){    
    double thisvalue;
    vector<double > order0, order0error, order0punish;
    vector<double > order1, order1error1, order1error2, order1punish;
    
    cal_order1punish (beta, tt, ttvalue, order1punish);

    thisvalue = runTraj (tt, ttvalue,
			 order0, order0error, order0punish,
			 dw, dt, gamma, kT, nst, x0, x1, beta, noneqTime, noneqCheckFeq, branchFeq, nstprint, seed);    
    if (rank == 0) printf ("# finish value at the point				\n");
    
    if (rank == 0) {
      FILE * fp;
      char tmpfilename[2000];
      sprintf (tmpfilename, "state.step%03d.out", iter+1);
      fp = fopen (tmpfilename, "w");
      for (unsigned ii = 0; ii < order0.size(); ++ii){
    	fprintf (fp, "%e   %e %e %e   %e\n",
    		 noneqCheckFeq * ii,
    		 order0[ii],
    		 order0punish[ii],
    		 order0[ii] + order0punish[ii],
		 2. * sqrt(order0error[ii] - order0[ii] * order0[ii])
	    );
      }
      fclose (fp);
    }
    
    order1.resize(nTimeFrame);
    order1error1.resize(nTimeFrame);
    order1error2.resize(nTimeFrame);
    for (unsigned ii = nTimeFrame - 1; ii > 0; --ii){
      if (rank == 0) printf ("# cal %d th grad					\n", ii);
      vector<double > order0a, order0errora, order0punisha;
      vector<double > order0b, order0errorb, order0punishb;
      vector<double > ttvaluea (ttvalue);
      vector<double > ttvalueb (ttvalue);
      ttvaluea[ii] += finiteDiffStep;
      ttvalueb[ii] -= finiteDiffStep;
      runTraj (tt, ttvaluea,
	       order0a, order0errora, order0punisha,
	       dw, dt, gamma, kT, nst, x0, x1, beta, noneqTime, noneqCheckFeq, branchFeq, nstprint, seed);
      runTraj (tt, ttvalueb,
	       order0b, order0errorb, order0punishb,
	       dw, dt, gamma, kT, nst, x0, x1, beta, noneqTime, noneqCheckFeq, branchFeq, nstprint, seed);
      order1[ii] = (order0a.back() - order0b.back()) / (2. * finiteDiffStep);
      order1error1[ii] = order0errora.back() - order0a.back() * order0a.back();
      order1error2[ii] = order0errorb.back() - order0b.back() * order0b.back();
      if (rank == 0) printf ("order1[%d] = %e ,  order1punish[%d] = %e\n", ii, order1[ii], ii, order1punish[ii]);
    }


    if (rank == 0) {
      FILE * fp;
      char tmpfilename[2000];
      sprintf (tmpfilename, "order1.step%03d.out", iter+1);
      fp = fopen (tmpfilename, "w");
      fprintf (fp, "# time   order1   order1punish   error1 error2\n");
      for (unsigned ii = 0; ii < order1.size(); ++ii){
    	fprintf (fp, "%e   %e   %e   %e %e %e\n",
    		 tt[ii],
    		 order1[ii],
    		 order1punish[ii],
		 2. * sqrt(order1error1[ii]),
		 2. * sqrt(order1error2[ii]),
		 2. * sqrt(order1error1[ii] + order1error2[ii]) / (2. * finiteDiffStep)
	    );
	  }
      fclose (fp);
    }

    for (unsigned ii = 0; ii < nTimeFrame; ++ii){
      ttvalue[ii] -= gradientDescentStep * (order1[ii] + order1punish[ii]);
    }    
	
    if (rank == 0) {
      FILE * fp;
      char tmpfilename[2000];
      sprintf (tmpfilename, "ctr.step%03d.out", iter+1);
      fp = fopen (tmpfilename, "w");
      for (unsigned ii = 0; ii < order1.size(); ++ii){
    	fprintf (fp, "%e   %e\n",
    		 tt[ii],
    		 ttvalue[ii]);
	  }
      fclose (fp);
    }
    // if (rank == 0){
    //   printf ("\n");
    //   printf ("step %d\n", iter);
    //   printf ("value of ctr:           ");
    //   for (unsigned ii = 0; ii < nTimeFrame; ++ii){
    // 	printf ("%e \t", ttvalue[ii]);
    //   }
    //   printf ("\n");
    //   printf ("endv %e \t endpunish: %e \t endtotal: %e\n",
    // 	      order0().back(), order0punish().back(),
    // 	      order0().back()+ order0punish().back()
    // 	      );
    //   printf ("value of end    order1: ");
    //   for (unsigned ii = 0; ii < order1().back().size(); ++ii){
    // 	printf ("%e \t", order1().back()[ii]);
    //   }
    //   printf ("\n");    
    //   printf ("value of punish order1: ");
    //   for (unsigned ii = 0; ii < resInfo.get_order1().back().size(); ++ii){
    // 	printf ("%e \t", order1punish[ii]);
    //   }
    //   printf ("\n");    
    // }
    
    // if (resInfo.get_order1().back().size() != nTimeFrame) {
    //   cerr << "rank: "
    // 	   << rank
    // 	   << ". problem: nTimeFrame and numMode do not match"
    // 	   << resInfo.get_order1().back().size() << " "
    // 	   << nTimeFrame
    // 	   << endl;
    // }


    
    // ttvalue.back() -= 0.5;

    // for (int ii = 0; ii < COMM_WORLD.Get_size(); ++ii){
    //   if (ii == rank){
    // 	printf ("rank: %03d value of ctr:", rank);
    // 	for (unsigned ii = 0; ii < nTimeFrame; ++ii){
    // 	  printf ("%e \t", ttvalue[ii]);
    // 	}
    // 	printf ("\n");
    //   }
    //   COMM_WORLD.Barrier();
    // }
    
    // if (rank == 0){
    //   printf ("value of ctr:           ");
    //   for (unsigned ii = 0; ii < nTimeFrame; ++ii){
    // 	printf ("%e \t", ttvalue[ii]);
    //   }
    //   printf ("\n");
    //   char tmpfilename[1024];
    //   sprintf (tmpfilename, "ctr.step%03d.out", iter+1);
    //   FILE * fp = fopen (tmpfilename, "w");
    //   for (unsigned ii = 0; ii < nTimeFrame; ++ii){
    // 	fprintf (fp, "%e %e\n", tt[ii], ttvalue[ii]);
    //   }
    //   fclose (fp);
    //   sprintf (tmpfilename, "state.step%03d.out", iter+1);
    //   fp = fopen (tmpfilename, "w");
    //   for (unsigned ii = 0; ii < resInfo.get_order0().size(); ++ii){
    // 	fprintf (fp, "%e   %e %e %e\n",
    // 		 noneqCheckFeq * ii,
    // 		 resInfo.get_order0()[ii],
    // 		 resInfo.get_order0punish()[ii],
    // 		 resInfo.get_order0()[ii] + resInfo.get_order0punish()[ii]);
    //   }
    //   fclose (fp);
    // }
  }
  

  sw_total.stop();

  if (rank == 0){
    cout << "# time static: user, real, sys" << endl;
    // cout << "eq inte:       "
    //      << sw_eq.user() << " "
    //      << sw_eq.real() << " "
    //      << sw_eq.system() << endl;
    // cout << "noneq inte:    "
    //      << sw_noneq.user() << " "
    //      << sw_noneq.real() << " "
    //      << sw_noneq.system() << endl;
    // cout << "quench inte:   "
    //      << sw_quench.user() << " "
    //      << sw_quench.real() << " "
    //      << sw_quench.system() << endl;
    // cout << "total inte:    "
    //      << sw_eq.user() + sw_noneq.user() + sw_quench.user() << " "
    //      << sw_eq.real() + sw_noneq.real() + sw_quench.real() << " "
    //      << sw_eq.system() + sw_noneq.system() + sw_quench.system() << endl;
    cout << "# total:    "
	 << sw_total.user() << " "
	 << sw_total.real() << " "
	 << sw_total.system() << endl;
  }
  
  Finalize ();
  return 0;
}


