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
#include "NoneqResponseInfo.h"
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


int main(int argc, char * argv[])
{
  MPI::Init (argc, argv);
  int rank = COMM_WORLD.Get_rank();  
  
  double gamma;
  double kT;
  double T;
  double dt;
  double pertSt0;
  double pertSt1;
  double sigma;
  double kk;
  double nst;
  unsigned nstprint = 1;
  double warmTime;
  
  double noneqTime;
  double noneqCheckFeq;
  double branchFeq;
  unsigned noneqCheckNumFeq;
  unsigned branchNumFeq;
  double x0, x1;
  int order;
  string sfile, lfile;
  string ifile;
  double rateLag;
  int nTrajRecord;
  
  unsigned long seed;

  // Stopwatch sw_eq, sw_noneq, sw_quench;
  Stopwatch sw_total;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("dt,d", po::value<double > (&dt)->default_value(0.0001), "time step [ps]")
      ("nst,n", po::value<double > (&nst)->default_value(10000), "number of time step")
      ("print-feq,p", po::value<unsigned > (&nstprint)->default_value(10000), "print frequency")
      ("gamma,g", po::value<double > (&gamma)->default_value(1.), "gamma [ps^-1]")
      ("temperature,t", po::value<double > (&T)->default_value(300.), "temperature [K]")
      ("single-well-k,k",po::value<double > (&kk)->default_value(8.0), "k parameter of the double single potential [kJ/(mol nm^2)]")
      ("gaussian-sigma,a", po::value<double > (&sigma)->default_value(0.16), "sigma of Gaussian perturbation [nm]")
      ("branch-feq", po::value<double > (&branchFeq)->default_value(1.), "branch frequency [ps]")
      ("noneq-check-feq", po::value<double > (&noneqCheckFeq)->default_value(10.), "non-equilibrium branch check frequency [ps]")
      ("noneq-time", po::value<double > (&noneqTime)->default_value(200.), "non-equilibrium simulation time [ps]")
      ("pert-strength0",po::value<double > (&pertSt0)->default_value(1.0), "perturbation strength 0 [kJ/(mol nm)]")
      ("pert-strength1",po::value<double > (&pertSt1)->default_value(1.0), "perturbation strength 1 [kJ/(mol nm)]")
      ("warm-time", po::value<double > (&warmTime)->default_value(100.), "warm up time of perturbation [ps]")
      ("order", po::value<int > (&order)->default_value(1), "order of response")
      ("save-corr", po::value<string > (&sfile), "save correlation")
      ("load-corr", po::value<string > (&lfile), "load saved correlation")
      ("rate-lag", po::value<double > (&rateLag)->default_value(1.0), "lag of rate calculation")
      ("x-low", po::value<double > (&x0)->default_value (-2.0), "the lower bound of x range considered")
      ("x-up",  po::value<double > (&x1)->default_value ( 2.0), "the upper bound of x range considered")
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
    std::cout << "# sigma: " << sigma << " [nm]" << std::endl;
    std::cout << "# gamma: " << gamma << " [ps^-1]" << std::endl;
    std::cout << "# pert st0 (reference)   : " << pertSt0 << " [kJ/(mol nm)]" << std::endl;
    std::cout << "# pert st1 (perturbation): " << pertSt1 << " [kJ/(mol nm)]" << std::endl;
    std::cout << "# branch every: " << branchNumFeq << " steps" << std::endl;
    std::cout << "# noneq check Feq: " << noneqCheckFeq << " [ps]" << std::endl;
    std::cout << "# noneq check every: " << noneqCheckNumFeq << " steps" << std::endl;
    std::cout << "# noneq time: " << noneqTime << " [ps]" << std::endl;
    std::cout << "# rateLag: " << rateLag << " [ps]" << std::endl;
    std::cout << "# xrange: [ " << x0 << " , " << x1 << " ] " << std::endl;
    std::cout << "###################################################" << std::endl;  
  }

  unsigned numNoneqCheck = int((noneqTime + 0.5 * noneqCheckFeq) / noneqCheckFeq) + 1;

  sw_total.start();

  // initial sets
  SingleWell dw (kk);

  int count = 0;
  int countBranch = 0;
  double time = 0;

  PertGaussian pert  (pertSt0, 0., sigma, warmTime);
  PertGaussian pert1 (pertSt1, 0., sigma, warmTime);

  EulerMaruyama inte (gamma, kT, dt,
		      NULL,
		      dynamic_cast<Force *> (&dw),
		      seed + rank*2);
  EulerMaruyama noneqInte (gamma, kT, dt,
			   dynamic_cast<Perturbation *> (&pert),
			   dynamic_cast<Force *> (&dw),
			   seed + rank*2 + 1);
  double inteSigma = noneqInte.getSigma();
  if (rank == 0){
    std::cout << "# inte sigma is: " << inteSigma << std::endl;
  }
  Dofs xx;
  xx.xx[0] = 0.;
  xx.vv[0] = 0.;
  nTrajRecord = int ((rateLag + 0.5 * dt) / dt);
  nTrajRecord ++;
  Traj traj (nTrajRecord);
    
  NoneqResponseInfo resInfo;
  resInfo.reinit (x0, x1, dt, noneqTime, noneqCheckFeq, pert);

  if (!vm.count("load-corr")){
    for (double ii = 0.; ii < nst+0.1; ii += 1.){
      // sw_eq.start();
      inte.step (xx, 0.);
      // sw_eq.stop();
      count ++;
      countBranch ++;
      time += dt;
      if (int(nstprint) == count && rank == 0){
	count = 0;
	printf ("# %f %f %f    \r", time, xx.xx[0], xx.vv[0]);
	fflush (stdout);
      }
      if (countBranch == int(branchNumFeq)){
	countBranch = 0;
	Dofs branchXX (xx);
	Dofs branchXX_old (xx);
	traj.clear ();
	traj.push_back (xx);
	resInfo.newTraj ();
	// branching
	for (double ttNoneq = 0.; ttNoneq < noneqTime-0.5*dt; ttNoneq += dt){
	  // printf ("%f %d %d\n", ttNoneq, countNoneqCheck, noneqCheckNumFeq);
	  branchXX_old = branchXX;
	  // sw_noneq.start();
	  noneqInte.step (branchXX, ttNoneq);
	  // sw_noneq.stop();
	  Dofs dw = noneqInte.getDw ();
	  traj.push_back (branchXX);
	  // if ( int((ttNoneq + 1.5 * dt) / dt) % int ( (rateLag+0.5 * dt) / dt) == 0){
	  //   printf ("tt: %f posi is %.16f   %.16f\n", ttNoneq+dt, branchXX.xx[0], branchXX_old.xx[0]);
	  // }
	  resInfo.depositMainTraj (traj, inteSigma, dw);
	}
      }
    }
    resInfo.average ();
    resInfo.collect ();
    
    if (vm.count("save-corr")){
      if (rank == 0) resInfo.save (sfile);
    }
  }
  else {
    resInfo.load (lfile);
  }  
    
  if (rank == 0){
    char tmpfilename[1024];
    sprintf (tmpfilename, "state.out");
    FILE *fp = fopen (tmpfilename, "w");
    for (unsigned ii = 0; ii < resInfo.get_order0().size(); ++ii){
      fprintf (fp, "%e   %e\n",
	       noneqCheckFeq * ii,
	       resInfo.get_order0()[ii]);
    }
    fclose (fp);
  }
  if (rank == 0){
    char tmpfilename[1024];
    sprintf (tmpfilename, "state.resp.out");
    FILE *fp = fopen (tmpfilename, "w");
    for (unsigned ii = 0; ii < numNoneqCheck; ++ii){
      double nowTime = ii * noneqCheckFeq;
      double rate;
      resInfo.calculate (nowTime, pert1, rate, order);
      fprintf (fp, "%e   %e\n",
	       nowTime, rate);
    }
    fclose (fp);
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


