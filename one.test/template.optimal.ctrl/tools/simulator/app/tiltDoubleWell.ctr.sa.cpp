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

double runTraj (const vector<double > & tt,		// size num mode
		const vector<double > & ttvalue,
		vector<double > & order0,		// size noneq check
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
  string ifile;
  double saNst;
  double saSigma;
  double saChangeMin;
  double saTmax;
  
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
      ("init,i", po::value<string > (&ifile), "initial guess of GD")
      ("saNst", po::value<double > (&saNst)->default_value (200), "")
      ("saTmax", po::value<double > (&saTmax)->default_value (1.0), "")
      ("saSigma", po::value<double > (&saSigma)->default_value (0.2), "")
      ("saChangeMin", po::value<double > (&saChangeMin)->default_value (0.05), "")
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

  double thisvalue;
  vector<double > order0, order0punish;
  thisvalue = runTraj (tt, ttvalue,
		       order0, order0punish,
		       dw, dt, gamma, kT, nst, x0, x1, beta, noneqTime, noneqCheckFeq, branchFeq, nstprint, seed);    

  for (unsigned iter = 0; iter < saNst + 0.5; ++iter){    

    double trialvalue;
    vector<double > trialttvalue (ttvalue);
    vector<double > trialOrder0, trialOrder0punish;

    double temperature = (1. - iter / double(saNst)) * saTmax;
    int changeII;
    double changeSize;
    if (rank == 0){
      changeII = int(RandomGenerator_MT19937::genrand_real2() * nTimeFrame);
      while ( (changeSize = RandomGenerator_MT19937::gaussian() * saSigma) < saChangeMin ) {
      }
    }
    COMM_WORLD.Bcast (&changeII, 1, MPI_INT, 0);
    COMM_WORLD.Bcast (&changeSize, 1, MPI_DOUBLE, 0);
    trialttvalue[changeII] += changeSize;
    
    trialvalue = runTraj (tt, trialttvalue,
			  trialOrder0, trialOrder0punish,
			  dw, dt, gamma, kT, nst, x0, x1, beta, noneqTime, noneqCheckFeq, branchFeq, nstprint, seed);    

    // if (rank == 0) printf ("# trialvalue: %e   %e %e\n", trialvalue, trialOrder0.back(), trialOrder0punish);
    printf ("# trialvalue: %e   %e %e\n", trialvalue, trialOrder0.back(), trialOrder0punish.back());
    int accept = 0;

    if (rank == 0){
      if (trialvalue < thisvalue){
	accept = 1;
      }
      else {
	if (RandomGenerator_MT19937::genrand_real1() < exp ( - (trialvalue - thisvalue) / temperature)){
	  accept = 1;
	}
      }
    }
    COMM_WORLD.Bcast (&accept, 1, MPI_INT, 0);

    if (accept){
      // if (rank == 0) printf ("# accepted trialvalue: %e   %e %e\n", trialvalue, trialOrder0.back(), trialOrder0punish);
      printf ("# accepted trialvalue: %e   %e %e\n", trialvalue, trialOrder0.back(), trialOrder0punish.back());
      ttvalue = trialttvalue;
      thisvalue = trialvalue;
      order0 = trialOrder0;
      order0punish = trialOrder0punish;
      
      if (rank == 0) {
	FILE * fp;
	char tmpfilename[2000];
	sprintf (tmpfilename, "state.step%03d.out", iter+1);
	fp = fopen (tmpfilename, "w");
	for (unsigned ii = 0; ii < order0.size(); ++ii){
	  fprintf (fp, "%e   %e %e %e\n",
		   noneqCheckFeq * ii,
		   order0[ii],
		   order0punish[ii],
		   order0[ii] + order0punish[ii]);
	}
	fclose (fp);
      }
      if (rank == 0) {
	FILE * fp;
	char tmpfilename[2000];
	sprintf (tmpfilename, "ctr.step%03d.out", iter+1);
	fp = fopen (tmpfilename, "w");
	for (unsigned ii = 0; ii < nTimeFrame; ++ii){
	  fprintf (fp, "%e   %e\n",
		   tt[ii],
		   ttvalue[ii]);
	}
	fclose (fp);
      }
    }    
  }
  

  sw_total.stop();

  if (rank == 0){
    cout << "# time static: user, real, sys" << endl;
    cout << "# total:    "
	 << sw_total.user() << " "
	 << sw_total.real() << " "
	 << sw_total.system() << endl;
  }
  
  Finalize ();
  return 0;
}


