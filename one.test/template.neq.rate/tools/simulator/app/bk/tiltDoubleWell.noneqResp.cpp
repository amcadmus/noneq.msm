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

namespace po = boost::program_options;

int main(int argc, char * argv[])
{
  double gamma;
  double kT;
  double T;
  double dt;
  double pertSt0, pertSt1;
  double aa;
  double kk;
  double nst;
  unsigned nstprint = 1;
  double warmTime;
  
  double noneqTime;
  double noneqCheckFeq;
  double branchFeq;
  double quenchTime;
  double quenchTemperture;
  double quenchkT;
  unsigned noneqCheckNumFeq;
  unsigned branchNumFeq;
  double x0, x1, v0, v1;
  unsigned nx, nv;
  int order;

  string sfile, lfile;
  
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
      ("double-well-k,k", po::value<double > (&kk)->default_value(8.0), "k parameter of the double well potential [kJ/(mol nm^2)]")
      ("double-well-a,a", po::value<double > (&aa)->default_value(1.0), "a parameter of the double well potential [nm]")
      ("branch-feq", po::value<double > (&branchFeq)->default_value(1.), "branch frequency [ps]")
      ("noneq-check-feq", po::value<double > (&noneqCheckFeq)->default_value(10.), "non-equilibrium branch check frequency [ps]")
      ("noneq-time", po::value<double > (&noneqTime)->default_value(200.), "non-equilibrium simulation time [ps]")
      ("quench-temperature", po::value<double > (&quenchTemperture)->default_value(150.), "quench temperature [K]")
      ("quench-time", po::value<double > (&quenchTime)->default_value(1.), "quench time [ps]")
      ("pert-strength0",po::value<double > (&pertSt0)->default_value(0.0), "perturbation strength 0 [kJ/(mol nm)]")
      ("pert-strength1",po::value<double > (&pertSt1)->default_value(1.0), "perturbation strength 1 [kJ/(mol nm)]")
      ("warm-time", po::value<double > (&warmTime)->default_value(100.), "warm up time of perturbation [ps]")
      ("order", po::value<int > (&order)->default_value(2), "order of response")
      ("save-corr", po::value<string > (&sfile), "save correlation")
      ("load-corr", po::value<string > (&lfile), "load saved correlation")
      ("x-low", po::value<double > (&x0)->default_value (-2.0), "the lower bound of x range considered")
      ("x-up",  po::value<double > (&x1)->default_value ( 2.0), "the upper bound of x range considered")
      ("v-low", po::value<double > (&v0)->default_value (-8.0), "the lower bound of v range considered")
      ("v-up",  po::value<double > (&v1)->default_value ( 8.0), "the upper bound of v range considered")
      ("x-grid", po::value<unsigned > (&nx)->default_value (50), "the number of grid point of x")
      ("v-grid", po::value<unsigned > (&nv)->default_value (50), "the number of grid point of v")
      ("seed",po::value<unsigned long > (&seed)->default_value(1), "random seed");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  RandomGenerator_MT19937::init_genrand (seed);
  kT = T * 1.38 * 6.02 * 1e-3;
  quenchkT = quenchTemperture * 1.38 * 6.02 * 1e-3;
  branchNumFeq = (branchFeq+0.5*dt) / dt;
  noneqCheckNumFeq = (noneqCheckFeq+0.5*dt) / dt;
  
  std::cout << "###################################################" << std::endl;
  std::cout << "# T: " << T << " [K]" << std::endl;
  std::cout << "# kT: " << kT << " [kJ/mol]" << std::endl;
  std::cout << "# k: " << kk << " [kJ/(mol nm^4)]" << std::endl;
  std::cout << "# a: " << aa << " [nm]" << std::endl;
  std::cout << "# barrier: " << 0.5 * kk * aa*aa*aa*aa << " [kJ/mol]" << std::endl;
  std::cout << "# gamma: " << gamma << " [ps^-1]" << std::endl;
  std::cout << "# pert st0 (reference)   : " << pertSt0 << " [kJ/(mol nm)]" << std::endl;
  std::cout << "# pert st1 (perturbation): " << pertSt1 << " [kJ/(mol nm)]" << std::endl;
  std::cout << "# warm time: " << warmTime << " [ps]" << std::endl;
  std::cout << "# quench T: " << quenchTemperture << " [K]" << std::endl;
  std::cout << "# quench kT: " << quenchkT << " [kJ/mol]" << std::endl;
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

  DoubleWell dw (kk, aa);
  PertConstTilt pert (pertSt0, warmTime);
  PertConstTilt pert1 (pertSt1, warmTime);
  // PertConstTilt pert (0.);
  
  EulerMaruyama inte (gamma, kT, dt,
		      NULL,
		      dynamic_cast<Force *> (&dw),
		      seed);
  EulerMaruyama quenchInte (gamma, quenchkT, dt,
			    // NULL,
			    dynamic_cast<Perturbation *> (&pert),
			    dynamic_cast<Force *> (&dw),
			    seed+1);
  EulerMaruyama noneqInte (gamma, kT, dt,
			   dynamic_cast<Perturbation *> (&pert),
			   dynamic_cast<Force *> (&dw),
			   seed+2);

  double inteSigma = noneqInte.getSigma();
  double quenchSigma = quenchInte.getSigma();
  std::cout << "inte sigma is: " << inteSigma << std::endl;
  std::cout << "quench sigma is: " << quenchSigma << std::endl;

  Dofs xx;
  xx.xx[0] = 0.;
  xx.vv[0] = 0.;

  unsigned numNoneqCheck = int((noneqTime + 0.5 * noneqCheckFeq) / noneqCheckFeq) + 1;
  Distribution_1d dist;
  Distribution_1d distQuench;
  vector<double > checkTimes(numNoneqCheck);
  vector<string > xFileNames(numNoneqCheck);
  vector<string > xvFileNames(numNoneqCheck);
  vector<string > xQuenchFileNames(numNoneqCheck);
  vector<string > xvQuenchFileNames(numNoneqCheck);
  vector<string > diffxFileNames(numNoneqCheck);
  vector<string > diffxvFileNames(numNoneqCheck);
  dist      .reinit (x0, x1, nx, v0, v1, nv);
  distQuench.reinit (x0, x1, nx, v0, v1, nv);
  for (unsigned ii = 0; ii < numNoneqCheck; ++ii){
    checkTimes[ii] = ii * noneqCheckFeq;
    int timeI = int(checkTimes[ii] + 0.005);
    int timeF = int(100 * (checkTimes[ii] + 0.005 - timeI));
    char name[2048];
    sprintf (name, "distrib.resp.x.%05d.%02d.out", timeI, timeF);
    xFileNames[ii] = string(name);
    sprintf (name, "distrib.resp.vx.%05d.%02d.out", timeI, timeF);
    xvFileNames[ii] = string(name);
    sprintf (name, "distrib.resp.quench.x.%05d.%02d.out", timeI, timeF);
    xQuenchFileNames[ii] = string(name);
    sprintf (name, "distrib.resp.quench.vx.%05d.%02d.out", timeI, timeF);
    xvQuenchFileNames[ii] = string(name);
    sprintf (name, "indicator.resp.x.%05d.%02d.out", timeI, timeF);
    diffxFileNames[ii] = string(name);
    sprintf (name, "indicator.resp.vx.%05d.%02d.out", timeI, timeF);
    diffxvFileNames[ii] = string(name);
  }
  
  int count = 0;
  int countBranch = 0;
  double time = 0;

  NoneqResponseInfo resInfo;
  resInfo.reinit (x0, x1, nx, v0, v1, nv, dt, noneqTime, noneqCheckFeq, quenchTime, pert);

  sw_total.start();
  if (!vm.count("load-corr")){
    for (double ii = 0.; ii < nst+0.1; ii += 1.){
      // sw_eq.start();
      inte.step (xx, 0.);
      // sw_eq.stop();
      count ++;
      countBranch ++;
      time += dt;
      if (int(nstprint) == count){
	count = 0;
	printf ("%f %f %f\n", time, xx.xx[0], xx.vv[0]);
	fflush (stdout);
      }
      if (countBranch == int(branchNumFeq)){
	countBranch = 0;
	Dofs branchXX (xx);
	Dofs branchXX_old (xx);
	resInfo.newTraj ();
	// at time 0, quench any way
	{
	  // quenching
	  resInfo.newQuenchTraj ();
	  Dofs quenchXX (branchXX);
	  Dofs quenchXX_old (branchXX);
	  for (double ttQuench = 0.; ttQuench < quenchTime-0.5*dt; ttQuench += dt){
	    quenchXX_old = quenchXX;
	    // sw_quench.start();
	    quenchInte.step(quenchXX, 0.);
	    // sw_quench.stop();
	    Dofs quenchDw = quenchInte.getDw();
	    resInfo.depositQuenchTraj (quenchXX_old, quenchXX, quenchSigma, quenchDw);
	  }
	}      
	int countNoneqCheck = 0;
	// branching
	for (double ttNoneq = 0.; ttNoneq < noneqTime-0.5*dt; ttNoneq += dt){
	  // printf ("%f %d %d\n", ttNoneq, countNoneqCheck, noneqCheckNumFeq);
	  branchXX_old = branchXX;
	  // sw_noneq.start();
	  noneqInte.step (branchXX, ttNoneq);
	  // sw_noneq.stop();
	  Dofs dw = noneqInte.getDw ();
	  resInfo.depositMainTraj (branchXX_old, branchXX, inteSigma, dw);
	  countNoneqCheck ++;
	  if (countNoneqCheck == int(noneqCheckNumFeq)){
	    // printf ("check non eq at time %f \n", ttNoneq);
	    resInfo.newQuenchTraj ();
	    countNoneqCheck = 0;
	    // quenching
	    Dofs quenchXX (branchXX);
	    Dofs quenchXX_old (branchXX);
	    for (double ttQuench = 0.; ttQuench < quenchTime-0.5*dt; ttQuench += dt){
	      quenchXX_old = quenchXX;
	      // sw_quench.start();
	      quenchInte.step(quenchXX, ttNoneq);
	      // sw_quench.stop();
	      Dofs quenchDw = quenchInte.getDw();
	      resInfo.depositQuenchTraj (quenchXX_old, quenchXX, quenchSigma, quenchDw);
	    }
	  }
	}
      }
    }
    resInfo.average ();
    if (vm.count("save-corr")){
      resInfo.save (sfile);
    }
  }
  else {
    resInfo.load (lfile);
  }  


  for (unsigned ii = 0; ii < numNoneqCheck; ++ii){
    double nowTime = ii * noneqCheckFeq;
    resInfo.calculate (nowTime, pert1, dist, distQuench, order);
    dist.print_x  (xFileNames[ii]);
    dist.print_xv (xvFileNames[ii]);
    distQuench.print_x  (xQuenchFileNames[ii]);
    distQuench.print_xv (xvQuenchFileNames[ii]);
    
    dist.substract (distQuench);
    dist.print_x  (diffxFileNames[ii]);
    dist.print_xv (diffxvFileNames[ii]);
  }
  sw_total.stop();

  cout << "time static: user, real, sys" << endl;
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
  cout << "total:    "
       << sw_total.user() << " "
       << sw_total.real() << " "
       << sw_total.system() << endl;
  
  return 0;
}


