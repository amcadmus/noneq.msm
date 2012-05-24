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
#include "Integrator.h"
#include "Distribution.h"
#include <boost/program_options.hpp>
#include "RandomGenerator.h"
#include "TimeCorrelation.h"

namespace po = boost::program_options;

int main(int argc, char * argv[])
{
  double gamma;
  double kT;
  double T;
  double dt;
  double pertSt;
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
  double corrTime;
  double corrStep;
  // unsigned noneqCheckNumFeq;
  unsigned branchNumFeq;
      
  unsigned long seed;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("dt,d", po::value<double > (&dt)->default_value(0.01), "time step [ps]")
      ("nst,n", po::value<double > (&nst)->default_value(1000.), "number of time step")
      ("print-feq,p", po::value<unsigned > (&nstprint)->default_value(10000), "print frequency")
      ("gamma,g", po::value<double > (&gamma)->default_value(1.), "gamma [ps^-1]")
      ("temperature,t", po::value<double > (&T)->default_value(300.), "temperature [K]")
      ("double-well-k,k", po::value<double > (&kk)->default_value(8.0), "k parameter of the double well potential [kJ/(mol nm^4)]")
      ("double-well-a,a", po::value<double > (&aa)->default_value(1.0), "a parameter of the double well potential [nm]")
      ("branch-feq", po::value<double > (&branchFeq)->default_value(1.), "branch frequency [ps]")
      ("noneq-check-feq", po::value<double > (&noneqCheckFeq)->default_value(10.), "non-equilibrium branch check frequency [ps]")
      ("noneq-time", po::value<double > (&noneqTime)->default_value(200.), "non-equilibrium simulation time [ps]")
      ("corr-time", po::value<double > (&corrTime)->default_value(10.), "correlation cal time [ps]")
      ("corr-step", po::value<double > (&corrStep)->default_value(.01), "step of correlation [ps]")
      ("quench-temperature", po::value<double > (&quenchTemperture)->default_value(150.), "quench temperature [K]")
      ("quench-time", po::value<double > (&quenchTime)->default_value(1.), "quench time [ps]")
      ("warm-time", po::value<double > (&warmTime)->default_value(100.), "warm up time of perturbation [ps]")
      ("seed",po::value<unsigned long > (&seed)->default_value(1), "random seed")      
      ("pert-strength,s",po::value<double > (&pertSt)->default_value(1.0), "perturbation strength [kJ/(mol nm)]");
      
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
  // noneqCheckNumFeq = (noneqCheckFeq+0.5*dt) / dt;
  
  std::cout << "###################################################" << std::endl;
  std::cout << "# T: " << T << " [K]" << std::endl;
  std::cout << "# kT: " << kT << " [kJ/mol]" << std::endl;
  std::cout << "# k: " << kk << " [kJ/(mol nm^4)]" << std::endl;
  std::cout << "# a: " << aa << " [nm]" << std::endl;
  std::cout << "# barrier: " << 0.5 * kk * aa*aa*aa*aa << " [kJ/mol]" << std::endl;
  std::cout << "# gamma: " << gamma << " [ps^-1]" << std::endl;
  std::cout << "# pert st: " << pertSt << " [kJ/(mol nm)]" << std::endl;
  std::cout << "# warm time: " << warmTime << " [ps]" << std::endl;
  std::cout << "# quench T: " << quenchTemperture << " [K]" << std::endl;
  std::cout << "# quench kT: " << quenchkT << " [kJ/mol]" << std::endl;
  std::cout << "# corr time: " << corrTime << " [ps]" << std::endl;
  std::cout << "# corr step: " << corrStep << " [ps]" << std::endl;
  // std::cout << "# branch Feq: " << branchFeq << " [ps]" << std::endl;
  // std::cout << "# branch every: " << branchNumFeq << " steps" << std::endl;
  std::cout << "# noneq check Feq: " << noneqCheckFeq << " [ps]" << std::endl;
  // std::cout << "# noneq check every: " << noneqCheckNumFeq << " steps" << std::endl;
  std::cout << "# noneq time: " << noneqTime << " [ps]" << std::endl;
  std::cout << "###################################################" << std::endl;  

  DoubleWell dw (kk, aa);
  PertConstTilt pert (pertSt, warmTime);
  DissipativeFlux flux (dynamic_cast<Perturbation *> (&pert),
			dynamic_cast<Force *> (&dw));
  // PertConstTilt pert (0.);
  
  EulerMaruyama inte (gamma, kT, dt,
		      NULL,
		      dynamic_cast<Force *> (&dw));
  EulerMaruyama quenchInte (gamma, quenchkT, dt,
			    NULL,
			    dynamic_cast<Force *> (&dw));
  
  Dofs xx;
  xx.xx[0] = aa;
  xx.vv[0] = 0.;

  TimeCorrelation tc;
  tc.reinit (corrStep, corrTime,
	     -2, 2, 50, -8, 8, 50,
	     &inte,
	     (quenchTime + .5*dt) / dt,
	     &quenchInte,
	     &flux);
  tc.calCorr (xx, nst*dt);
  // tc.print();
  
  int count = 0;
  int countBranch = 0;
  double time = 0;
  Distribution_1d dist       (-2, 2, 50, -8, 8, 50);
  Distribution_1d distQuench (-2, 2, 50, -8, 8, 50);

  for (double ii = 0.; ii < nst+0.1; ii += 1.){
    inte.step (xx, 0.);
    count ++;
    countBranch ++;
    time += dt;
    if (int(nstprint) == count){
      count = 0;
      printf ("%f %f %f\n", time, xx.xx[0], xx.vv[0]);
    }
    if (countBranch == int(branchNumFeq)){
      countBranch = 0;
      dist.deposite (xx);
      Dofs quenchXX (xx);
      // branching
      for (double ttQuench = 0.; ttQuench <= quenchTime+0.5*dt; ttQuench += dt){
	quenchInte.step(quenchXX, 0.);
      }
      distQuench.deposite (quenchXX);
    }
  }
  dist.average();
  distQuench.average();
  dist.substract (distQuench);
  
  vector<vector<vector<double > > > timeNew;
  tc.calIndicator (dist.values,
		   noneqTime, noneqCheckFeq, 1./kT, pert, timeNew);
  for (double tt = 0; tt < noneqTime + 0.5 * noneqCheckFeq; tt += noneqCheckFeq){
    unsigned ii = (tt + 0.5 * noneqCheckFeq) / noneqCheckFeq;
    dist.values = timeNew[ii];
    int timeI = int(tt);
    int timeF = int(100 * (tt - timeI));
    char name[2048];
    sprintf (name, "indicator.linear.x.%05d.%02d.out", timeI, timeF);
    dist.print_x (name);
    sprintf (name, "indicator.linear.vx.%05d.%02d.out", timeI, timeF);
    dist.print_xv (name);
  }
  
  return 0;
}

