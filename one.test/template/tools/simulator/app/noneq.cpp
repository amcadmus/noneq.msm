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
  unsigned noneqCheckNumFeq;
  unsigned branchNumFeq;
      
  unsigned long seed;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("dt,d", po::value<double > (&dt)->default_value(0.01), "time step [ps]")
      ("nst,n", po::value<double > (&nst)->default_value(1000), "number of time step")
      ("nstprint,p", po::value<unsigned > (&nstprint)->default_value(10), "print frequency")
      ("gamma,g", po::value<double > (&gamma)->default_value(1.), "gamma [ps^-1]")
      ("temperture,t", po::value<double > (&T)->default_value(300.), "temerature [K]")
      ("double-well-k,k", po::value<double > (&kk)->default_value(8.0), "k parameter of the double well potential [kJ/(mol nm^4)]")
      ("double-well-a,a", po::value<double > (&aa)->default_value(1.0), "a parameter of the double well potential [nm]")
      ("branch-feq", po::value<double > (&branchFeq)->default_value(1.), "branch frequency [ps]")
      ("noneq-check-feq", po::value<double > (&noneqCheckFeq)->default_value(10.), "non-equilibrium branch check frequency [ps]")
      ("noneq-time", po::value<double > (&noneqTime)->default_value(200.), "non-equilibrium simulation time [ps]")
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
  noneqCheckNumFeq = (noneqCheckFeq+0.5*dt) / dt;
  
  std::cout << "###################################################" << std::endl;
  std::cout << "# T: " << T << " [K]" << std::endl;
  std::cout << "# kT: " << kT << " [kJ/mol]" << std::endl;
  std::cout << "# k: " << kk << " [kJ/mol]" << std::endl;
  std::cout << "# a: " << aa << " [nm]" << std::endl;
  std::cout << "# barrier: " << 0.5 * kk * aa*aa*aa*aa << " [kJ/mol]" << std::endl;
  std::cout << "# gamma: " << gamma << " [ps^-1]" << std::endl;
  std::cout << "# pert st: " << pertSt << " [kJ/(mol nm)]" << std::endl;
  std::cout << "# quench T: " << quenchTemperture << " [K]" << std::endl;
  std::cout << "# quench kT: " << quenchkT << " [kJ/mol]" << std::endl;
  std::cout << "# branch Feq: " << branchFeq << " [ps]" << std::endl;
  std::cout << "# branch every: " << branchNumFeq << " steps" << std::endl;
  std::cout << "# noneq check Feq: " << noneqCheckFeq << " [ps]" << std::endl;
  std::cout << "# noneq check every: " << noneqCheckNumFeq << " steps" << std::endl;
  std::cout << "# noneq time: " << noneqTime << " [ps]" << std::endl;
  std::cout << "###################################################" << std::endl;  

  DoubleWell dw (kk, aa);
  PertConstTilt pert (pertSt, warmTime);
  // PertConstTilt pert (0.);
  
  EulerMaruyama inte (gamma, kT, dt,
		      NULL,
		      dynamic_cast<Force *> (&dw));
  EulerMaruyama quenchInte (gamma, quenchkT, dt,
			    NULL,
			    dynamic_cast<Force *> (&dw));
  EulerMaruyama noneqInte (gamma, kT, dt,
			   dynamic_cast<Perturbation *> (&pert),
			   dynamic_cast<Force *> (&dw));
  
  Dofs xx;
  xx.xx[0] = aa;
  xx.vv[0] = 0.;

  unsigned numNoneqCheck = noneqTime / noneqCheckFeq;
  numNoneqCheck ++;
  vector<Distribution_1d > dists(numNoneqCheck);
  vector<Distribution_1d > distsQuench(numNoneqCheck);
  vector<double > checkTimes(numNoneqCheck);
  vector<string > xFileNames(numNoneqCheck);
  vector<string > xvFileNames(numNoneqCheck);
  vector<string > xQuenchFileNames(numNoneqCheck);
  vector<string > xvQuenchFileNames(numNoneqCheck);
  vector<string > diffxFileNames(numNoneqCheck);
  vector<string > diffxvFileNames(numNoneqCheck);
  for (unsigned ii = 0; ii < numNoneqCheck; ++ii){
    checkTimes[ii] = ii * noneqCheckFeq;
    int timeI = int(checkTimes[ii]);
    int timeF = int(100 * (checkTimes[ii] - timeI));
    dists[ii].reinit (-2, 2, 50, -8, 8, 50);
    distsQuench[ii].reinit (-2, 2, 50, -8, 8, 50);
    char name[2048];
    sprintf (name, "distrib.x.%05d.%02d.out", timeI, timeF);
    xFileNames[ii] = string(name);
    sprintf (name, "distrib.vx.%05d.%02d.out", timeI, timeF);
    xvFileNames[ii] = string(name);
    sprintf (name, "distrib.quench.x.%05d.%02d.out", timeI, timeF);
    xQuenchFileNames[ii] = string(name);
    sprintf (name, "distrib.quench.vx.%05d.%02d.out", timeI, timeF);
    xvQuenchFileNames[ii] = string(name);
    sprintf (name, "indicator.x.%05d.%02d.out", timeI, timeF);
    diffxFileNames[ii] = string(name);
    sprintf (name, "indicator.xv.%05d.%02d.out", timeI, timeF);
    diffxvFileNames[ii] = string(name);
  }
  
  int count = 0;
  int countBranch = 0;
  double ttNoneq;
  double ttQuench;
  unsigned iiCheck = 0;
  double time = 0;

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
      Dofs branchXX (xx);
      iiCheck = 0;
      // branching
      int countNoneqCheck = noneqCheckNumFeq;
      for (ttNoneq = 0.; ttNoneq <= noneqTime+0.5*dt; ttNoneq += dt){
	// printf ("%f %d %d\n", ttNoneq, countNoneqCheck, noneqCheckNumFeq);
 	if (countNoneqCheck == int(noneqCheckNumFeq)){
	  // printf ("check non eq at time %f \n", ttNoneq);
	  countNoneqCheck = 0;
	  dists[iiCheck].deposite (branchXX);
	  // quenching
	  Dofs quenchXX (branchXX);
	  for (ttQuench = 0.; ttQuench <= quenchTime+0.5*dt; ttQuench += dt){
	    quenchInte.step(quenchXX, ttNoneq);
	  }
	  distsQuench[iiCheck].deposite (quenchXX);
	  iiCheck ++;
	}
	noneqInte.step (branchXX, ttNoneq);
	countNoneqCheck ++;
      }
    }
  }
  
  // for (double ii = 0.; ii < nst; ii += 1.){
  //   inte.step (xx, 0.);
  //   time += dt;
  //   count ++;
  //   countBranch ++;
  //   if (int(nstprint) == count){
  //     count = 0;
  //     printf ("%f %f %f\n", time, xx.xx[0], xx.vv[0]);
  //     dist.deposite (xx);
  //   }
  //   if (int(branchFeq) == countBranch){
  //     // std::cout << "start branch" << std::endl;
  //     countBranch = 0;
  //     Dofs branchPosi (xx);
  //     for (double branchTime = 0; branchTime <= branchLength; branchTime += dt){
  // 	branchInte.step (branchPosi, 0.);
  // 	// printf ("%f %f %f\n", branchTime, branchPosi.xx[0], branchPosi.vv[0]);
  //     }
  //     branchDist.deposite (branchPosi);
  //   }  
  // }

  for (unsigned ii = 0; ii < numNoneqCheck; ++ii){
    dists[ii].average();
    distsQuench[ii].average();
    dists[ii].print_x  (xFileNames[ii]);
    dists[ii].print_xv (xvFileNames[ii]);
    distsQuench[ii].print_x  (xQuenchFileNames[ii]);
    distsQuench[ii].print_xv (xvQuenchFileNames[ii]);
    
    dists[ii].substract (distsQuench[ii]);
    dists[ii].print_x  (diffxFileNames[ii]);
    dists[ii].print_xv (diffxvFileNames[ii]);
  }

  return 0;
}
