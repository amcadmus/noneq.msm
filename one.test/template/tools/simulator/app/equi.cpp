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

namespace po = boost::program_options;

int main(int argc, char * argv[])
{
  double gamma;
  double kT;
  double T;
  double dt;
  // double pertSt;
  double aa;
  double kk;
  double nst;
  unsigned nstprint = 1;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("dt,d",po::value<double > (&dt)->default_value(0.01), "time step [ps]")
      ("nst,n",po::value<double > (&nst)->default_value(1000), "number of time step")
      ("nstprint,p",po::value<unsigned > (&nstprint)->default_value(10), "print frequency")
      ("gamma,g", po::value<double > (&gamma)->default_value(1.), "gamma [ps^-1]")
      ("temperature,t",po::value<double > (&T)->default_value(300.), "temperature [K]")
      ("double-well-k,k",po::value<double > (&kk)->default_value(8.0), "k parameter of the double well potential [kJ/(mol nm^4)]")
      ("double-well-a,a",po::value<double > (&aa)->default_value(1.0), "a parameter of the double well potential [nm]");
      // ("pert-strength,s",po::value<double > (&pertSt)->default_value(1.0), "perturbation strength [kJ/(mol nm)]");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }
  kT = T * 1.38 * 6.02 * 1e-3;
  
  std::cout << "###################################################" << std::endl;
  std::cout << "# T: " << T << " [K]" << std::endl;
  std::cout << "# kT: " << kT << " [kJ/mol]" << std::endl;
  std::cout << "# k: " << kk << " [kJ/(mol nm^4)]" << std::endl;
  std::cout << "# a: " << aa << " [nm]" << std::endl;
  std::cout << "# barrier: " << 0.5 * kk * aa*aa*aa*aa << " [kJ/mol]" << std::endl;
  std::cout << "# gamma: " << gamma << " [ps^-1]" << std::endl;
  // std::cout << "# pert st: " << pertSt << " [kJ/(mol nm)]" << std::endl;
  std::cout << "###################################################" << std::endl;  

  DoubleWell dw (kk, aa);
  // PertConstTilt pert (pertSt);
  
  EulerMaruyama inte (gamma, kT, dt,
		      NULL,
		      dynamic_cast<DoubleWell *> (&dw));
  Dofs xx;
  xx.xx[0] = aa;
  xx.vv[0] = 0.;
  double time = 0.;

  Distribution_1d dist(-2, 2, 50, -8, 8, 50);
  
  printf ("%f %f %f\n", time, xx.xx[0], xx.vv[0]);
  int count = 0;
  for (double ii = 0.; ii < nst; ii += 1.){
    inte.step (xx, time);
    time += dt;
    count ++;
    if (int(nstprint) <= count){
      printf ("%f %f %f\n", time, xx.xx[0], xx.vv[0]);
      dist.deposite (xx);
      count = 0;
    }
  }

  dist.average();

  FILE * fp = fopen ("distrib.vx.out", "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << std::endl;
    return 1;
  }
  dist.print_xv (fp);
  fclose (fp);
  fp = fopen ("distrib.x.out", "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << std::endl;
    return 1;
  }
  dist.print_x (fp);
  fclose (fp);

  return 0;
}
