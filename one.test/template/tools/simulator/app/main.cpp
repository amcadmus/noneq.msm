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
#include <boost/program_options.hpp>

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
  unsigned nst;
  unsigned nstprint = 1;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("dt,d",po::value<double > (&dt)->default_value(0.01), "time step [ps]")
      ("nst,n",po::value<unsigned > (&nst)->default_value(1000), "number of time step")
      ("nstprint,p",po::value<unsigned > (&nstprint)->default_value(10), "print frequency")
      ("gamma,g", po::value<double > (&gamma)->default_value(1.), "gamma [ps^-1]")
      ("temperture,t",po::value<double > (&T)->default_value(300.), "temerature [K]")
      ("double-well-k,k",po::value<double > (&kk)->default_value(8.0), "k parameter of the double well potential [kJ/(mol nm^4)]")
      ("double-well-a,a",po::value<double > (&aa)->default_value(1.0), "a parameter of the double well potential [nm]")
      ("pert-strength,s",po::value<double > (&pertSt)->default_value(1.0), "perturbation strength [kJ/(mol nm)]");
      
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
  std::cout << "# k: " << kk << " [kJ/mol]" << std::endl;
  std::cout << "# a: " << aa << " [nm]" << std::endl;
  std::cout << "# barrier: " << 0.5 * kk * aa*aa*aa*aa << " [kJ/mol]" << std::endl;
  std::cout << "# gamma: " << gamma << " [ps^-1]" << std::endl;
  std::cout << "# pert st: " << pertSt << " [kJ/(mol nm)]" << std::endl;
  std::cout << "###################################################" << std::endl;  

  DoubleWell dw (kk, aa);
  PertConstTilt pert (pertSt);
  
  EulerMaruyama inte (gamma, kT, dt,
		      NULL,
		      dynamic_cast<DoubleWell *> (&dw));
  Dofs xx;
  xx.xx[0] = aa;
  xx.vv[0] = 0.;
  double time = 0.;
  
  printf ("%f %f %f\n", time, xx.xx[0], xx.vv[0]);
  for (unsigned ii = 0; ii < nst; ++ii){
    inte.step (xx, time);
    time += dt;
    if ((ii+1) % nstprint == 0){
      printf ("%f %f %f\n", time, xx.xx[0], xx.vv[0]);
    }
  }

  return 0;
}
