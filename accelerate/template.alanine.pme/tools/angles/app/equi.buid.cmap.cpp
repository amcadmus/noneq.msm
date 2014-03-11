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

#include "Defines.h"
#include "Distribution.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;
#define MaxLineLength 2048

bool myread (FILE * fp,
	     float & time,
	     ValueType & phi,
	     ValueType & psi)
{
  size_t rv;
  rv = fread (&time, sizeof(float), 1, fp);
  if (rv != 1){
    // cerr << "error read time or reach EOF" << endl;
    return false;
  }
  rv = fread (&phi, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error read phi " << endl;
    exit(1);
  }
  rv = fread (&psi, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error read psi " << endl;
    exit(1);
  }
  return true;
}


int main(int argc, char * argv[])
{
  std::string ifile, ofile;
  double refh;
  int nbin;
  
  po::options_description desc ("Allow options");
  desc.add_options()
    ("help,h", "print this message")
    ("nbin,n",  po::value<int > (&nbin)->default_value(24), "number of bins")
    ("output,o", po::value<std::string > (&ofile)->default_value ("avg.angle.out"), "the output of average angle")
    ("input,f",  po::value<std::string > (&ifile)->default_value ("traj.name"), "the file of traj names");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  refh = 360. / double(nbin);
  float time;
  ValueType phi, psi;
  // Distribution_1d  dist (-180 - 0.0 * refh, 180 - 0.0 * refh, nbin,
  // 			 -180 - 0.0 * refh, 180 - 0.0 * refh, nbin);
  Distribution_1d  dist (-180 - 0.5 * refh, 180 - 0.5 * refh, nbin,
  			 -180 - 0.5 * refh, 180 - 0.5 * refh, nbin);

  ifstream fpname (ifile.c_str());
  if (!fpname){
    std::cerr << "cannot open file " << ifile << std::endl;
    return 1;
  }
  char nameline [MaxLineLength];
  while (fpname.getline(nameline, MaxLineLength)){
    if (nameline[0] == '#') continue;
    FILE *fp = fopen (nameline, "r");
    cout << "reading file " << nameline << endl;
    if (fp == NULL){
      std::cerr << "cannot open file " << ifile << std::endl;
      return 1;
    }
    while (myread(fp, time, phi, psi)){
      dist.deposite (psi, phi);
    }
    fclose (fp);
  }
  
  dist.average ();
    
  char frameFileName [MaxLineLength];
  sprintf (frameFileName, "%s", ofile.c_str());
  dist.print_xv (frameFileName);

  dist.print_x ("avg.phi.out");
  dist.print_v ("avg.psi.out");
  
  return 0;
}
