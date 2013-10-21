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
  
  po::options_description desc ("Allow options");
  desc.add_options()
    ("help,h", "print this message")
    ("refh,r",  po::value<double > (&refh)->default_value(10), "size of bin (deg.)")
    ("output,o", po::value<std::string > (&ofile)->default_value ("avg.angle.out"), "the output of average angle")
    ("input,f",  po::value<std::string > (&ifile)->default_value ("angle.dat"), "the file of file names");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  int nbin = 360 / refh;
  float time;
  ValueType phi, psi;
  Distribution_1d  dist (-180, 180, nbin, -180, 180, nbin);

  FILE *fp = fopen (ifile.c_str(), "r");
  cout << "reading file " << ifile << endl;
  if (fp == NULL){
    std::cerr << "cannot open file " << ifile << std::endl;
    return 1;
  }
  while (myread(fp, time, phi, psi)){
    dist.deposite (psi, phi);
  }

  char frameFileName [MaxLineLength];
  sprintf (frameFileName, "%s", ofile.c_str());
  dist.average ();
  dist.print_xv (frameFileName);
  
  return 0;
}
