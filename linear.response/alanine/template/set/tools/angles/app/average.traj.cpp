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
  unsigned numBlock = 20;
  double refh;
  float time_prec = .01;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("refh,r",  po::value<double > (&refh)->default_value(10), "size of bin (deg.)")
      ("output,o", po::value<std::string > (&ofile)->default_value ("avg.angle"), "the output of average angle")
      ("input,f",  po::value<std::string > (&ifile)->default_value ("angle.name"), "the file of file names");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  ifstream fpname (ifile.c_str());
  if (!fpname){
    std::cerr << "cannot open file " << ifile << std::endl;
    return 1;
  }
  char nameline [MaxLineLength];
  int nbin = 360 / refh;
  float time;
  ValueType phi, psi;
  vector<float > times;
  vector<Distribution_1d>  (dists);
  unsigned countFile = 0;

  while (fpname.getline(nameline, MaxLineLength)){
    if (nameline[0] == '#') continue;
    FILE *fp = fopen (nameline, "r");
    cout << "reading file " << nameline << endl;
    if (fp == NULL){
      std::cerr << "cannot open file " << nameline << std::endl;
      return 1;
    }
    countFile ++;
    if (countFile == 1){
      while (myread(fp, time, phi, psi)){
	times.push_back (time);
	Distribution_1d tmpdist (-180, 180, nbin, -180, 180, nbin);
	tmpdist.deposite (psi, phi);
	dists.push_back (tmpdist);
      }
    }
    else {
      unsigned countFrame = 0;
      while (myread(fp, time, phi, psi)){
	if (countFrame >= times.size()){
	  cerr << "inconsistent frames" << endl;
	  return 1;
	}
	dists[countFrame].deposite (psi, phi);
	countFrame ++;
      }
    }
    fclose (fp);
  }

  for (unsigned ii = 0; ii < times.size(); ++ii){
    int inttime = int (times[ii] + time_prec * 0.5);
    int dectime = int ((times[ii] + time_prec * 0.5 - inttime) * 100);
    char frameFileName [MaxLineLength];
    sprintf (frameFileName, "%s.%05d.%02d", ofile.c_str(), inttime, dectime);
    dists[ii].average ();
    dists[ii].print_xv (frameFileName);
  }
  
  return 0;
}
