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
#include "StringSplit.h"
#include "Distribution.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;

#define MaxLineLength 2048

bool myread (FILE * fp,
	     float & time,
	     ValueType & angle)
{
  char line [MaxLineLength];
  
  while (1){
    char* rv = fgets (line, MaxLineLength, fp);
    if (rv == NULL) return false;
    if (line[0] == '#' || line[0] == '@') {
      continue;
    }
    else {
      vector<string> words;
      StringOperation::split (string(line), words);
      if (words.size() < 2) {
	cerr << "error of reading a line!! exit" << endl;
	exit(1);
      }
      time = atof(words[0].c_str());
      angle = atof(words[1].c_str());
      break;
    }
  }
  return true;
}

bool myreadGxs (FILE * fp,
		float & time,
		ValueType & gxs)
{
  size_t rv;
  rv = fscanf (fp, "%f", &time);
  if (rv != 1){
    // cerr << "error read time or reach EOF" << endl;
    return false;
  }
  rv = fscanf (fp, "%lf", &gxs);  
  if (rv != 1){
    cerr << "error read gxs " << endl;
    exit(1);
  }
  return true;
}



int main(int argc, char * argv[])
{
  std::string idfile, iffile, ofile, gfile;
  double tol, ds;
  float time_prec = .01;

  po::options_description desc ("Allow options");
  desc.add_options()
    ("help,h", "print this message")
    ("tol,t", po::value<double > (&tol)->default_value (30.), "the default value of tolrence of trans")
    ("ds,s", po::value<double > (&ds)->default_value (1.), "value of ds")
    ("output,o", po::value<std::string > (&ofile)->default_value ("prob.resp.out"), "the output of probability distribution of the trans conformation")
    ("gxs-file,g",  po::value<std::string > (&gfile)->default_value ("gxs.name"), "the file of file names: gxs")
    ("input-dir,d",  po::value<std::string > (&idfile)->default_value ("success.dir.name"), "the file of successful dirs")
    ("input-file,f",  po::value<std::string > (&iffile)->default_value ("angaver.xvg"), "the file name");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  ifstream fpname (idfile.c_str());
  if (!fpname){
    std::cerr << "cannot open file " << idfile << std::endl;
    return 1;
  }
  ifstream fpname1 (gfile.c_str());
  if (!fpname1){
    std::cerr << "cannot open file " << gfile << std::endl;
    return 1;
  }
  char nameline [MaxLineLength];
  char nameline2 [MaxLineLength];
  float time, time1;
  ValueType angle, gxs;
  vector<float > times;
  vector<double > counts;
  vector<double > dcounts;
  unsigned countFile = 0;

  while (fpname.getline(nameline, MaxLineLength)&&
	 fpname1.getline(nameline2, MaxLineLength)){
    if (nameline[0] == '#') continue;
    if (nameline[0] == '@') continue;
    if (nameline2[0] == '#') continue;
    if (nameline2[0] == '@') continue;
    char nameline1 [MaxLineLength];
    sprintf (nameline1, "%s/%s", nameline, iffile.c_str());
    FILE *fp = fopen (nameline1, "r");
    FILE *fp1 = fopen (nameline2, "r");
    cout << "reading file " << nameline1
	 << " and " << nameline2 << endl;
    if (fp == NULL){
      std::cerr << "cannot open file " << nameline1 << std::endl;
      return 1;
    }
    if (fp1 == NULL){
      std::cerr << "cannot open file " << nameline2 << std::endl;
      return 1;
    }
    countFile ++;
    if (countFile == 1){
      while (myread(fp, time, angle)){
	if (!myreadGxs(fp1, time1, gxs) || fabs(time - time1) > time_prec){
	  cerr << "inconsistent files" << endl;
	  return (1);
	}
	times.push_back (time);
	if (angle < -180 + tol || angle > 180 - tol) {
	  counts.push_back (1.);
	  dcounts.push_back (1. * ds * gxs);
	}
	else {
	  counts.push_back (0.);
	  dcounts.push_back (0. * ds * gxs);
	}
      }
    }
    else {
      unsigned countFrame = 0;
      while (myread(fp, time, angle)){
	if (countFrame >= times.size()){
	  cerr << "inconsistent frames" << endl;
	  return 1;
	}
	if (!myreadGxs(fp1, time1, gxs) || fabs(time - time1) > time_prec){
	  cerr << "inconsistent files" << endl;
	  return (1);
	}
	if (angle < -180 + tol || angle > 180 - tol) {
	  counts[countFrame] += (1.);
	  dcounts[countFrame] += (1.) * ds * gxs;
	}
	countFrame ++;
      }
    }
    fclose (fp);
    fclose (fp1);
  }

  FILE * fp = fopen (ofile.c_str(), "w");
  fprintf (fp, "# print the result of ds: %f\n", ds);
  for (unsigned ii = 0; ii < times.size(); ++ii){
    fprintf (fp, "%f   %f %f %f\n",
	     times[ii],
	     (counts[ii] + dcounts[ii]) / double(countFile),
	     counts[ii] / double(countFile),
	     dcounts[ii] / double(countFile)
	     );
  }
  fclose (fp);
  
  return 0;
}
