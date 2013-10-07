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



int main(int argc, char * argv[])
{
  std::string idfile, iffile, ofile;
  double tol;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("tol,t", po::value<double > (&tol)->default_value (30.), "the default value of tolrence of trans")
      ("output,o", po::value<std::string > (&ofile)->default_value ("prob.out"), "the output of probability distribution of the trans conformation")
      ("input-dir,f",  po::value<std::string > (&idfile)->default_value ("success.dir.name"), "the file of successful dirs")
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
  char nameline [MaxLineLength];
  float time;
  ValueType angle;
  vector<float > times;
  vector<double > counts;
  unsigned countFile = 0;

  while (fpname.getline(nameline, MaxLineLength)){
    if (nameline[0] == '#') continue;
    if (nameline[0] == '@') continue;
    char nameline1 [MaxLineLength];
    sprintf (nameline1, "%s/%s", nameline, iffile.c_str());
    FILE *fp = fopen (nameline1, "r");
    cout << "reading file " << nameline1 << endl;
    if (fp == NULL){
      std::cerr << "cannot open file " << nameline1 << std::endl;
      return 1;
    }
    countFile ++;
    if (countFile == 1){
      while (myread(fp, time, angle)){
	times.push_back (time);
	if (angle < -180 + tol || angle > 180 - tol) {
	  counts.push_back (1.);
	}
	else {
	  counts.push_back (0.);
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
	if (angle < -180 + tol || angle > 180 - tol) {
	  counts[countFrame] += (1.);
	}
	countFrame ++;
      }
    }
    fclose (fp);
  }

  FILE * fp = fopen (ofile.c_str(), "w");
  for (unsigned ii = 0; ii < times.size(); ++ii){
    fprintf (fp, "%f %f\n", times[ii], counts[ii] / double(countFile));
  }
  fclose (fp);
  
  return 0;
}
