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
#include <boost/program_options.hpp>
#include "StringSplit.h"
#define MaxLineLength 2048

namespace po = boost::program_options;
using namespace std;

int main(int argc, char * argv[])
{
  std::string ifiles, slists, ofile;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("input-list", po::value<string > (&ifiles)->default_value ("table_d0.xvg"), "input tables")
      ("scale-list", po::value<string > (&slists)->default_value ("1.0"),   "scale of the tables")
      ("output,o", po::value<string > (&ofile)->default_value ("out.xvg"), "the output table");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }
  
  // analyze input
  vector<string > fileNames;
  vector<string > scaleNames;
  unsigned nFile = 0;
  StringOperation::split (ifiles, fileNames);
  StringOperation::split (slists, scaleNames);
  if (fileNames.size() != scaleNames.size()){
    cerr << "in consistent scale and file names size" << endl;
    return 1;
  }
  nFile = fileNames.size();
  vector<double > scales (nFile);
  for (unsigned ii = 0; ii < nFile; ++ii){
    scales[ii] = atof (scaleNames[ii].c_str());
  }
  cout << "# here are tables and there sizes" << endl;
  for (unsigned ii = 0; ii < nFile; ++ii){
    cout << "# " << fileNames[ii] << "  " << scales[ii] << endl;
  }
  
  vector<double > xx, vv, dd;

  for (unsigned ii = 0; ii < nFile; ++ii){
    ifstream fpname (fileNames[ii].c_str());
    if (!fpname){
      cerr << "cannot open file " << fileNames[ii] << endl;
      return 1;
    }
    char nameline [MaxLineLength];
    int countLine = 0;
    while (fpname.getline(nameline, MaxLineLength)){
      if (nameline[0] == '#' || nameline[0] == '@') continue;
      vector<string > words;
      StringOperation::split (string(nameline), words);
      if (words.size() < 3) {
	cerr << "wrong file format of " << fileNames[ii] << endl;
	exit (1);
      }
      if (ii == 0){
	xx.push_back (atof(words[0].c_str()));
	vv.push_back (scales[ii] * atof(words[1].c_str()));
	dd.push_back (scales[ii] * atof(words[2].c_str()));
      }
      else {
	if (countLine > int(xx.size()) - 1) {
	  cerr << "inconsistent number of lines of file " << fileNames[ii] << endl;
	  return 1;
	}
	if (xx[countLine] != atof(words[0].c_str())){
	  cerr << "inconsistent x value of tables " << fileNames[0] << " and " << fileNames[ii] << endl;
	  return 1;
	}
	vv[countLine] += (scales[ii] * atof(words[1].c_str()));
	dd[countLine] += (scales[ii] * atof(words[2].c_str()));	
      }
      countLine ++;
    }
  }

  FILE * fout = fopen (ofile.c_str(), "w");
  if (fout == NULL){
    cout << "cannot open file " << ofile<< endl;
    exit (1);
  }
  fprintf (fout, "# combined table from:\n");
  for (unsigned ii = 0; ii < nFile; ++ii){
    fprintf (fout, "# %f %s\n", scales[ii], fileNames[ii].c_str());
  }
  for (unsigned ii = 0; ii < xx.size(); ++ii){
    fprintf (fout, "%f %e %e\n", xx[ii], vv[ii], dd[ii]);
  }
  
  fclose (fout);
  
  return 0;
}
