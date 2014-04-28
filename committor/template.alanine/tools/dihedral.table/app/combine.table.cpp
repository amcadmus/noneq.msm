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
  std::string ifile0, ifile1, ofile;
  double a0, a1;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("input-0", po::value<string > (&ifile0)->default_value ("table_d0.xvg"), "the first input table")
      ("input-1", po::value<string > (&ifile1)->default_value ("j.xvg"), "the second input table")
      ("scale-0",  po::value<double > (&a0)->default_value (1.),   "scale of first input table")
      ("scale-1",  po::value<double > (&a1)->default_value (2.),   "scale of second input table")
      ("output,o", po::value<string > (&ofile)->default_value ("out.xvg"), "the output table");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  vector<double > xx0, vv0, dd0;
  vector<double > xx1, vv1, dd1;
  
  ifstream fpname0 (ifile0.c_str());
  if (!fpname0){
    cerr << "cannot open file " << ifile0 << endl;
    return 1;
  }
  char nameline [MaxLineLength];
  while (fpname0.getline(nameline, MaxLineLength)){
    if (nameline[0] == '#' || nameline[0] == '@') continue;
    vector<string > words;
    StringOperation::split (string(nameline), words);
    if (words.size() < 3) {
      cerr << "wrong file format of " << ifile0 << endl;
      exit (1);
    }
    xx0.push_back (atof(words[0].c_str()));
    vv0.push_back (atof(words[1].c_str()));
    dd0.push_back (atof(words[2].c_str()));
  }

  ifstream fpname1 (ifile1.c_str());
  if (!fpname1){
    cerr << "cannot open file " << ifile1 << endl;
    return 1;
  }
  while (fpname1.getline(nameline, MaxLineLength)){
    if (nameline[0] == '#' || nameline[0] == '@') continue;
    vector<string > words;
    StringOperation::split (string(nameline), words);
    if (words.size() < 3) {
      cerr << "wrong file format of " << ifile1 << endl;
      exit (1);
    }
    xx1.push_back (atof(words[0].c_str()));
    vv1.push_back (atof(words[1].c_str()));
    dd1.push_back (atof(words[2].c_str()));
  }

  if (xx0.size() != xx1.size()){
    cerr << "unconsistent files, exit" << endl;
    exit(1);
  }

  FILE * fout = fopen (ofile.c_str(), "w");
  if (fout == NULL){
    cout << "cannot open file " << ofile<< endl;
    exit (1);
  }
  fprintf (fout, "# combined table from %f x %s and %f x %s \n",
	   a0, ifile0.c_str(),
	   a1, ifile1.c_str());
  for (unsigned ii = 0; ii < xx0.size(); ++ii){
    fprintf (fout, "%f %e %e\n",
	     xx0[ii],
	     vv0[ii] * a0 + vv1[ii] * a1,
	     dd0[ii] * a0 + dd1[ii] * a1);
  }
  
  fclose (fout);
  
  return 0;
}
