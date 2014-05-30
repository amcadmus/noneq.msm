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

#define MaxLineLength 65536

namespace po = boost::program_options;
using namespace std;

int main(int argc, char * argv[])
{
  std::string ifile, ofile;
  double bin;
  int nuse;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("bin,r",  po::value<double > (&bin)->default_value (.1),   "the bin size")
      ("num-use,n",  po::value<int > (&nuse)->default_value (0),   "the number of base used, 0 is all")
      ("input,f", po::value<std::string > (&ifile)->default_value ("base.k"), "the input base file")
      ("output,o", po::value<std::string > (&ofile)->default_value ("table_dih.xvg"), "the output of dihedral potential");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  ifstream fpname (ifile.c_str());
  if (!fpname){
    std::cerr << "\n cannot open file " << ifile << std::endl;
    return 1;
  }
  unsigned ndata = (360. + 0.01 * bin) / bin;
  bin = 360. / double(ndata);
  ndata ++;
  cout << "# print table with bin size: " << bin << endl;
  cout << "# ndata: " << ndata << endl;
  vector<double > xx (ndata);
  vector<double > vv (ndata, 0.);
  for (unsigned ii = 0; ii < ndata; ++ii){
    xx[ii] = -180. + ii * bin;
  }

  int count = 0;
  
  char valueline [MaxLineLength];
  while (fpname.getline(valueline, MaxLineLength)){
    if (valueline[0] == '#') continue;
    ++count;
    if (nuse != 0 && count > nuse) break;
    cout << "# process " << count << "th line: " << valueline << endl;
    vector<string > words;
    StringOperation::split (string(valueline), words);
    if (words.size() != 3){
      cerr << "wrong format, exit " << endl;
      return 1;
    }
    double mode = atof (words[1].c_str());
    double magn = atof (words[2].c_str());
    if (words[0] == string("sin")){
      for (unsigned ii = 0; ii < ndata; ++ii){
	vv[ii] += magn * sin (mode * xx[ii] * M_PI / 180.);
      }
    }
    else if (words[0] == string("cos")){
      for (unsigned ii = 0; ii < ndata; ++ii){
	vv[ii] += magn * cos (mode * xx[ii] * M_PI / 180.);
      }
    }
  }
    
  FILE * fout = fopen (ofile.c_str(), "w");
  if (fout == NULL){
    cout << "cannot open file " << ofile<< endl;
    exit (1);
  }
  for (unsigned ii = 0; ii < ndata; ++ii){
    fprintf (fout, "%f %e\n", xx[ii], vv[ii]);
  }
  fclose (fout);
  
  return 0;
}
