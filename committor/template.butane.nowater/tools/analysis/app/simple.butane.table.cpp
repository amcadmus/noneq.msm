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
#include "BlockAverage.h"
#include "SetTable.h"

#define MaxLineLength 2048

namespace po = boost::program_options;
using namespace std;

int main(int argc, char * argv[])
{
  std::string ofile;
  double width, hbin;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("width,w", po::value<double > (&width)->default_value (40), "width of the set")
      ("hbin,b", po::value<double > (&hbin)->default_value (5), "width of the set")
      ("output,o", po::value<std::string > (&ofile)->default_value ("coreset.dat"), "file of coresets");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  int nbin = (360. + 0.5 * hbin)/(double(hbin));

  FILE * fp = fopen (ofile.c_str(), "w");
  for (int ii = 0; ii < nbin; ++ii){
    double xx = -180 + 0.5 * hbin + ii * hbin;
    if (xx < -180 + 0.5 * width || xx > 180 - 0.5 * width){
      fprintf (fp, "%d\n", 1);
    }
    else if (xx > -60 - 0.5 * width && xx < -60 + 0.5 * width){
      fprintf (fp, "%d\n", -1);
    }
    else if (xx >  60 - 0.5 * width && xx <  60 + 0.5 * width){
      fprintf (fp, "%d\n", -1);
    }
    else {
      fprintf (fp, "%d\n", 0);
    }  
  }
  
  return 0;
}
