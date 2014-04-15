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
  std::string icfile;
  double angle;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("input-angle,a", po::value<double > (&angle)->default_value (0), "the angle in degree")
      ("input-coreset-file,f", po::value<std::string > (&icfile)->default_value ("coreset.dat"), "the file of core sets");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  AngleSetTable1D stable;
  stable.reinit (icfile);

  cout << stable.calIndicate (angle) << endl;

  return 0;
}

