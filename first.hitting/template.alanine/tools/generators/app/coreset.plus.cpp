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

#include "SetTable.h"

#define MaxLineLength 2048

namespace po = boost::program_options;
using namespace std;

int main(int argc, char * argv[])
{
  std::string i1file, i2file, ofile;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("input-coreset-1", po::value<std::string > (&i1file)->default_value ("coreset.1.dat"), "file of coresets")
      ("input-coreset-2", po::value<std::string > (&i2file)->default_value ("coreset.2.dat"), "file of coresets")
      ("output,o", po::value<std::string > (&ofile)->default_value ("coreset.dat"), "file of coresets");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  AngleSetTable2D table1, table2;
  table1.reinit (i1file);
  table2.reinit (i2file);

  table1.plus (table2);

  table1.print (ofile);
  
  return 0;
}

