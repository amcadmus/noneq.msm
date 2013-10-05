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
#define MaxLineLength 10240

#include "StringSplit.h"
#include "BlockAverage.h"
namespace po = boost::program_options;
using namespace std;

int main(int argc, char * argv[])
{
  std::string ifile, ofile;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("output,o", po::value<std::string > (&ofile)->default_value ("meta.flux.inte.out"), "the integrated file")
      ("input,f",  po::value<std::string > (&ifile)->default_value ("meta.flux.out"), "the file to integrate");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  std::ifstream fpname (ifile.c_str());
  if (!fpname){
    std::cerr << "cannot open file " << ifile << std::endl;
    return 1;
  }
  std::vector<std::vector<double > > values;
  std::ifstream fptraj (ifile.c_str());
  char line [MaxLineLength];

  while (fptraj.getline(line, MaxLineLength)){
    if (line[0] == '#') continue;
    std::string thisline (line);
    std::vector<std::string > words;
    StringOperation::split (thisline, words);
    std::vector<double > lineValues;
    for (unsigned ii = 0; ii < words.size(); ++ii){
      lineValues.push_back(atof(words[ii].c_str()));
    }
    values.push_back(lineValues);
  }
  fptraj.close();

  std::vector<std::vector<double > > sums (values.size());
  for (unsigned ii = 0; ii < values.size(); ++ii){
    for (unsigned jj = 0; jj < values[ii].size(); ++jj){
      sums[ii].resize (values[ii].size());
      if (ii == 0){
	sums[ii][jj] = values[ii][jj];
      }
      else{
	if (jj == 0){
	  sums[ii][jj] = values[ii][jj];
	}
	else{
	  sums[ii][jj] = sums[ii-1][jj] + values[ii][jj];
	}
      }
    }
  }
  
  
  FILE * fout = fopen (ofile.c_str(), "w");
  for (unsigned ii = 0; ii < sums.size(); ++ii){
    for (unsigned jj = 0; jj < sums[ii].size(); ++jj){
      if (jj == 0){
	fprintf (fout, "%f  ", sums[ii][jj]);
      }
      else {
	fprintf (fout, "%f  ", sums[ii][jj] * (values[1][0] - values[0][0]));
      }
    }
    fprintf (fout, "\n");
  }
  
  fclose (fout);
  
  return 0;
}
