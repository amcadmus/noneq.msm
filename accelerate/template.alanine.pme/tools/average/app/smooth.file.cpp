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
      ("output,o", po::value<std::string > (&ofile)->default_value ("meta.flux.smooth.out"), "smoothed file")
      ("input,f",  po::value<std::string > (&ifile)->default_value ("meta.flux.out"), "the file to smoothify");
      
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

  
  FILE * fout = fopen (ofile.c_str(), "w");
  for (unsigned ii = 0; ii < values.size(); ++ii){
    if (ii == 0 || ii == values.size()-1){
      for (unsigned jj = 0; jj < values[ii].size(); ++jj){
	fprintf (fout, "%f  ", values[ii][jj]);
      }
      fprintf (fout, "\n");
    }
    else {
      fprintf (fout, "%f  ", values[ii][0]);
      for (unsigned jj = 1; jj < values[ii].size(); ++jj){
	fprintf (fout, "%f  ", 0.25 * (values[ii][jj] * 2. + values[ii-1][jj] + values[ii+1][jj]));
      }
      fprintf (fout, "\n");
    }
  }
  
  fclose (fout);
  
  return 0;
}
