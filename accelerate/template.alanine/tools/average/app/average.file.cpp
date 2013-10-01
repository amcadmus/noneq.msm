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
  unsigned numBlock = 20;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("output,o", po::value<std::string > (&ofile)->default_value ("avg.h.bond.count.out"), "the output of count of h-bond")
      ("input,f",  po::value<std::string > (&ifile)->default_value ("h.count.name"), "the file of file names");
      
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
  char nameline [MaxLineLength];
  std::vector<std::vector<double > > values;
  vector<vector<vector<double > > > data;
  unsigned countFile = 0;
  
  while (fpname.getline(nameline, MaxLineLength)){
    if (nameline[0] == '#') continue;
    std::ifstream fptraj (nameline);
    if (!fptraj){
      std::cerr << "cannot open file " << nameline << std::endl;
      return 1;
    }
    countFile ++;
    if (countFile == 1){
      char line [MaxLineLength];
      while (fptraj.getline(line, MaxLineLength)){
	if (line[0] == '#') continue;
	std::string thisline (line);
	std::vector<std::string > words;
	StringOperation::split (thisline, words);
	std::vector<double > lineValues;
	vector<vector<double > > lineData;
	for (unsigned ii = 0; ii < words.size(); ++ii){
	  lineValues.push_back(atof(words[ii].c_str()));
	  lineData.push_back(vector<double >(1, lineValues.back()));
	}
	values.push_back (lineValues);
	data.push_back (lineData);
      }
    }
    else {
      unsigned lineCount = 0;
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
	if (lineCount >= values.size()){
	  std::cerr << "wrong number of lines of file "
		    << nameline << std::endl;
	  return 1;
	}
	if (lineValues.size() != values[lineCount].size()){
	  std::cerr << "broken line "
		    << lineCount
		    << " of file "
		    << nameline << std::endl;
	  return 1;
	}
	for (unsigned ii = 1; ii < lineValues.size(); ++ii){
	  values[lineCount][ii] += lineValues[ii];
	  data[lineCount][ii].push_back(lineValues[ii]);
	}
	lineCount ++;
      }
    }
  }

  BlockAverage ba;
  
  FILE * fout = fopen (ofile.c_str(), "w");
  for (unsigned ii = 0; ii < values.size(); ++ii){
    fprintf (fout, "%.3f  ", values[ii][0]);
    for (unsigned jj = 1; jj < values[ii].size(); ++jj){
      values[ii][jj] /= countFile;
      ba.processData (data[ii][jj], numBlock);
      fprintf (fout, "%f %f   ", values[ii][jj], ba.getAvgError());
    }
    fprintf (fout, "\n");
  }
  fclose (fout);
  
  return 0;
}
