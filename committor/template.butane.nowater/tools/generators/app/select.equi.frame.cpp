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
#include "StringSplit.h"

#define MaxLineLength 2048

namespace po = boost::program_options;
using namespace std;

int main(int argc, char * argv[])
{
  std::string icfile, iafile, ofile;
  double begin;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("begin,b",  po::value<double > (&begin)->default_value (0.), "start of the traj")
      ("input-coreset,c", po::value<string > (&icfile)->default_value ("coreset.dat"), "input selection map")
      ("input-angle,a", po::value<string > (&iafile)->default_value ("angaver.xvg"), "input angle")
      ("output,o", po::value<string > (&ofile)->default_value ("equi.frame"), "output equi frames");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  AngleSetTable1D stable ;
  stable.reinit (icfile);

  ifstream angle (iafile.c_str());
  if (!angle){
    cerr << "cannot open angle file " << angle << endl;
    return 1;
  }

  char valueline_ang [MaxLineLength];
  double prev_time = 0;
  int countTraj = 0;
  FILE * fp = fopen (ofile.c_str(), "w");

  while (angle.getline(valueline_ang, MaxLineLength)){
    if (valueline_ang[0] == '#' || valueline_ang[0] == '@'){
      continue;
    }
    vector<string > words;
    StringOperation::split (string(valueline_ang), words);
    if (words.size() < 2) {
      cerr << "wrong file format of " << iafile << endl;
      exit (1);
    }
    double time, anglev;
    time   = (atof(words[0].c_str()));
    anglev = (atof(words[1].c_str()));

    if (time < begin - 1e-3){
      continue;
    }
    else {
      if (countTraj %1000 == 0){
	printf ("read traj at %f     \r", time);
	fflush (stdout);
      }
    }
    countTraj ++;
    if (countTraj == 1){
      prev_time = time;
      continue;
    }

    if (stable.calIndicate (anglev) == 2){
      fprintf (fp, "%08d\t%f\n", countTraj, time - 0.5 * (time - prev_time));
    }
    
  }
  
  if (fclose(fp) != 0) {
    std::cerr << "ERROR: " << " calling fclose()"
                << " at " << __FILE__ << ":" << __LINE__
                << std::endl << std::flush;
      exit(1);
  }
  
      
  return 0;
}
