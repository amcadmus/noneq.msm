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

namespace po = boost::program_options;
using namespace std;

int main(int argc, char * argv[])
{
  std::string ifile, ofile;
  double lower, upper, start;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("lower,l",  po::value<double > (&lower)->default_value (-100.),   "the lower bond of starting angle ")
      ("upper,u",  po::value<double > (&upper)->default_value ( 100.),   "the upper bond of starting angle ")
      ("start,s",  po::value<double > (&start)->default_value ( 100.),   "the start time of equilibrium ")
      ("input,f",  po::value<string > (&ifile)->default_value ("angaver.xvg"), "input file of angles")
      ("output,o", po::value<string > (&ofile)->default_value ("equi.frame"), "the output of equi frames");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  FILE * fin = fopen (ifile.c_str(), "r");
  if (fin == NULL){
    cout << "cannot open file " << ifile<< endl;
    exit (1);
  }
  FILE * fout = fopen (ofile.c_str(), "w");
  if (fout == NULL){
    cout << "cannot open file " << ofile<< endl;
    exit (1);
  }
  double time, time0, time1, angle;
  if (2 != fscanf (fin, "%lf %lf", &time0, &angle) ||
      2 != fscanf (fin, "%lf %lf", &time1, &angle)){
    cerr << "wrong format!" << endl;
    return 1;
  }
  double dt = time1 - time0;

  int count = 1;
  while (2 == fscanf (fin, "%lf %lf", &time, &angle)){
    if (angle >= lower && angle <= upper && time >= start){
      fprintf (fout, "%08d \t %.12e\n", count, time - 0.5 * dt);
      count++;
    }
  }
  
  fclose (fin);
  fclose (fout);
  
  return 0;
}
