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
  std::string ofile;
  double kk, phi, nn, bin;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("kphi,k", po::value<double > (&kk)->default_value (5.92), "the k value")
      ("phi,p",  po::value<double > (&phi)->default_value (0),   "the \\phi value")
      ("multi,n",  po::value<double > (&nn)->default_value (3),   "the n value")
      ("bin,r",  po::value<double > (&bin)->default_value (.1),   "the bin size")
      ("output,o", po::value<std::string > (&ofile)->default_value ("table_d0.xvg"), "the output of count of h-bond");
  
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  unsigned ndata = (360. + 0.01 * bin) / bin;
  bin = 360. / double(ndata);
  ndata ++;
  cout << "# print table with bin size: " << bin << endl;
  cout << "# ndata: " << ndata << endl;

  FILE * fout = fopen (ofile.c_str(), "w");
  if (fout == NULL){
    cout << "cannot open file " << ofile<< endl;
    exit (1);
  }
  fprintf (fout, "# table of dihedral gmx funct1. k: %f phi: %f n: %f\n",
	   kk, phi, nn);
  for (unsigned ii = 0; ii < ndata; ++ii){
    double myangle = -180 + ii * bin;
    double tmpa = (nn * myangle - phi) / 180. * M_PI;
    double mv = kk * (1. + cos(tmpa));
    double md = nn * kk * sin(tmpa) / 180. * M_PI;
    fprintf (fout, "%f %e %e\n", myangle, mv, md);
  }
  
  fclose (fout);
  
  return 0;
}
