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
  double kk, uu, ss, bin;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("max,k", po::value<double > (&kk)->default_value (10.), "the maximum value")
      ("mu,u",  po::value<double > (&uu)->default_value (60.),   "the center of the potential")
      ("sigma,s",  po::value<double > (&ss)->default_value (25.),   "the value of sigma")
      ("bin,r",  po::value<double > (&bin)->default_value (.1),   "the bin size")
      ("output,o", po::value<std::string > (&ofile)->default_value ("table_d0.xvg"), "the output of dihedral potential");
      
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
  fprintf (fout, "# table of gaussian shaped potential max: %f mu: %f sigma: %f\n",
	   kk, uu, ss);
  for (unsigned ii = 0; ii < ndata; ++ii){
    double myangle = -180 + ii * bin;
    double dist = myangle - uu;
    if (fabs(dist + 360.) < fabs(dist)) {
      dist = dist + 360.;
    }
    else if (fabs(dist - 360) < fabs(dist)) {
      dist = dist - 360.;
    }
    double tmp = exp ( - (dist * dist) / (2. * ss * ss));
    double mv = kk * tmp;
    double md = kk * tmp * (dist / (ss * ss));
    fprintf (fout, "%f %e %e\n", myangle, mv, md);
  }
  
  fclose (fout);
  
  return 0;
}
