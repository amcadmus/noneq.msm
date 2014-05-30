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

#include "Distribution.h"
#include "MetastableSet.h"

#define MaxLineLength 2048

namespace po = boost::program_options;
using namespace std;

int main(int argc, char * argv[])
{
  std::string phifile, psifile, ofile;
  double hbinx, hbiny;
  double threshold, begin, end;
  double phil, phiu, psil, psiu;
  int printvalue;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print coreset from Ramanchandran plot")
      ("begin,b",  po::value<double > (&begin)->default_value (0.), "start of the traj")
      ("end,e",  po::value<double > (&end)->default_value (0.), "end of the traj")
      ("threshold,t",  po::value<double > (&threshold)->default_value (0.01), "start of the traj")
      ("print-indicator,i",  po::value<int > (&printvalue)->default_value (1), "printed value in the core set file")
      ("phi-low", po::value<double > (&phil)->default_value (-180), "upper range of phi")
      ("phi-up",  po::value<double > (&phiu)->default_value ( 180), "lower range of phi")
      ("psi-low", po::value<double > (&psil)->default_value (-180), "upper range of psi")
      ("psi-up",  po::value<double > (&psiu)->default_value ( 180), "lower range of psi")
      ("input-angle-phi", po::value<std::string > (&phifile)->default_value ("angaver.phi.xvg"), "the angle file name")
      ("input-angle-psi", po::value<std::string > (&psifile)->default_value ("angaver.psi.xvg"), "the angle file name")
      ("bin-size-phi", po::value<double > (&hbinx)->default_value (15), "width of the set")
      ("bin-size-psi", po::value<double > (&hbiny)->default_value (15), "width of the set")
      ("output,o", po::value<std::string > (&ofile)->default_value ("coreset.dat"), "file of coresets");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  int nbinx = (360. + 0.5 * hbinx)/(double(hbinx));
  int nbiny = (360. + 0.5 * hbiny)/(double(hbiny));

  Distribution_1d dist (-180, 180, nbinx, -180, 180, nbiny);

  FILE * fpphi = fopen (phifile.c_str(), "r");
  if (fpphi == NULL) {
    cerr << "cannot open file " << phifile << endl;
    return 1;
  }
  FILE * fppsi = fopen (psifile.c_str(), "r");
  if (fppsi == NULL) {
    cerr << "cannot open file " << psifile << endl;
    return 1;
  }

  double timephi, anglephi;
  double timepsi, anglepsi;
  int count = 0;
  
  while ((2 == fscanf (fpphi, "%lf %lf", &timephi, &anglephi)) &&
	 (2 == fscanf (fppsi, "%lf %lf", &timepsi, &anglepsi)) )  {
    if (timephi != timepsi){
      cerr << "inconsistent time of phi and psi, exit" << endl;
      return 1;
    }
    if (timephi < begin - 1e-3){
      continue;
    }
    if (end > 0 && timephi > end + 1e-3){
      break;
    }
    if (count++ %1000 == 0){
      printf ("read angle at %f     \r", timephi);
      fflush (stdout);
    }
    dist.deposite (anglephi, anglepsi);
  }
  printf ("\nread until %f\n", timephi);

  dist.average();
  dist.print_xv ("equi.dist.out");
  
  double sum = 0;
  for (int ii = 0; ii < nbinx; ++ii){
    for (int jj = 0; jj < nbiny; ++jj){
      sum += dist.values[ii][jj] * dist.hx * dist.hv;
    }
  }
  printf ("sum is %f, number is %f, count is %d\n", sum, dist.nframe, count);

  MetastableSet ms (phil, phiu, psil, psiu);
  
  double sum1 = 0.;
  for (int ii = 0; ii < nbinx; ++ii){
    for (int jj = 0; jj < nbiny; ++jj){
      if (ms.inSet (dist.gridx[ii], dist.gridv[jj])){
	sum1 += dist.values[ii][jj] * dist.hx * dist.hv;
      }
    }
  }
  
  
  FILE * fp = fopen (ofile.c_str(), "w");
  
  double sum2 = 0;
  for (int ii = 0; ii < nbinx; ++ii){
    for (int jj = 0; jj < nbiny; ++jj){
      int ind ;
      if (ms.inSet (dist.gridx[ii], dist.gridv[jj])){
	// if (phil <= dist.gridx[ii] && dist.gridx[ii] < phiu &&
	// 	  psil <= dist.gridv[jj] && dist.gridv[jj] < psiu ) {
	sum2 += dist.values[ii][jj] * dist.hx * dist.hv / sum1;
	if (dist.values[ii][jj] * dist.hx * dist.hv / sum1 > threshold) {
	  ind = printvalue;
	}
	else {
	  ind = 0;
	}
      }
      else {
	ind  = 0;
      }
      // if (ind!=0){
      // 	// fprintf (fp, "%d  ", ind);
      // 	fprintf (fp, "%f  ", dist.values[ii][jj] * dist.hx * dist.hv / sum1);
      // }
      // else{
      fprintf (fp, "%d  ", ind);
      // }
    }
    fprintf (fp, "\n");
  }
  printf ("sum2 is %f\n", sum2);
  
  return 0;
}

