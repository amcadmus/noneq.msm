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
  std::string phifile, psifile, icfile, ofile;
  unsigned nDataBlock;
  double begin;
  int toind, notinind;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("begin,b",  po::value<double > (&begin)->default_value (0.), "start of the traj")
      ("num-data-block,n",  po::value<unsigned > (&nDataBlock)->default_value (1), "number of data in block")
      ("input-angle-phi", po::value<std::string > (&phifile)->default_value ("angaver.phi.xvg"), "the angle file name")
      ("input-angle-psi", po::value<std::string > (&psifile)->default_value ("angaver.psi.xvg"), "the angle file name")
      ("input-coreset-file,c", po::value<std::string > (&icfile)->default_value ("coreset.dat"), "the file of core sets")
      ("to-indicator",   po::value<int > (&toind  )->default_value (3), "to core set with this indicator")
      ("notin-indicator",   po::value<int > (&notinind)->default_value (0), "indicator for not in any core set")
      ("output,o", po::value<std::string > (&ofile)->default_value ("committor.out"), "the output of the committor");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  AngleSetTable2D stable;
  stable.reinit (icfile);
  
  int nbinphi = stable.sizePhi();
  int nbinpsi = stable.sizePsi();
  double hbinphi = 360./(double(nbinphi));
  double hbinpsi = 360./(double(nbinpsi));

  vector<BlockAverage_acc > bas (nbinphi * nbinpsi);
  for (int ii = 0; ii < nbinphi * nbinpsi; ++ii){
    bas[ii].reinit (nDataBlock);
  }

  ifstream anglephi (phifile.c_str());
  if (!anglephi){
    cerr << "cannot open angle file " << phifile << endl;
    return 1;
  }
  ifstream anglepsi (psifile.c_str());
  if (!anglepsi){
    cerr << "cannot open angle file " << psifile << endl;
    return 1;
  }
  
  char valueline_angphi [MaxLineLength];
  char valueline_angpsi [MaxLineLength];
  int prev_indicate = -10;
  vector<int > index_traj;
  unsigned countTraj = 0;
  
  while (anglephi.getline(valueline_angphi, MaxLineLength) &&
	 anglepsi.getline(valueline_angpsi, MaxLineLength)){
    if (valueline_angphi[0] == '#' || valueline_angphi[0] == '@' ||
	valueline_angpsi[0] == '#' || valueline_angpsi[0] == '@' ){
      // cout << "angle file format problem..." << endl;
      // return 1;
      continue;
    }
    double timephi, anglevphi;
    double timepsi, anglevpsi;
    vector<string > words;
    StringOperation::split (string(valueline_angphi), words);
    if (words.size() < 2) {
      cerr << "wrong file format of " << phifile << endl;
      exit (1);
    }
    timephi   = (atof(words[0].c_str()));
    anglevphi = (atof(words[1].c_str()));
    StringOperation::split (string(valueline_angpsi), words);
    if (words.size() < 2) {
      cerr << "wrong file format of " << psifile << endl;
      exit (1);
    }
    timepsi   = (atof(words[0].c_str()));
    anglevpsi = (atof(words[1].c_str()));
    if (fabs(timephi - timepsi) > 1e-3){
      cerr << "inconsistent time in the phi and the psi angle file, something goes wrong, exit" << endl;
      return 1;
    }
    
    if (timephi < begin - 1e-3){
      continue;
    }
    else {
      if (countTraj %1000 == 0){
	printf ("read traj at %f     \r", timephi);
	fflush (stdout);
      }
    }
    countTraj ++;
    
    int index = stable.calIndex (anglevphi, anglevpsi);
    int indicate = stable.calIndicate (index);
    // printf ("%f %d   ", timephi, indicate);
    // printf ("%f %d \n", timephi, indicate);

    if (indicate == notinind){
      index_traj.push_back (index);
    }
    else {
      if (prev_indicate == notinind){ // end of a nontrivial traj
	if (indicate == toind){ // hit!
	  // printf (" hit ");
	  for (unsigned ii = 0; ii < index_traj.size(); ++ii){
	    // printf (" %d ", stable.calIndicate(index_traj[ii]));
	    bas[index_traj[ii]].deposite (1.);
	  }
	}
	if (indicate != toind){ // goes back or hit a wrong set
	  for (unsigned ii = 0; ii < index_traj.size(); ++ii){
	    bas[index_traj[ii]].deposite (0.);
	  }
	}
	// clear traj
	index_traj.clear ();
	index_traj.reserve (1000);
      }
    }
    // printf ("\n");
    
    prev_indicate = indicate;
  }
  printf ("\n");

  for (int ii = 0; ii < nbinphi * nbinpsi; ++ii){
    bas[ii].calculate ();
  }

  FILE * fp = fopen (ofile.c_str(), "w");
  for (int ii = 0; ii < nbinphi; ++ii){
    double xx = -180 + 0.5 * hbinphi + ii * hbinphi;
    for (int jj = 0; jj < nbinpsi; ++jj){
      int index = ii * nbinpsi + jj;
      double yy = -180 + 0.5 * hbinpsi + jj * hbinpsi;
      double vv, ee;
      double var;
      if (stable.calIndicate(index) != toind && stable.calIndicate(index) != notinind){
	vv = 0;
	ee = 0;
	var = 0;
      }
      else if (stable.calIndicate(index) == toind){
	vv = 1;
	ee = 0;
	var = 0;
      }
      else {
	vv = bas[index].getAvg();
	ee = bas[index].getAvgError();
	var = bas[index].getVar();
      }
      fprintf (fp, "%f %f  %e %e %e\n", xx, yy, vv, var, ee);
    }
    fprintf (fp, "\n");
  }
  fclose (fp);
  
  return 0;
}

