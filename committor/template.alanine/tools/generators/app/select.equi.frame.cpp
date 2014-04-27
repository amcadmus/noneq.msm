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
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_trr.h"

#define MaxLineLength 2048

namespace po = boost::program_options;
using namespace std;

int main(int argc, char * argv[])
{
  std::string icfile, phifile, psifile, itfile, ofile;
  double begin;
  int targetInd;
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("begin,b",  po::value<double > (&begin)->default_value (0.), "start of the traj")
      ("target-indicator,i", po::value<int > (&targetInd)->default_value (9), "the indicator to select")
      ("input-coreset,c", po::value<string > (&icfile)->default_value ("coreset.dat"), "input selection map")
      ("input-phi", po::value<string > (&phifile)->default_value ("angaver.phi.xvg"), "input angle phi")
      ("input-psi", po::value<string > (&psifile)->default_value ("angaver.psi.xvg"), "input angle psi")
      ("input-traj,f", po::value<string > (&itfile)->default_value ("traj.trr"), "input trajectory")
      ("output,o", po::value<string > (&ofile)->default_value ("equi.frame.trr"), "output equi frames");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  AngleSetTable2D stable ;
  stable.reinit (icfile);

  XDRFILE * xdrfp = xdrfile_open (itfile.c_str(), "r");
  if (xdrfp == NULL) {
    cerr << "cannot open traj file " << itfile << endl;
    return 1;
  }
  int natoms;
  char itfilename [MaxLineLength];
  strcpy (itfilename, itfile.c_str());
  read_trr_natoms (itfilename, &natoms);
  int step;
  float trajtime;
  float lambda;
  matrix box;
  rvec *xx = (rvec *) malloc (sizeof(rvec) * natoms);
  rvec *vv = (rvec *) malloc (sizeof(rvec) * natoms);
  rvec *ff = (rvec *) malloc (sizeof(rvec) * natoms);
  
  ifstream anglephi (phifile.c_str());
  if (!anglephi) {
    cerr << "cannot open angle file " << phifile << endl;
    return 1;
  }
  ifstream anglepsi (psifile.c_str());
  if (!anglepsi) {
    cerr << "cannot open angle file " << psifile << endl;
    return 1;
  }

  char valueline_angphi [MaxLineLength];
  char valueline_angpsi [MaxLineLength];
  int countTraj = 0;
  int countSel = 0;
  XDRFILE * fpo = xdrfile_open (ofile.c_str(), "w");

  while (anglephi.getline(valueline_angphi, MaxLineLength) &&
	 anglepsi.getline(valueline_angpsi, MaxLineLength)){
    if (valueline_angphi[0] == '#' || valueline_angphi[0] == '@'){
      continue;
    }
    if (valueline_angpsi[0] == '#' || valueline_angpsi[0] == '@'){
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

    int rv = read_trr (xdrfp, natoms, &step, &trajtime, &lambda, box, xx, vv, ff);
    if (rv != exdrOK){
      cerr << "error of reading trr file, return " << endl;
      return 1;
    }
    
    if (fabs(trajtime - timephi) > 1e-3){
      cerr << "inconsistent time in the traj file and the angle file, something goes wrong, exit" << endl;
      return 1;
    }

    if (timephi < begin - 1e-3){
      continue;
    }
    else {
      if (countTraj %1000 == 0){
	// printf ("read traj at %f     \r", timephi);
	// fflush (stdout);
      }
    }

    if (stable.calIndicate (anglevphi, anglevpsi) == targetInd){
      countSel ++;
      printf ("select time %f, angle: %f %f, print as %f\n",
	      timephi,
	      anglevphi, anglevpsi,
	      float(countSel));
      rv = write_trr (fpo, natoms, countSel, float(countSel), lambda, box, xx, vv, ff);
      if (rv != exdrOK){
	cerr << "error of writing trr file, return " << endl;
	return 1;
      }
    }
  }  

  xdrfile_close (xdrfp);
  xdrfile_close (fpo);
  free (xx);
  free (vv);
  free (ff);
      
  return 0;
}
