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

#include "Defines.h"
#include "Distribution.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;
#define MaxLineLength 2048

bool myread (FILE * fp,
	     float & time,
	     ValueType & phi,
	     ValueType & psi)
{
  size_t rv;
  rv = fread (&time, sizeof(float), 1, fp);
  if (rv != 1){
    // cerr << "error read time or reach EOF" << endl;
    return false;
  }
  rv = fread (&phi, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error read phi " << endl;
    exit(1);
  }
  rv = fread (&psi, sizeof(double), 1, fp);
  if (rv != 1){
    cerr << "error read psi " << endl;
    exit(1);
  }
  return true;
}

bool myreadGxs (FILE * fp,
		float & time,
		ValueType & gxs)
{
  size_t rv;
  rv = fscanf (fp, "%f", &time);
  if (rv != 1){
    // cerr << "error read time or reach EOF" << endl;
    return false;
  }
  rv = fscanf (fp, "%lf", &gxs);  
  if (rv != 1){
    cerr << "error read gxs " << endl;
    exit(1);
  }
  return true;
}

void depositMetastable (const double & phi,
			const double & psi,
			const vector<MetastableSet> & sets,
			vector<double > & counts)
{
  int cc = 0;
  counts.resize (sets.size());
  for (unsigned ii = 0; ii < sets.size(); ++ii){
    counts[ii] = 0;
    if (sets[ii].inSet(phi, psi)){
      counts[ii] = 1;
      cc ++;
    }
  }
  
  if (cc != 1) {
    if (cc == 0){
      cerr << "not in any set!!!" << endl;
    }
    else {
      cerr << "in multi sets" << endl;
    }
  }
}


int main(int argc, char * argv[])
{
  std::string ifile, ofile, gfile;
  double ds;
  // unsigned numBlock = 20;
  // double refh;
  float time_prec = .01;
  double setA_psi_b = 128, setA_psi_e = 13;
  double setA_phi_b =-125, setA_phi_e = 74;
  double setB_psi_b = 128, setB_psi_e = 13;
  double setB_phi_b = 74,  setB_phi_e =-125;
  double setC_psi_b = 13,  setC_psi_e = 128;
  double setC_phi_b = -180,setC_phi_e = 180;
  
  double setA1_psi_b =-134, setA1_psi_e = 13;
  double setA1_phi_b =-125, setA1_phi_e = 74;
  double setA2_psi_b = 128, setA2_psi_e =-134;
  double setA2_phi_b =-125, setA2_phi_e = 74;

  double setB1_psi_b =-110, setB1_psi_e = 13;
  double setB1_phi_b = 74,  setB1_phi_e =-125;
  double setB2_psi_b = 128, setB2_psi_e =-110;
  double setB2_phi_b = 74,  setB2_phi_e =-125;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("ds,s", po::value<double > (&ds)->default_value (1.), "value of ds")
      ("input,f",  po::value<std::string > (&ifile)->default_value ("angle.name"), "the file of file names: angle")
      ("gxs-file,g",  po::value<std::string > (&gfile)->default_value ("gxs.name"), "the file of file names: gxs")
      ("output,o", po::value<std::string > (&ofile)->default_value ("metastable.resp.out"), "the output of metastable propulation");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  // MetastableSet setA 
  // MetastableSet setB (setB_psi_b, setB_psi_e, setB_phi_b, setB_phi_e);
  // MetastableSet setC (setC_psi_b, setC_psi_e, setC_phi_b, setC_phi_e);
  vector<MetastableSet > sets;
  // sets.push_back (MetastableSet(setA_psi_b, setA_psi_e, setA_phi_b, setA_phi_e));
  // sets.push_back (MetastableSet(setB_psi_b, setB_psi_e, setB_phi_b, setB_phi_e));
  sets.push_back (MetastableSet(setA1_psi_b, setA1_psi_e, setA1_phi_b, setA1_phi_e));
  sets.push_back (MetastableSet(setA2_psi_b, setA2_psi_e, setA2_phi_b, setA2_phi_e));
  sets.push_back (MetastableSet(setB1_psi_b, setB1_psi_e, setB1_phi_b, setB1_phi_e));
  sets.push_back (MetastableSet(setB2_psi_b, setB2_psi_e, setB2_phi_b, setB2_phi_e));
  sets.push_back (MetastableSet(setC_psi_b, setC_psi_e, setC_phi_b, setC_phi_e));

  ifstream fpname (ifile.c_str());
  if (!fpname){
    std::cerr << "cannot open file " << ifile << std::endl;
    return 1;
  }
  ifstream fpname1 (gfile.c_str());
  if (!fpname1){
    std::cerr << "cannot open file " << gfile << std::endl;
    return 1;
  }
  char nameline [MaxLineLength];
  char nameline1 [MaxLineLength];
  float time, time1;
  ValueType phi, psi, gxs;
  vector<float > times;
  vector<vector<double > > counts;
  vector<vector<double > > dcounts;
  unsigned countFile = 0;

  while (fpname.getline(nameline, MaxLineLength) &&
	 fpname1.getline(nameline1, MaxLineLength)){
    if (nameline[0] == '#') continue;
    if (nameline1[0] == '#') continue;
    FILE *fp = fopen (nameline, "r");
    FILE *fp1 = fopen (nameline1, "r");
    cout << "reading file " << nameline
	 << " " << nameline1<< endl;
    if (fp == NULL){
      std::cerr << "cannot open file " << nameline << std::endl;
      return 1;
    }
    if (fp1 == NULL){
      std::cerr << "cannot open file " << nameline1 << std::endl;
      return 1;
    }
    countFile ++;
    vector<double > tmpcount;
    if (countFile == 1){
      while (myread(fp, time, phi, psi)){
	if (!myreadGxs(fp1, time1, gxs) || fabs(time - time1) > time_prec){
	  cerr << "inconsistent files" << endl;
	  return (1);
	}
	times.push_back (time);
	depositMetastable (psi, phi, sets, tmpcount);
	counts.push_back (tmpcount);
	for (unsigned dd = 0; dd < tmpcount.size(); ++dd){
	  tmpcount[dd] *= ds * gxs;
	}
	dcounts.push_back (tmpcount);
      }	
    }
    else {
      unsigned countFrame = 0;
      while (myread(fp, time, phi, psi)){
	if (countFrame >= times.size()){
	  cerr << "inconsistent frames" << endl;
	  return 1;
	}
	if (!myreadGxs(fp1, time1, gxs) || fabs(time - time1) > time_prec){
	  cerr << "inconsistent files" << endl;
	  return (1);
	}
	depositMetastable (psi, phi, sets, tmpcount);
	for (unsigned dd = 0; dd < tmpcount.size(); ++dd){
	  counts[countFrame][dd] += tmpcount[dd];
	  dcounts[countFrame][dd] += tmpcount[dd] * ds * gxs;
	}
	countFrame ++;
      }
    }
    fclose (fp);
    fclose (fp1);
  }

  FILE * fp = fopen (ofile.c_str(), "w");
  for (unsigned ii = 0; ii < times.size(); ++ii){
    fprintf (fp, "%f ", times[ii]);
    for (unsigned dd = 0; dd < sets.size(); ++dd){
      fprintf (fp, "%f %f %f   ",
	       (counts[ii][dd] + dcounts[ii][dd]) / double(countFile),
	       (dcounts[ii][dd]) / double(countFile),
	       ( counts[ii][dd]) / double(countFile)
	  );
    }
    fprintf (fp, "\n");
  }
  fclose (fp);
  
  return 0;
}
