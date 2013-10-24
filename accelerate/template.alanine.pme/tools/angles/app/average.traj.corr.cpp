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
#include "Traj.h"
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

void calCorr (const vector<double> & counts0,
	      const vector<double> & counts1,
	      vector<vector<double > > & countCorr)
{
  countCorr.resize (counts0.size());
  for (unsigned ii = 0; ii < countCorr.size(); ++ii){
    countCorr[ii].resize (counts1.size());
    for (unsigned jj = 0; jj < countCorr[ii].size(); ++jj){
      if (counts1[ii] != 0 && counts0[jj] != 0){
	countCorr[ii][jj] = 1.;
      }
      else {
	countCorr[ii][jj] = 0.;
      }
    }
  }
}

int main(int argc, char * argv[])
{
  std::string ifile, ofile, ofwfile, obwfile, ofluxfile, odirfile;
  double lagTime;
  // unsigned numBlock = 20;
  // double refh;
  // float time_prec = .01;

  // double setA_psi_b = 128, setA_psi_e = 13;
  // double setA_phi_b =-125, setA_phi_e = 74;
  // double setB_psi_b = 128, setB_psi_e = 13;
  // double setB_phi_b = 74,  setB_phi_e =-125;
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
      ("corr-lag-time,c", po::value<double > (&lagTime)->default_value (10.), "the lag time of corrlation")
      ("output-meta,o", po::value<std::string > (&ofile)->default_value ("metastable.out"), "the output of metastable propulation")
      ("output-meta-corr-forward", po::value<std::string > (&ofwfile)->default_value ("meta.corr.fw.out"), "the output of metastable propulation forward correlation")
      ("output-meta-corr-backward", po::value<std::string > (&obwfile)->default_value ("meta.corr.bw.out"), "the output of metastable propulation backward correlation")
      ("output-meta-flux", po::value<std::string > (&ofluxfile)->default_value ("meta.flux.out"), "the output of metastable propulation flux")
      ("output-meta-direct-flux", po::value<std::string > (&odirfile)->default_value ("meta.direct.flux.out"), "the output of metastable propulation flux")
      ("input,f",  po::value<std::string > (&ifile)->default_value ("angle.name"), "the file of file names");

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
    std::cerr << "cannot open input file " << ifile << std::endl;
    std::cerr << std::endl;
    return 1;
  }
  char nameline [MaxLineLength];
  float time;
  double dt;
  int nLagTime;
  ValueType phi, psi;
  vector<float > times;
  vector<vector<double > > counts;
  vector<vector<vector<double > > > corrsBw;
  vector<vector<vector<double > > > corrsFw;
  unsigned countFile = 0;

  while (fpname.getline(nameline, MaxLineLength)){
    if (nameline[0] == '#') continue;
    FILE *fp = fopen (nameline, "r");
    cout << "reading file " << nameline << endl;
    if (fp == NULL){
      std::cerr << "cannot open file online " << nameline << std::endl;
      std::cerr << std::endl;
      return 1;
    }
    float time0, time1;
    myread(fp, time0, phi, psi);
    myread(fp, time1, phi, psi);
    dt = time1 - time0;
    nLagTime = int ((lagTime + 0.5 * dt) / dt) + 1;
    if (nLagTime == 1) nLagTime ++;
    Traj traj (nLagTime);
    fclose (fp);
    fp = fopen (nameline, "r");
    if (fp == NULL){
      std::cerr << "cannot open file online " << nameline << std::endl;
      std::cerr << std::endl;
      return 1;
    }
    
    countFile ++;

    if (countFile == 1){
      vector<double > tmpcount (sets.size(), 0.0);
      vector<vector<double > > tmpcorr (sets.size(), tmpcount);
      while (myread(fp, time, phi, psi)){
	times.push_back (time);
	depositMetastable (psi, phi, sets, tmpcount);
	counts.push_back (tmpcount);
	traj.push_back (tmpcount);
	if (traj.full ()){
	  calCorr (traj.front(), traj.back(), tmpcorr);
	}	
	corrsBw.push_back (tmpcorr);
      }
    }
    else {
      unsigned countFrame = 0;
      vector<double > tmpcount (sets.size(), 0.0);
      vector<vector<double > > tmpcorr (sets.size(), tmpcount);
      while (myread(fp, time, phi, psi)){
	if (countFrame >= times.size()){
	  cerr << "inconsistent frames" << endl;
	  return 1;
	}
	if (fabs (time - times[countFrame]) > 0.001 * time) {
	  cerr << "inconsistent time " << endl;
	  return 1;
	}
	depositMetastable (psi, phi, sets, tmpcount);
	for (unsigned dd = 0; dd < tmpcount.size(); ++dd){
	  counts[countFrame][dd] += tmpcount[dd];
	}
	traj.push_back (tmpcount);
	if (traj.full ()){
	  calCorr (traj.front(), traj.back(), tmpcorr);
	}	
	for (unsigned dd = 0; dd < tmpcorr.size(); ++dd){
	  for (unsigned mm = 0; mm < tmpcorr[dd].size(); ++mm){
	    corrsBw[countFrame][dd][mm] += tmpcorr[dd][mm];
	  }
	}
	countFrame ++;
      }
    }
    fclose (fp);
  }

  FILE * fp = fopen (ofile.c_str(), "w");
  FILE * fp1 = fopen (obwfile.c_str(), "w");
  FILE * fp2 = fopen (ofwfile.c_str(), "w");
  FILE * fp3 = fopen (ofluxfile.c_str(), "w");
  FILE * fp4 = fopen (odirfile.c_str(), "w");

  for (unsigned ii = 0; ii < times.size(); ++ii){
    fprintf (fp, "%f ", times[ii]);
    for (unsigned dd = 0; dd < sets.size(); ++dd){
      counts[ii][dd] = counts[ii][dd] / double(countFile);
      fprintf (fp, "%f ", counts[ii][dd]);
    }
    fprintf (fp, "\n");
  }

  for (unsigned ii = 0; ii < times.size(); ++ii){
    if (int(ii) < nLagTime-1) continue;
    fprintf (fp1, "%f ", times[ii]);
    for (unsigned dd = 0; dd < sets.size(); ++dd){
      if (counts[ii][dd] != 0){
	for (unsigned mm = 0; mm < sets.size(); ++mm){
	  fprintf (fp1, "%f ", corrsBw[ii][dd][mm] / counts[ii][dd] / double (countFile));
	}
      }
      else {
	for (unsigned mm = 0; mm < sets.size(); ++mm){
	  fprintf (fp1, "%f ", 0.);
	}
      }
      fprintf (fp1, "  ");
    }
    fprintf (fp1, "\n");
  }

  for (unsigned ii = 0; ii < times.size()-nLagTime+1; ++ii){
    fprintf (fp2, "%f ", times[ii]);
    for (unsigned dd = 0; dd < sets.size(); ++dd){
      if (counts[ii+nLagTime-1][dd] != 0){
	for (unsigned mm = 0; mm < sets.size(); ++mm){
	    fprintf (fp2, "%f ", corrsBw[ii+nLagTime-1][mm][dd] / counts[ii][dd] / double (countFile));
	}
      }
      else{
	for (unsigned mm = 0; mm < sets.size(); ++mm){
	  fprintf (fp2, "%f ", 0.);
	}
      }
      fprintf (fp2, "  ");
      }
    fprintf (fp2, "\n");
  }
  

  for (unsigned ii = nLagTime - 1; ii < times.size(); ++ii){
    fprintf (fp3, "%f ", times[ii]);
    for (unsigned dd = 0; dd < sets.size(); ++dd){
      for (unsigned mm = 0; mm < sets.size(); ++mm){
	fprintf (fp3, "%f ", (corrsBw[ii][dd][mm] - corrsBw[ii][mm][dd]) / double (countFile) / lagTime);
      }
      fprintf (fp3, "  ");
    }
    fprintf (fp3, "\n");
  }
  
  for (unsigned ii = nLagTime - 1; ii < times.size(); ++ii){
    fprintf (fp4, "%f ", times[ii]);
    for (unsigned dd = 0; dd < sets.size(); ++dd){
      for (unsigned mm = 0; mm < sets.size(); ++mm){
	fprintf (fp4, "%f ", (corrsBw[ii][dd][mm]) / double (countFile) / lagTime);
      }
      fprintf (fp4, "  ");
    }
    fprintf (fp4, "\n");
  }
  


  fclose (fp);
  fclose (fp1);
  fclose (fp2);
  fclose (fp3);
  fclose (fp4);
  
  return 0;
}
