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
#include "StringSplit.h"
#include "Distribution.h"
#include "Traj.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;

#define MaxLineLength 2048

static bool myread (FILE * fp,
		    float & time,
		    ValueType & angle)
{
  char line [MaxLineLength];
  while (1){
    char* rv = fgets (line, MaxLineLength, fp);
    if (rv == NULL) return false;
    if (line[0] == '#' || line[0] == '@') {
      continue;
    }
    else {
      vector<string> words;
      StringOperation::split (string(line), words);
      if (words.size() < 2) {
	cerr << "error of reading a line!! exit" << endl;
	exit(1);
      }
      time = atof(words[0].c_str());
      angle = atof(words[1].c_str());
      break;
    }
  }
  return true;
}

static void depositMetastable (const double & angle,
			       const double & tol,
			       vector<double > & counts)
{
  counts.resize (2, 0.);
  if (angle < -180 + tol || angle > 180 - tol) {
    counts[0] = 1;
    counts[1] = 0;
  }
  else {
    counts[0] = 0;
    counts[1] = 1;
  }
}

static void calCorr (const vector<double> & counts0,
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
  std::string idfile, iffile, igxsfile, ofile, ofluxfile;
  double tol, ds, lagTime;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("tol,t", po::value<double > (&tol)->default_value (30.), "the default value of tolrence of trans")
      ("ds,s", po::value<double > (&ds)->default_value (1.), "value of ds")
      ("corr-lag-time,c", po::value<double > (&lagTime)->default_value (1.), "the lag time of corrlation")
      ("output,o", po::value<std::string > (&ofile)->default_value ("prob.resp.out"), "the output of probability distribution of the trans conformation")
      ("output-flux", po::value<std::string > (&ofluxfile)->default_value ("prob.flux.resp.out"), "the output of probability flux of the trans conformation")
      ("input-dir,d",  po::value<std::string > (&idfile)->default_value ("success.dir.name"), "the file of successful dirs")
      ("input-file,f",  po::value<std::string > (&iffile)->default_value ("angaver.xvg"), "the file name")
      ("input-gxs",  po::value<std::string > (&igxsfile)->default_value ("gxs.out"), "the gxs file name");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  ifstream fpname (idfile.c_str());
  if (!fpname){
    std::cerr << "cannot open file " << idfile << std::endl;
    return 1;
  }
  char nameline [MaxLineLength];
  float time, time2;
  double dt;
  int nLagTime;
  ValueType angle;
  ValueType gxs;
  vector<float > times;
  vector<vector<double > > counts;
  vector<vector<double > > dcounts;
  vector<vector<vector<double > > > corrsBw;
  vector<vector<vector<double > > > dcorrsBw;
  unsigned countFile = 0;

  while (fpname.getline(nameline, MaxLineLength)){
    if (nameline[0] == '#') continue;
    if (nameline[0] == '@') continue;
    char nameline1 [MaxLineLength];
    char nameline2 [MaxLineLength];
    sprintf (nameline1, "%s/%s", nameline, iffile.c_str());
    sprintf (nameline2, "%s/%s", nameline, igxsfile.c_str());
    FILE *fp = fopen (nameline1, "r");
    cout << "reading file " << nameline1 << endl;
    if (fp == NULL){
      std::cerr << "cannot open file " << nameline1 << std::endl;
      return 1;
    }
    float time0, time1;
    myread(fp, time0, angle);
    myread(fp, time1, angle);
    dt = time1 - time0;
    nLagTime = int ((lagTime + 0.5 * dt) / dt) + 1;
    if (nLagTime == 1) nLagTime ++;
    Traj traj (nLagTime);
    fclose (fp);
    fp = fopen (nameline1, "r");
    FILE *fp2 = fopen (nameline2, "r");
    if (fp == NULL){
      std::cerr << "cannot open file online " << nameline1 << std::endl;
      std::cerr << std::endl;
      return 1;
    }
    if (fp2 == NULL){
      std::cerr << "cannot open file online " << nameline2 << std::endl;
      std::cerr << std::endl;
      return 1;
    }

    countFile ++;
    
    if (countFile == 1){
      vector<double > tmpcount (2, 0.0);
      vector<vector<double > > tmpcorr (2, tmpcount);
      while (myread(fp, time, angle)){
	bool success = myread(fp2, time2, gxs);
	if (!success || fabs(time - time2) > 0.001) {
	  cerr << "inconsistent gxs and angle files" << endl;
	  return 1;
	}
	times.push_back (time);
	depositMetastable (angle, tol, tmpcount);
	counts.push_back (tmpcount);
	traj.push_back (tmpcount);
	if (traj.full ()){
	  calCorr (traj.front(), traj.back(), tmpcorr);
	}	
	corrsBw.push_back (tmpcorr);
	for (unsigned dd = 0; dd < tmpcount.size(); ++dd){
	  tmpcount[dd] *= ds * gxs;
	}
	for (unsigned dd = 0; dd < tmpcorr.size(); ++dd){
	  for (unsigned mm = 0; mm < tmpcorr[dd].size(); ++mm){
	    tmpcorr[dd][mm] *= ds * gxs;
	  }
	}
	dcounts.push_back (tmpcount);
	dcorrsBw.push_back (tmpcorr);    
      }
    }
    else {
      unsigned countFrame = 0;
      vector<double > tmpcount (2, 0.0);
      vector<vector<double > > tmpcorr (2, tmpcount);
      while (myread(fp, time, angle)){
	if (countFrame >= times.size()){
	  cerr << "inconsistent frames" << endl;
	  return 1;
	}
	if (fabs (time - times[countFrame]) > 0.001 * time) {
	  cerr << "inconsistent time " << endl;
	  return 1;
	}
	bool success = myread(fp2, time2, gxs);
	if (!success || fabs(time - time2) > 0.001) {
	  cerr << "inconsistent gxs and angle files" << endl;
	  return 1;
	}	
	depositMetastable (angle, tol, tmpcount);
	for (unsigned dd = 0; dd < tmpcount.size(); ++dd){
	  counts[countFrame][dd] += tmpcount[dd];
	  dcounts[countFrame][dd] += tmpcount[dd] * ds * gxs;
	}
	traj.push_back (tmpcount);
	if (traj.full ()){
	  calCorr (traj.front(), traj.back(), tmpcorr);
	}	
	for (unsigned dd = 0; dd < tmpcorr.size(); ++dd){
	  for (unsigned mm = 0; mm < tmpcorr[dd].size(); ++mm){
	    corrsBw[countFrame][dd][mm] += tmpcorr[dd][mm];
	    dcorrsBw[countFrame][dd][mm] += tmpcorr[dd][mm] * ds * gxs;
	  }
	}
	countFrame ++;
      }
    }
    fclose (fp);
    fclose (fp2);
  }

  FILE * fp = fopen (ofile.c_str(), "w");
  FILE * fp3 = fopen (ofluxfile.c_str(), "w");

  for (unsigned ii = 0; ii < times.size(); ++ii){
    for (unsigned dd = 0; dd < 2; ++dd){
      counts[ii][dd] /= double(countFile);
      dcounts[ii][dd] /= double(countFile);
    }
  }
  for (unsigned ii = nLagTime - 1; ii < times.size(); ++ii){
    for (unsigned dd = 0; dd < 2; ++dd){
      for (unsigned mm = 0; mm < 2; ++mm){
	corrsBw[ii][dd][mm] /= double(countFile);
	dcorrsBw[ii][dd][mm] /= double(countFile);
      }
    }
  }

  for (unsigned ii = 0; ii < times.size(); ++ii){
    fprintf (fp, "%f ", times[ii]);
    for (unsigned dd = 0; dd < 2; ++dd){
      fprintf (fp, "%f %f  ",
	       counts[ii][dd],
	       dcounts[ii][dd]
	  );
    }
    fprintf (fp, "\n");
  }
  // for (unsigned ii = 0; ii < times.size(); ++ii){
  //   fprintf (fp, "%f %f\n", times[ii], counts[ii] / double(countFile));
  // }
  // fclose (fp);
  
  for (unsigned ii = nLagTime - 1; ii < times.size(); ++ii){
    fprintf (fp3, "%f ", times[ii]);
    for (unsigned dd = 0; dd < 2; ++dd){
      for (unsigned mm = 0; mm < 2; ++mm){
	fprintf (fp3, "%f %f  ",
		 (corrsBw[ii][dd][mm] - corrsBw[ii][mm][dd]) / lagTime,
		 (dcorrsBw[ii][dd][mm] - dcorrsBw[ii][mm][dd]) / lagTime
	    );
      }
      fprintf (fp3, "  ");
    }
    fprintf (fp3, "\n");
  }

  fclose (fp);
  fclose (fp3);

  return 0;
}
