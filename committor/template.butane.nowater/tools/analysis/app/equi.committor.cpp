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
  std::string itfile, icfile, ofile;
  unsigned nDataBlock;
  double begin;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("begin,b",  po::value<double > (&begin)->default_value (0.), "start of the traj")
      ("num-data-block,n",  po::value<unsigned > (&nDataBlock)->default_value (1), "number of data in block")
      ("input-angle,a", po::value<std::string > (&itfile)->default_value ("angaver.xvg"), "the angle file name")
      ("input-coreset-file,f", po::value<std::string > (&icfile)->default_value ("coreset.dat"), "the file of core sets")
      ("output,o", po::value<std::string > (&ofile)->default_value ("committor.out"), "the output of the committor");
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  AngleSetTable1D stable;
  stable.reinit (icfile);
  
  int nbin = stable.size();
  double hbin = 360./(double(nbin));

  vector<BlockAverage_acc> bas (nbin);
  for (int ii = 0; ii < nbin; ++ii){
    bas[ii].reinit (nDataBlock);
  }

  ifstream angle (itfile.c_str());
  if (!angle){
    cerr << "cannot open angle file " << angle << endl;
    return 1;
  }
  
  char valueline_ang [MaxLineLength];
  int prev_indicate = -10;
  vector<int > index_traj;
  unsigned countTraj = 0;
  
  while (angle.getline(valueline_ang, MaxLineLength)){
    if (valueline_ang[0] == '#' || valueline_ang[0] == '@'){
      continue;
    }
    vector<string > words;
    StringOperation::split (string(valueline_ang), words);
    if (words.size() < 2) {
      cerr << "wrong file format of " << itfile << endl;
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
    
    int index = stable.calIndex (anglev);
    int indicate = stable.calIndicate (index);

    if (indicate == 0){
      index_traj.push_back (index);
    }
    else {
      if (prev_indicate == 0){ // end of a nontrivial traj
	if (indicate == 1){ // hit!
	  for (unsigned ii = 0; ii < index_traj.size(); ++ii){
	    bas[index_traj[ii]].deposite (1.);
	  }
	}
	if (indicate == -1){ // goes back
	  for (unsigned ii = 0; ii < index_traj.size(); ++ii){
	    bas[index_traj[ii]].deposite (0.);
	  }
	}
	// clear traj
	index_traj.clear ();
	index_traj.reserve (1000);
      }
    }
    
    prev_indicate = indicate;
  }

  for (int ii = 0; ii < nbin; ++ii){
    bas[ii].calculate ();
  }

  FILE * fp = fopen (ofile.c_str(), "w");
  for (int ii = 0; ii < nbin; ++ii){
    double xx = -180 + 0.5 * hbin + ii * hbin;
    double vv, ee;
    double var;
    if (stable.calIndicate(ii) == -1){
      vv = 0;
      ee = 0;
      var = 0;
    }
    else if (stable.calIndicate(ii) == 1){
      vv = 1;
      ee = 0;
      var = 0;
    }
    else {
      vv = bas[ii].getAvg();
      ee = bas[ii].getAvgError();
      var = bas[ii].getVar();
    }
    fprintf (fp, "%f %e %e %e\n", xx, vv, var, ee);
  }
  fclose (fp);
  
  return 0;
}

