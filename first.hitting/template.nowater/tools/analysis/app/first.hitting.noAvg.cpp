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
#include "MetastableSet.h"
#include "StringSplit.h"
#include "BlockAverage.h"

#define MaxLineLength 2048

namespace po = boost::program_options;
using namespace std;

inline int angleIdx (const double & myangle,
		     const double & bin)
{
  return int((myangle + 180) / bin);
}


int main(int argc, char * argv[])
{
  std::string iffile, idfile, ofile;
  double center, width;
  double gate;
  unsigned nDataInBlock;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("gate,g",  po::value<double > (&gate)->default_value (5.0),   "the probability of first hitting time smaller than the gate")
      ("meta-center,c",  po::value<double > (&center)->default_value (180),   "center of the metastable set")
      ("meta-width,w",  po::value<double > (&width)->default_value (30),   "width of the metastable set")
      ("num-in-block,n",  po::value<unsigned > (&nDataInBlock)->default_value (1),   "number of data in a block")
      ("input-dir,d", po::value<std::string > (&idfile)->default_value ("success.dir.name"), "file including angle file names")
      ("input-file-name,f", po::value<std::string > (&iffile)->default_value ("angaver.xvg"), "file including traj file names")
      ("output,o", po::value<std::string > (&ofile)->default_value ("fmpt.out"), "the output of first mean passage time");
  
      
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  // unsigned nbin = (360. + 0.01 * bin) / bin;
  // bin = 360. / double(nbin);
  // cout << "# bin size: " << bin << endl;
  // cout << "# nbin: " << nbin << endl;
  // // record
  // vector<vector<double > > record (nbin);
  // vector<vector<double > > record_j (nbin);
  // cout << "# j table is " << jfile << endl;
  // cout << "# angle names is " << ifile << endl;
  // cout << "# traj names  is " << tfile << endl;

  double metal, metau;
  metal = center - width;
  if (metal < -180) metal += 360.;
  metau = center + width;
  if (metau > 180) metau -= 360.;
  cout << "# metastable set: [ " << metal << " , " << metau << " ]" << endl;
  // metastable set
  MetastableSet ms (metal, metau);
  // average
  BlockAverage_acc ba (nDataInBlock);
  BlockAverage_acc ba_fht (nDataInBlock);
  
  // analyze traj
  ifstream fpname (idfile.c_str());
  if (!fpname){
    std::cerr << "cannot open file " << idfile << std::endl;
    return 1;
  }
  char nameline [MaxLineLength];
  int printCount = 0;
  int countFile = 0;
  int countFound = 0;
  int countUnFound = 0;
  
  while (fpname.getline(nameline, MaxLineLength)){
    if (nameline[0] == '#') continue;
    string filename (nameline);
    filename += string("/") + iffile;
    ifstream angname (filename.c_str());
    if (!angname){
      std::cerr << "cannot open file " << filename << std::endl;
      return 1;
    }
    if (printCount == 10) {
      printf ("# reading file %s       \r", nameline);
      printCount = 0;
    }
    printCount ++;
    countFile ++;
    double times = 0.;
    double anglev;
    char valueline [MaxLineLength];
    bool find = false;
    while (angname.getline(valueline, MaxLineLength)){
      if (valueline[0] == '#' || valueline[0] == '@') continue;
      vector<string > words;
      StringOperation::split (string(valueline), words);
      if (words.size() < 2) {
	cerr << "wrong file format of " << filename << endl;
	exit (1);
      }
      times  = (atof(words[0].c_str()));
      anglev = (atof(words[1].c_str()));
      if (ms.inSet(anglev)) {
	find = true;
	break;
      }
      if (times > gate) {
	find = false;
	break;
      }
    }
    if (find == true && times <= gate){
      ba.deposite (1.0);
      countFound ++ ;
    }
    else {
      ba.deposite (0.0);
      countUnFound ++;
    }
  }
  
  printf ("\n");
  printf ("# read %d files, %d ( %.1f %% ) hit meta, %d ( %.1f %% ) do not hit\n",
	  countFile,
	  countFound, ((double)(countFound))/((double)(countFile)) * 100.,
	  countUnFound, ((double)(countUnFound))/((double)(countFile)) * 100.);

  ba.calculate ();
  
  printf ("# prob. first hitting time smaller than %f\n", gate);
  printf ("%e   %e\n", ba.getAvg(), ba.getAvgError());
  

  //   double sum_j = 0.;
  //   if (ms.inSet(anglev[0])) {
  //     record[angleIdx(anglev[0], bin)].push_back (0.0);
  //     record_j[angleIdx(anglev[0], bin)].push_back (0.0);
  //     countFound ++;
  //   }
  //   else {
  //     bool found = true;
  //     unsigned nextstab = 0;
  //     do {
  // 	nextstab ++;
  // 	if (nextstab >= anglev.size()){
  // 	  found = false;
  // 	  break;
  // 	}
  //     } while (! ms.inSet(anglev[nextstab]));
  //     if (found) {
  // 	countFound ++;
  // 	record [angleIdx(anglev[0], bin)].push_back (-(times[0] - times[nextstab]));
  // 	sum_j = 0.;
  // 	for (unsigned ii = 0; ii < nextstab; ++ii){
  // 	  AngleCalculator ac (boxes[ii]);
  // 	  double tmp = - pp.value_periodic (anglev[ii]) ;
  // 	  double tmpAngle = ac.calAngle (c4Coord[ii]);
  // 	  vector<vector<double > > tmpGrad;
  // 	  ac.calGrad (c4Coord[ii], angleDx, tmpGrad);
  // 	  // vector<vector<double > > yy (4);
  // 	  // vector<double > yy1(3, 0.);
  // 	  // yy[3] = yy[2] = yy[1] = yy[0] = yy1;
  // 	  // yy[0][0] = 0.1;
  // 	  // yy[2][1] = 0.1;
  // 	  // yy[3][1] = yy[3][2] = 0.1;
  // 	  // ac.calGrad (yy, angleDx, tmpGrad);
  // 	  double sumGrad = 0.;
  // 	  for (unsigned kk = 0; kk < tmpGrad.size(); ++kk){
  // 	    for (unsigned ll = 0; ll < 3; ++ll){
  // 	      sumGrad += tmpGrad[kk][ll] * tmpGrad[kk][ll];
  // 	    }
  // 	  }
  // 	  // printf ("time %f, angle %f %f, sumGrad %e, dj is %f\n", times[ii], anglev[ii], tmpAngle, sumGrad, tmp);
  // 	  sum_j += 0.5 * tmp * tmp * sumGrad * (times[1] - times[0]);
  // 	}
  // 	record_j [angleIdx(anglev[0], bin)].push_back (sum_j);
  //     }
  //     else {
  // 	countUnFound ++;
  //     }
  //   }
  // }
  // printf ("\n");
  // printf ("# read %d files, %d ( %.1f %% ) hit meta, %d ( %.1f %% ) do not hit\n",
  // 	  countFile,
  // 	  countFound, ((double)(countFound))/((double)(countFile)) * 100.,
  // 	  countUnFound, ((double)(countUnFound))/((double)(countFile)) * 100.);

  // BlockAverage ba;  
  // BlockAverage ba_j;  
  // FILE * fout = fopen (ofile.c_str(), "w");
  // if (fout == NULL){
  //   cout << "cannot open file " << ofile<< endl;
  //   exit (1);
  // }
  // for (unsigned ii = 0; ii < record.size(); ++ii){
  //   double myangle = -180 + bin * (ii + 0.5);
  //   if (record[ii].size() == 0){
  //     int prevIdx = int(ii);
  //     do {
  // 	prevIdx --;
  //     } while (prevIdx >= 0 && record[prevIdx].size() == 0);
  //     if (prevIdx >= 0){
  // 	fprintf (fout, "%f %e %e\n",
  // 		 myangle,
  // 		 sigma * ba.getAvg() + xi * ba_j.getAvg(), sigma * ba.getAvgError() + xi * ba_j.getAvgError()
  // 		 );
  //     }
  //     else{
  // 	fprintf (fout, "%f\n", myangle);
  //     }
  //   }
  //   else if (record[ii].size() < nBlock * 5){
  //     ba.processData (record[ii], record[ii].size());
  //     ba_j.processData (record_j[ii], record_j[ii].size());
  //     fprintf (fout, "%f %e %e    %e %e    %e %e\n",
  // 	       myangle,
  // 	       sigma * ba.getAvg() + xi * ba_j.getAvg(), sigma * ba.getAvgError() + xi * ba_j.getAvgError(),
  // 	       ba.getAvg(), ba.getAvgError(),
  // 	       ba_j.getAvg(), ba_j.getAvgError()
  // 	  );
  //   }
  //   else {
  //     ba.processData (record[ii], nBlock);
  //     ba_j.processData (record_j[ii], nBlock);
  //     fprintf (fout, "%f %e %e    %e %e    %e %e\n",
  // 	       myangle,
  // 	       sigma * ba.getAvg() + xi * ba_j.getAvg(), sigma * ba.getAvgError() + xi * ba_j.getAvgError(),
  // 	       ba.getAvg(), ba.getAvgError(),
  // 	       ba_j.getAvg(), ba_j.getAvgError()
  // 	  );
  //     // fprintf (fout, "%f %e %e\n", myangle, ba.getAvg(), ba.getAvgError());
  //   }
  // }  
  // fclose (fout);

  // string ofile1 = ofile + string (".f");
  // fout = fopen (ofile.c_str(), "w");
  // if (fout == NULL){
  //   cout << "cannot open file " << ofile<< endl;
  //   exit (1);
  // }
  // for (unsigned ii = 0; ii < record.size(); ++ii){
  //   double myangle = -180 + bin * (ii + 0.5);
  //   if (record[ii].size() == 0){
  //     fprintf (fout, "%f\n", myangle);
  //   }
  //   else if (record[ii].size() < nBlock * 5){
  //     ba.processData (record[ii], record[ii].size());
  //     ba_j.processData (record_j[ii], record_j[ii].size());
  //     fprintf (fout, "%f %e %e    %e %e    %e %e\n",
  // 	       myangle,
  // 	       sigma * ba.getAvg() + xi * ba_j.getAvg(), sigma * ba.getAvgError() + xi * ba_j.getAvgError(),
  // 	       ba.getAvg(), ba.getAvgError(),
  // 	       ba_j.getAvg(), ba_j.getAvgError()
  // 	  );
  //   }
  //   else {
  //     ba.processData (record[ii], nBlock);
  //     ba_j.processData (record_j[ii], nBlock);
  //     fprintf (fout, "%f %e %e  \n",
  // 	       myangle,
  // 	       sigma * ba.getAvg() + xi * ba_j.getAvg(), sigma * ba.getAvgError() + xi * ba_j.getAvgError(),
  // 	       ba.getAvg(), ba.getAvgError(),
  // 	       ba_j.getAvg(), ba_j.getAvgError()
  // 	  );
  //     // fprintf (fout, "%f %e %e\n", myangle, ba.getAvg(), ba.getAvgError());
  //   }
  // }  
  // fclose (fout);
  
  
  return 0;
}
