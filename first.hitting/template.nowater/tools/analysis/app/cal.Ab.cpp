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

int cal_n_base (const int in)
{
  int value = 0;
  bool find = false;
  if (in == 0) return 0;
  while ((value + value * value) != in && value < 100000){
    value ++;
    find = true;
  }
  if (find){
    return value;
  }
  else {
    return -1;
  }
}


int main(int argc, char * argv[])
{
  std::string iffile, igxsfile, idfile, ofile;
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
      ("input-dir,d", po::value<std::string > (&idfile)->default_value ("success.dir.name"), "file including successful dir names")
      ("input-angle-name", po::value<std::string > (&iffile)->default_value ("angaver.xvg"), "the angle file name")
      ("input-gxs-name", po::value<std::string > (&igxsfile)->default_value ("gxs12.out"), "the gxs file name")
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
  vector<BlockAverage_acc > vecb;
  vector<vector<BlockAverage_acc > > matA;
  
  // analyze traj
  ifstream fpname (idfile.c_str());
  if (!fpname){
    std::cerr << "\n cannot open file " << idfile << std::endl;
    return 1;
  }
  char nameline [MaxLineLength];
  int printCount = 0;
  int countFile = 0;
  int countFound = 0;
  int countUnFound = 0;
  int nBase = -1;
  
  while (fpname.getline(nameline, MaxLineLength)){
    if (nameline[0] == '#') continue;
    string filename_ang (nameline);
    filename_ang += string("/") + iffile;
    ifstream angname (filename_ang.c_str());
    if (!angname){
      std::cerr << "\n cannot open file " << filename_ang << std::endl;
      return 1;
    }
    string filename_gxs (nameline);
    filename_gxs += string("/") + igxsfile;
    ifstream gxsname (filename_gxs.c_str());
    if (!gxsname){
      std::cerr << "\n cannot open file " << filename_gxs << std::endl;
      return 1;
    }
    if (printCount == 100) {
      printf ("# reading file %s and %s      \r", filename_ang.c_str(), filename_gxs.c_str());
      fflush (stdout);
      printCount = 0;
    }
    printCount ++;
    countFile ++;
    vector<double > times;
    vector<double > anglev;
    vector<vector<double > > gxs1;
    vector<vector<vector<double > > > gxs2;
    char valueline_ang [MaxLineLength];
    char valueline_gxs [MaxLineLength];
    bool find = false;
    while (angname.getline(valueline_ang, MaxLineLength) &&
	   gxsname.getline(valueline_gxs, MaxLineLength)){
      if (valueline_ang[0] == '#' || valueline_ang[0] == '@' ||
	  valueline_gxs[0] == '#' || valueline_gxs[0] == '@') {
	cerr << "data files should not contain any line starting with # or @\n" << endl;
	return 1;
      }
      vector<string > words;
      StringOperation::split (string(valueline_ang), words);
      if (words.size() < 2) {
	cerr << "wrong file format of " << filename_ang << endl;
	exit (1);
      }
      times .push_back (atof(words[0].c_str()));
      anglev.push_back (atof(words[1].c_str()));
      StringOperation::split (string(valueline_gxs), words);
      if (nBase < 0) {
	nBase = cal_n_base (words.size() - 1);
	if (nBase == -1){
	  cerr << "invalid input line of " << filename_gxs << ", may be more than 100000 bases? " << endl;
	  return 1;
	}
	vecb.resize (nBase);
	matA.resize (nBase);
	for (int ii = 0; ii < nBase; ++ii){
	  matA[ii].resize (nBase);
	}
      }
      vector<double > tmp1 (nBase, 0.);
      vector<vector<double > > tmp2(nBase, vector<double >(nBase, 0.)) ;
      for (int ii = 0; ii < nBase; ++ii){
	tmp1[ii] = atof (words[1+ii].c_str());
	for (int jj = 0; jj < nBase; ++jj){
	  tmp2[ii][jj] = atof (words[1 + nBase + ii * nBase + jj].c_str());
	}
      }
      gxs1.push_back (tmp1);
      gxs2.push_back (tmp2);

      if (ms.inSet(anglev.back())) {
	find = true;
	break;
      }
    }
    if (find == true){
      double sum1 = 0.;
      double sum2 = 0.;
      for (int ii = 0; ii < nBase; ++ii){
	sum1 += gxs1.back()[ii];
	for (int jj = 0; jj < nBase; ++jj){
	  sum2 += gxs2.back()[ii][jj];
	}
      }
      if (times.back() < gate){
	ba.deposite (1.0 * exp( - sum1 - 0.5 * sum2));
	// ba.deposite (1.0 * (1. - sum1));
	// ba.deposite (1.0 * exp( - sum1));
	// ba.deposite (1.0);
	for (int ii = 0; ii < nBase; ++ii){
	  vecb[ii].deposite (1.0 * gxs1.back()[ii] * exp( - sum1 - 0.5 * sum2));
	  for (int jj = 0; jj < nBase; ++jj){
	    matA[ii][jj].deposite(1.0 * gxs2.back()[ii][jj]  * exp( - sum1 - 0.5 * sum2));
	  }
	}
      }
      else {
	ba.deposite (0.0);
	for (int ii = 0; ii < nBase; ++ii){
	  vecb[ii].deposite (0.0);
	  for (int jj = 0; jj < nBase; ++jj){
	    matA[ii][jj].deposite(0.0);
	  }
	}
      }
      ba_fht.deposite (times.back() * exp( - sum1 - 0.5 * sum2));
      // ba_fht.deposite (times.back() * (1. - sum1));
      // ba_fht.deposite (times.back() * exp( - sum1));
      // ba_fht.deposite (times.back());
      countFound ++ ;
    }
    else {
      countUnFound ++;
    }
  }
  
  printf ("\n");
  printf ("# read %d files, %d ( %.1f %% ) hit meta, %d ( %.1f %% ) do not hit\n",
	  countFile,
	  countFound, ((double)(countFound))/((double)(countFile)) * 100.,
	  countUnFound, ((double)(countUnFound))/((double)(countFile)) * 100.);

  ba.calculate ();
  ba_fht.calculate ();
  
  printf ("# avg. first hitting time\n");
  printf ("%f   %f\n", ba_fht.getAvg(), ba_fht.getAvgError());
  printf ("# prob. first hitting time smaller than %f\n", gate);
  printf ("%e   %e\n", ba.getAvg(), ba.getAvgError());

  printf ("# vect b\n");
  for (int ii = 0; ii < nBase; ++ii){
    vecb[ii].calculate();
    printf ("%e   ", vecb[ii].getAvg());
  }
  printf ("# vect b error\n");
  for (int ii = 0; ii < nBase; ++ii){
    printf ("%e   ", vecb[ii].getAvgError());
  }
  printf ("#;\n");
  printf ("# mat A error\n");
  for (int ii = 0; ii < nBase; ++ii){
    for (int jj = 0; jj < nBase; ++jj){
      matA[ii][jj].calculate();
      printf ("%e   ", matA[ii][jj].getAvg());
    }
    printf (";\n");
  }


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
