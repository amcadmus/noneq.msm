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

void
read_base_info (const string & ibfile,
		vector<double > & kk)
{
  double tmp;
  kk.clear(); 
  FILE * fp = fopen (ibfile.c_str(), "r");
  if (fp == NULL){
    cerr << "cannot open file " << ibfile << endl;
    exit(1);
  }
  while (fscanf (fp, "%lf", &tmp) == 1){
    kk.push_back (tmp);
  }
  fclose (fp);
}


int main(int argc, char * argv[])
{
  std::string iffile, igxsfile, idfile, ibfile, ofile, icfile;
  unsigned nDataInBlock;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("num-in-block,n",  po::value<unsigned > (&nDataInBlock)->default_value (1),   "number of data in a block")
      ("input-coreset,c", po::value<std::string > (&icfile)->default_value ("coreset.dat"), "the input base info")
      ("input-base-info,b", po::value<std::string > (&ibfile)->default_value ("base.info"), "the input base info")
      ("input-dir,d", po::value<std::string > (&idfile)->default_value ("success.dir.name"), "file including successful dir names")
      ("input-angle-name,a", po::value<std::string > (&iffile)->default_value ("angaver.xvg"), "the angle file name")
      ("input-gxs-name,x", po::value<std::string > (&igxsfile)->default_value ("gxs12.out"), "the gxs file name")
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

  int nBase = -1;
  vector<double > baseKK;
  read_base_info (ibfile, baseKK);
  nBase = int(baseKK.size());

  AngleSetTable1D stable;
  stable.reinit (icfile);

  // average
  BlockAverage_acc ba (nDataInBlock);
  vector<double > gxs1;
  vector<vector<double > > gxs2;
  vector<BlockAverage_acc > vecb;
  vector<vector<BlockAverage_acc > > matA;
  gxs1.resize (nBase);
  gxs2.resize (nBase);
  vecb.resize (nBase);
  matA.resize (nBase);
  for (int ii = 0; ii < nBase; ++ii){
    gxs2[ii].resize (nBase);
    matA[ii].resize (nBase);
  }
  
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
  int countNumInGate = 0;
  
  while (fpname.getline(nameline, MaxLineLength)){
    if (nameline[0] == '#') continue;
    string filename_ang (nameline);
    filename_ang += string("/") + iffile;
    countFile ++;
    ifstream angname (filename_ang.c_str());
    if (!angname){
      ba.deposite (0.0);
      for (int ii = 0; ii < nBase; ++ii){
	vecb[ii].deposite (0.0);
	for (int jj = 0; jj < nBase; ++jj){
	  matA[ii][jj].deposite(0.0);
	}
      }
      countUnFound ++;      
    }
    else {
      string filename_gxs (nameline);
      filename_gxs += string("/") + igxsfile;
      ifstream gxsname (filename_gxs.c_str());
      if (!gxsname){
	std::cerr << "\n cannot open file " << filename_gxs << std::endl;
	return 1;
      }
      if (printCount == 100) {
	// printf ("# reading file %s and %s      \r", filename_ang.c_str(), filename_gxs.c_str());
	// fflush (stdout);
	printCount = 0;
      }
      printCount ++;
      double times;
      double anglev;
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
	times  = (atof(words[0].c_str()));
	anglev = (atof(words[1].c_str()));
	StringOperation::split (string(valueline_gxs), words);
	for (int ii = 0; ii < nBase; ++ii){
	  gxs1[ii] = atof (words[1+ii].c_str());
	  for (int jj = 0; jj < nBase; ++jj){
	    gxs2[ii][jj] = atof (words[1 + nBase + ii * nBase + jj].c_str());
	  }
	}

	if (stable.calIndicate(anglev) == 1) {
	  find = true;
	  break;
	}
      }
    
      if (find){
	double sum1 = 0.;
	double sum2 = 0.;
	for (int ii = 0; ii < nBase; ++ii){
	  sum1 += gxs1[ii] * baseKK[ii];
	  for (int jj = 0; jj < nBase; ++jj){
	    sum2 += gxs2[ii][jj] * baseKK[ii] * baseKK[jj];
	  }
	}
	double tmpexp = exp( - sum1 - 0.5 * sum2);
	countNumInGate ++;
	ba.deposite (1.0 * tmpexp);
	printf ("file %s time %f deposited: %e\n", filename_gxs.c_str(), times, tmpexp);
	for (int ii = 0; ii < nBase; ++ii){
	  vecb[ii].deposite (1.0 * gxs1[ii] * tmpexp);
	  for (int jj = 0; jj < nBase; ++jj){
	    matA[ii][jj].deposite(1.0 * gxs2[ii][jj]  * tmpexp);
	  }
	}
	countFound ++ ;
      }
      else {
	ba.deposite (0.0);
	for (int ii = 0; ii < nBase; ++ii){
	  vecb[ii].deposite (0.0);
	  for (int jj = 0; jj < nBase; ++jj){
	    matA[ii][jj].deposite(0.0);
	  }
	}
	countUnFound ++;
      }
    }
  }
  
  printf ("\n");
  printf ("# read %d files, %d ( %.1f %% ) hit meta, %d ( %.1f %% ) do not hit\n",
	  countFile,
	  countFound, ((double)(countFound))/((double)(countFile)) * 100.,
	  countUnFound, ((double)(countUnFound))/((double)(countFile)) * 100.);
  // printf ("# read %d files, %d ( %.1f %% ) hit meta in gate\n",
  // 	  countFile,
  // 	  countNumInGate, ((double)(countNumInGate))/((double)(countFile)) * 100.);

  ba.calculate ();
  
  printf ("# committor\n");
  printf ("# value   var    error\n");
  printf ("%e   %e   %e\n", ba.getAvg(), ba.getVar(), ba.getAvgError());

  printf ("# vect b\n");
  for (int ii = 0; ii < nBase; ++ii){
    vecb[ii].calculate();
    printf ("%e   ", vecb[ii].getAvg());
  }
  printf ("\n# vect b error\n");
  for (int ii = 0; ii < nBase; ++ii){
    printf ("%e   ", vecb[ii].getAvgError());
  }
  printf ("\n#;\n");
  printf ("# mat A\n");
  for (int ii = 0; ii < nBase; ++ii){
    for (int jj = 0; jj < nBase; ++jj){
      matA[ii][jj].calculate();
      printf ("%e   ", matA[ii][jj].getAvg());
    }
    printf (";\n");
  }
  printf ("# mat A error\n");
  for (int ii = 0; ii < nBase; ++ii){
    for (int jj = 0; jj < nBase; ++jj){
      matA[ii][jj].calculate();
      printf ("%e   ", matA[ii][jj].getAvgError());
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
