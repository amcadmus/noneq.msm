#include "Analyzer.h"
#include <string.h>
#include <fstream>
#include <iostream>
#include <algorithm>

TopInfo::
TopInfo ()
    : numAtomOnALA (22),
      comIndexALA (0),
      numAtomOnH2o (3),
      OIndexH2o (0),
      H1IndexH2o (1),
      H2IndexH2o (2)
{
}

void TrajLoader_xtc::
clear ()
{
  if (inited){
    free (xx);
    xdrfile_close (xd);
    inited = false;
  }
}


TrajLoader_xtc::
TrajLoader_xtc ()
    : inited (false)
{
}

TrajLoader_xtc::
~TrajLoader_xtc ()
{
  clear ();
}

TrajLoader_xtc::
TrajLoader_xtc (const char * filename,
		const TopInfo info)
    : inited (false)
{
  reinit (filename, info);
}  

void TrajLoader_xtc::
reinit (const char * filename,
	const TopInfo info)
{
  char tmpname[2048];
  strncpy (tmpname, filename, 2047);
  
  xd = xdrfile_open (filename, "r");
  if (xd == NULL){
    std::cerr << "cannot open file " << filename << std::endl;
    return ;
  }
  read_xtc_natoms (tmpname, &natoms);
  step = 0;
  time = 0.;
  box = VectorType (0.,0.,0.);
  xx = (rvec *) malloc (sizeof(rvec) * natoms);
  prec = 1000.;

  tinfo = info;
  numMolALA = 1;
  numAtomALA = tinfo.numAtomOnALA * numMolALA;
  numAtomH2o = natoms - numAtomALA;
  if (numAtomH2o % tinfo.numAtomOnH2o != 0){
    std::cerr << "inconsistent num atom and mol size, exit" << std::endl;
    exit(1);
  }
  numMolH2o = numAtomH2o / tinfo.numAtomOnH2o;
  
  inited = true;

  load ();

  xdrfile_close (xd);
  xd = xdrfile_open (filename, "r");  
  if (xd == NULL){
    std::cerr << "cannot open file " << filename << std::endl;
    return ;
  }
}

bool TrajLoader_xtc::
load ()
{
  if (inited){
    matrix tmpBox;
    int st = read_xtc (xd, natoms, &step, &time, tmpBox, xx, &prec);
    box.x = tmpBox[0][0];
    box.y = tmpBox[1][1];
    box.z = tmpBox[2][2];
    if (st == exdrOK) return true;
    else return false;
  }
  else {
    return false;
  }
}

void TrajLoader_xtc::
formCoords (std::vector<std::vector<ValueType > > & ala,
	    std::vector<std::vector<ValueType > > & h2o)
{
  
  ala.resize(numMolALA * tinfo.numAtomOnALA);
  h2o.resize(numMolH2o * 3);

  for (unsigned ii = 0; ii < numMolALA; ++ii){
    for (unsigned jj = 0; jj < tinfo.numAtomOnALA; ++jj){
      ala[ii * tinfo.numAtomOnALA + jj].resize(3);
      for (unsigned dd = 0; dd < 3; ++dd){
	ala[ii * tinfo.numAtomOnALA + jj][dd] = xx[ii * tinfo.numAtomOnALA + jj][dd];
      }
    }
  }
  for (unsigned ii = 0; ii < numMolH2o; ++ii){
    h2o[ii*3+0].resize(3);
    h2o[ii*3+1].resize(3);
    h2o[ii*3+2].resize(3);
    if (int(numAtomALA + ii * tinfo.numAtomOnH2o + tinfo.H2IndexH2o) > int(natoms)){
      std::cerr << "wrong index of water! exit" << std::endl;
      exit(1);
    }
    for (unsigned dd = 0; dd < 3; ++dd){
      h2o[ii*3+0][dd] = xx[numAtomALA + ii * tinfo.numAtomOnH2o + tinfo.OIndexH2o][dd];
      h2o[ii*3+1][dd] = xx[numAtomALA + ii * tinfo.numAtomOnH2o + tinfo.H1IndexH2o][dd];
      h2o[ii*3+2][dd] = xx[numAtomALA + ii * tinfo.numAtomOnH2o + tinfo.H2IndexH2o][dd];
    }
  }
}


// void TrajLoader_trr::
// clear ()
// {
//   if (inited){
//     free (xx);
//     free (vv);
//     free (ff);
//     xdrfile_close (xd);
//     inited = false;
//   }
// }


// TrajLoader_trr::
// TrajLoader_trr ()
//     : inited (false)
// {
// }

// TrajLoader_trr::
// ~TrajLoader_trr ()
// {
//   clear ();
// }

// TrajLoader_trr::
// TrajLoader_trr (const char * filename,
// 		const TopInfo info)
//     : inited (false)
// {
//   reinit (filename, info);
// }  

// void TrajLoader_trr::
// reinit (const char * filename,
// 	const TopInfo info)
// {
//   char tmpname[2048];
//   strncpy (tmpname, filename, 2047);
  
//   xd = xdrfile_open (filename, "r");
//   if (xd == NULL){
//     std::cerr << "cannot open file " << filename << std::endl;
//     return ;
//   }
//   read_trr_natoms (tmpname, &natoms);
//   step = 0;
//   time = 0.;
//   box = VectorType (0.,0.,0.);
//   xx = (rvec *) malloc (sizeof(rvec) * natoms);
//   vv = (rvec *) malloc (sizeof(rvec) * natoms);
//   ff = (rvec *) malloc (sizeof(rvec) * natoms);
//   lambda = 1;
  
//   tinfo = info;
//   numMolALA = 1;
//   numAtomALA = tinfo.numAtomOnALA * numMolALA;
//   numAtomH2o = natoms - numAtomALA;
//   if (numAtomH2o % tinfo.numAtomOnH2o != 0){
//     std::cerr << "inconsistent num atom and mol size, exit" << std::endl;
//     exit(1);
//   }
//   numMolH2o = numAtomH2o / tinfo.numAtomOnH2o;
  
//   inited = true;

//   load ();

//   xdrfile_close (xd);
//   xd = xdrfile_open (filename, "r");  
//   if (xd == NULL){
//     std::cerr << "cannot open file " << filename << std::endl;
//     return ;
//   }
// }

// bool TrajLoader_trr::
// load ()
// {
//   if (inited){
//     matrix tmpBox;
//     int st = read_trr (xd, natoms, &step, &time, &lambda, tmpBox, xx, vv, ff);
//     box.x = tmpBox[0][0];
//     box.y = tmpBox[1][1];
//     box.z = tmpBox[2][2];
//     if (st == exdrOK) return true;
//     else return false;
//   }
//   else {
//     return false;
//   }
// }

// void TrajLoader_trr::
// formCoords (std::vector<std::vector<ValueType > > & alaxx,
// 	    std::vector<std::vector<ValueType > > & alavv,
// 	    std::vector<std::vector<ValueType > > & alaff,
// 	    std::vector<std::vector<ValueType > > & h2oxx,
// 	    std::vector<std::vector<ValueType > > & h2ovv,
// 	    std::vector<std::vector<ValueType > > & h2off)
// {
  
//   alaxx.resize(numMolALA);
//   alavv.resize(numMolALA);
//   alaff.resize(numMolALA);
//   h2oxx.resize(numMolH2o * 3);
//   h2ovv.resize(numMolH2o * 3);
//   h2off.resize(numMolH2o * 3);

//   for (unsigned ii = 0; ii < numMolALA; ++ii){
//     alaxx[ii].resize(3);
//     alavv[ii].resize(3);
//     alaff[ii].resize(3);
//     for (unsigned dd = 0; dd < 3; ++dd){
//       alaxx[ii][dd] = xx[ii * tinfo.numAtomOnALA + tinfo.comIndexALA][dd];
//       alavv[ii][dd] = vv[ii * tinfo.numAtomOnALA + tinfo.comIndexALA][dd];
//       alaff[ii][dd] = ff[ii * tinfo.numAtomOnALA + tinfo.comIndexALA][dd];
//     }
//   }
//   for (unsigned ii = 0; ii < numMolH2o; ++ii){
//     h2oxx[ii*3+0].resize(3);
//     h2oxx[ii*3+1].resize(3);
//     h2oxx[ii*3+2].resize(3);
//     h2ovv[ii*3+0].resize(3);
//     h2ovv[ii*3+1].resize(3);
//     h2ovv[ii*3+2].resize(3);
//     h2off[ii*3+0].resize(3);
//     h2off[ii*3+1].resize(3);
//     h2off[ii*3+2].resize(3);
//     if (int(numAtomALA + ii * tinfo.numAtomOnH2o + tinfo.H2IndexH2o) > int(natoms)){
//       std::cerr << "wrong index of water! exit" << std::endl;
//       exit(1);
//     }
//     for (unsigned dd = 0; dd < 3; ++dd){
//       h2oxx[ii*3+0][dd] = xx[numAtomALA + ii * tinfo.numAtomOnH2o + tinfo.OIndexH2o][dd];
//       h2oxx[ii*3+1][dd] = xx[numAtomALA + ii * tinfo.numAtomOnH2o + tinfo.H1IndexH2o][dd];
//       h2oxx[ii*3+2][dd] = xx[numAtomALA + ii * tinfo.numAtomOnH2o + tinfo.H2IndexH2o][dd];
//       h2ovv[ii*3+0][dd] = vv[numAtomALA + ii * tinfo.numAtomOnH2o + tinfo.OIndexH2o][dd];
//       h2ovv[ii*3+1][dd] = vv[numAtomALA + ii * tinfo.numAtomOnH2o + tinfo.H1IndexH2o][dd];
//       h2ovv[ii*3+2][dd] = vv[numAtomALA + ii * tinfo.numAtomOnH2o + tinfo.H2IndexH2o][dd];
//       h2off[ii*3+0][dd] = ff[numAtomALA + ii * tinfo.numAtomOnH2o + tinfo.OIndexH2o][dd];
//       h2off[ii*3+1][dd] = ff[numAtomALA + ii * tinfo.numAtomOnH2o + tinfo.H1IndexH2o][dd];
//       h2off[ii*3+2][dd] = ff[numAtomALA + ii * tinfo.numAtomOnH2o + tinfo.H2IndexH2o][dd];
//     }
//   }
// }



