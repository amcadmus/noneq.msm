#include "Trajectory.h"
#include <stdlib.h>
#include <string.h>
#include <iostream>

void XtcLoader::
clear ()
{
  if (inited){
    free (xx);
    xdrfile_close (xd);
    time = 0.;
    step = 0;
    inited = false;
  }
}

XtcLoader::
XtcLoader ()
    : time (0.),
      box (vector<double > (3, 0.)),
      inited (false)
{
}

XtcLoader::
~XtcLoader ()
{
  clear ();
}

XtcLoader::
XtcLoader (const char * filename)
    : time (0.),
      box (vector<double > (3, 0.)),
      inited (false)
{
  reinit (filename);
}  

bool XtcLoader::
reinit (const char * filename)
{
  char tmpname[2048];
  strncpy (tmpname, filename, 2047);
  
  xd = xdrfile_open (filename, "r");
  if (xd == NULL){
    std::cerr << "cannot open file " << filename << std::endl;
    return false;
  }
  read_xtc_natoms (tmpname, &natoms);
  step = 0;
  time = 0.;
  box = vector<double > (3, 0.);
  xx = (rvec *) malloc (sizeof(rvec) * natoms);
  prec = 1000.;

  inited = true;

  load ();

  xdrfile_close (xd);
  xd = xdrfile_open (filename, "r");  
  if (xd == NULL){
    std::cerr << "cannot open file " << filename << std::endl;
    clear ();
    return false;
  }
  
  return true;
}

bool XtcLoader::
load ()
{
  if (inited){
    matrix tmpBox;
    int st = read_xtc (xd, natoms, &step, &time, tmpBox, xx, &prec);
    box[0] = tmpBox[0][0];
    box[1] = tmpBox[1][1];
    box[2] = tmpBox[2][2];
    if (st == exdrOK) return true;
    else return false;
  }
  else {
    std::cerr << "not initiated, do nothing." << std::endl;
    return false;
  }
}

void XtcLoader::
getFrame (vector<vector<double > > & frame)
{
  frame.resize (natoms);
  for (int ii = 0; ii < natoms; ++ii){
    frame[ii].resize(3);
    frame[ii][0] = xx[ii][0];
    frame[ii][1] = xx[ii][1];
    frame[ii][2] = xx[ii][2];
  }
}




void TrrLoader::
clear ()
{
  if (inited){
    free (xx);
    free (vv);
    free (ff);
    xdrfile_close (xd);
    time = 0.;
    step = 0;
    inited = false;
  }
}

TrrLoader::
TrrLoader ()
    : time (0.),
      box (vector<double > (3, 0.)),
      inited (false)
{
}

TrrLoader::
~TrrLoader ()
{
  clear ();
}

TrrLoader::
TrrLoader (const char * filename)
    : time (0.),
      box (vector<double > (3, 0.)),
      inited (false)
{
  reinit (filename);
}  

bool TrrLoader::
reinit (const char * filename)
{
  char tmpname[2048];
  strncpy (tmpname, filename, 2047);
  
  xd = xdrfile_open (filename, "r");
  if (xd == NULL){
    std::cerr << "cannot open file " << filename << std::endl;
    return false;
  }
  read_trr_natoms (tmpname, &natoms);
  step = 0;
  time = 0.;
  box = vector<double > (3, 0.);
  xx = (rvec *) malloc (sizeof(rvec) * natoms);
  vv = (rvec *) malloc (sizeof(rvec) * natoms);
  ff = (rvec *) malloc (sizeof(rvec) * natoms);

  inited = true;

  if (load ()) {
    xdrfile_close (xd);
    xd = xdrfile_open (filename, "r");  
    if (xd == NULL){
      std::cerr << "cannot open file " << filename << std::endl;
      clear ();
      return false;
    }
  }
  else {
    std::cerr << "fail to load the 1st frame" << endl;
    return false;
  }
  
  return true;
}

bool TrrLoader::
load ()
{
  if (inited){
    matrix tmpBox;
    int st = read_trr (xd, natoms, &step, &time, &lambda, tmpBox, xx, vv, ff);
    box[0] = tmpBox[0][0];
    box[1] = tmpBox[1][1];
    box[2] = tmpBox[2][2];
    if (st == exdrOK) return true;
    else return false;
  }
  else {
    std::cerr << "not initiated, do nothing." << std::endl;
    return false;
  }
}

void TrrLoader::
getFrame (vector<vector<double > > & xx_,
	  vector<vector<double > > & vv_,
	  vector<vector<double > > & ff_)
{
  xx_.resize (natoms);
  vv_.resize (natoms);
  ff_.resize (natoms);
  for (int ii = 0; ii < natoms; ++ii){
    xx_[ii].resize(3);
    xx_[ii][0] = xx[ii][0];
    xx_[ii][1] = xx[ii][1];
    xx_[ii][2] = xx[ii][2];
    vv_[ii].resize(3);
    vv_[ii][0] = vv[ii][0];
    vv_[ii][1] = vv[ii][1];
    vv_[ii][2] = vv[ii][2];
    ff_[ii].resize(3);
    ff_[ii][0] = ff[ii][0];
    ff_[ii][1] = ff[ii][1];
    ff_[ii][2] = ff[ii][2];
  }
}


