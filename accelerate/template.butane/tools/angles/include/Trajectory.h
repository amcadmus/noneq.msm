#ifndef __MDFileManager_Trajectory_h_wanghan__
#define __MDFileManager_Trajectory_h_wanghan__

// #include "Defines.h"
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_xtc.h"
#include "xdrfile/xdrfile_trr.h"
#include <vector> 

using namespace std;

class XtcLoader 
{
  XDRFILE *xd;
  int natoms;
  int step;
  float time;
  vector<double > box;
  rvec * xx;
  float prec;
  bool inited;
  void clear ();
public:
  XtcLoader ();
  ~XtcLoader ();
  XtcLoader (const char * filename);
  bool reinit (const char * filename);
  bool load ();
public:
  const vector<double > & getBox () const {return box;}
  float getTime () const {return time;}
public:
  void getFrame (vector<vector<double > > & frame);
}
    ;


class TrrLoader 
{
  XDRFILE *xd;
  int natoms;
  int step;
  float time;
  float lambda;
  vector<double > box;
  rvec * xx;
  rvec * vv;
  rvec * ff;
  bool inited;
  void clear ();
public:
  TrrLoader ();
  ~TrrLoader ();
  TrrLoader (const char * filename);
  bool reinit (const char * filename);
  bool load ();
public:
  const vector<double > & getBox () const {return box;}
  float getTime () const {return time;}
public:
  void getFrame (vector<vector<double > > & xx,
		 vector<vector<double > > & vv,
		 vector<vector<double > > & ff);
}
    ;

#endif
