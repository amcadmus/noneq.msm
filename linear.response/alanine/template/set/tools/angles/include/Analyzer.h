#ifndef __Analyzer_h_wanghan__
#define __Analyzer_h_wanghan__

#include "Defines.h"
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_xtc.h"
#include "xdrfile/xdrfile_trr.h"
#include <vector>

using namespace std;

struct TopInfo
{
  unsigned numAtomOnALA;
  unsigned comIndexALA;
  unsigned numAtomOnH2o;
  unsigned OIndexH2o;
  unsigned H1IndexH2o;
  unsigned H2IndexH2o;
  TopInfo ();
}
    ;

class TrajLoader_xtc 
{
  XDRFILE *xd;
  int natoms;
  unsigned numMolALA, numMolH2o;
  unsigned numAtomALA, numAtomH2o;
  int step;
  float time;
  VectorType box;
  rvec * xx;
  float prec;
  bool inited;
  void clear ();
  TopInfo tinfo;
public:
  TrajLoader_xtc ();
  ~TrajLoader_xtc ();
  TrajLoader_xtc (const char * filename,
		  const TopInfo info = TopInfo());
  void reinit (const char * filename,
	       const TopInfo info = TopInfo());
  bool load ();
public:
  // int getNumAtom () const {return natoms;}
  // int getNumAtomALA () const {return numAtomALA;}
  // int getNumAtomH2o () const {return numAtomH2o;}
  VectorType getBox () const {return box;}
  float getTime () const {return time;}
  int getNatoms () const {return natoms;}
public:
  void formCoords (std::vector<std::vector<ValueType > > & ala,
		   std::vector<std::vector<ValueType > > & h2o);
}
    ;


// class TrajLoader_trr
// {
//   XDRFILE *xd;
//   int natoms;
//   unsigned numMolALA, numMolH2o;
//   unsigned numAtomALA, numAtomH2o;
//   int step;
//   float time;
//   VectorType box;
//   rvec * xx;
//   rvec * vv;
//   rvec * ff;
//   float lambda;
//   bool inited;
//   void clear ();
//   TopInfo tinfo;
// public:
//   TrajLoader_trr ();
//   ~TrajLoader_trr ();
//   TrajLoader_trr (const char * filename,
// 		  const TopInfo info = TopInfo());
//   void reinit (const char * filename,
// 	       const TopInfo info = TopInfo());
//   bool load ();
// public:
//   // int getNumAtom () const {return natoms;}
//   // int getNumAtomALA () const {return numAtomALA;}
//   // int getNumAtomH2o () const {return numAtomH2o;}
//   VectorType getBox () const {return box;}
//   float getTime () const {return time;}
// public:
//   void formCoords (std::vector<std::vector<ValueType > > & alaxx,
// 		   std::vector<std::vector<ValueType > > & alavv,
// 		   std::vector<std::vector<ValueType > > & alaff,
// 		   std::vector<std::vector<ValueType > > & h2oxx,
// 		   std::vector<std::vector<ValueType > > & h2ovv,
// 		   std::vector<std::vector<ValueType > > & h2off);
// }
//     ;



#endif
