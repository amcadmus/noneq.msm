#ifndef __NONEQ_MD_Traj_h_wanghan__
#define __NONEQ_MD_Traj_h_wanghan__

#include "Defines.h"
#include <vector>

using namespace std;

class Traj
{
  vector<vector<double > > data;
  int posi;
  int nValid;
  int normalizedIdx (const int & idx) const;
public :
  Traj ();
  Traj (const int & ndata);
public:
  void push_back (const vector<double > & value);
  void reinit (const int & ndata);
  void clear ();
  bool full () const {return nValid == int(data.size());}
  const int & getNValid () const {return nValid;}
  const int tailPosi () const {return normalizedIdx(posi-1);}
  const vector<double > & getValue (const int & idx) const {return data[normalizedIdx(idx)];}
  const vector<double > & front () const {return getValue(posi);}
  const vector<double > & back  () const {return getValue(posi-1);}
}
    ;

inline int Traj::
normalizedIdx (const int & idx) const
{
  int newidx (idx);
  if (idx >= int(data.size())) newidx -= data.size();
  else if (idx < 0) newidx += data.size();
  return newidx;
}

inline void Traj::
push_back (const vector<double > & value)
{
  data[posi] = value;
  posi++;
  posi = normalizedIdx (posi);
  if (nValid < int(data.size())) ++nValid;
}

inline void Traj::
clear ()
{
  posi = 0;
  nValid = 0;
}


#endif
