#ifndef __Traj_h_NONEQMSM_wanghan__
#define __Traj_h_NONEQMSM_wanghan__

#include "Defines.h"

using namespace std;

class Traj
{
  vector<Dofs > data;
  int posi;
  int nValid;
  int normalizedIdx (const int & idx) const;
public :
  Traj (const int & ndata);
public:
  void push_back (const Dofs & value);
  void clear ();
  bool full () const {return nValid == int(data.size());}
  const int & getNValid () const {return nValid;}
  const int tailPosi () const {return normalizedIdx(posi-1);}
  const Dofs & getValue (const int & idx) const {return data[normalizedIdx(idx)];}
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
push_back (const Dofs & value)
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
