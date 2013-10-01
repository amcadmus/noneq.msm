#include "Traj.h"

Traj::
Traj (const int & ndata)
    : data (ndata), posi(0), nValid(0)
{
}

Traj::
Traj () 
    : posi(0), nValid(0)
{
}

void Traj::
reinit (const int & ndata) 
{
  data.clear ();
  data.resize (ndata);
  posi = 0;
  nValid = 0;
}
