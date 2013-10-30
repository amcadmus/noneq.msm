#include "Traj.h"

Traj::
Traj (const int & ndata)
    : data (ndata, Dofs()), posi(0), nValid(0)
{
}
