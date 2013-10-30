#ifndef __Defines_h_NONEQMSM_wanghan__
#define __Defines_h_NONEQMSM_wanghan__

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

#define NUMDOFS	1

using namespace std;

struct Dofs
{
  vector<double > xx;
  vector<double > vv;
  Dofs ();
}
    ;

inline Dofs::
Dofs ()
    : xx(NUMDOFS), vv(NUMDOFS)
{
}


#endif
