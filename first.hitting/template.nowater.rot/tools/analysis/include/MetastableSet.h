#ifndef __MetastableSet_h_wanghan__
#define __MetastableSet_h_wanghan__

#include <vector>
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

using namespace std;

class MetastableSet 
{
private:
  double phi_b, phi_e;
public:
  MetastableSet (const double & phb,
		 const double & phe);
  bool inSet (const double & phi) const;
}
    ;

#endif

