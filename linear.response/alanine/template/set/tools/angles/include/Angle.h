#ifndef __NONEQ_MD_Angle_wanghan_h__
#define __NONEQ_MD_Angle_wanghan_h__

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
#include "Defines.h"

using namespace std;

class AngleCalculator
{
  VectorType box;
public:
  AngleCalculator (const VectorType & box);
  void calPhiPsi (const vector<vector<ValueType> > & ala,
		  ValueType & phi,
		  ValueType & psi);
}
    ;

#endif
