#ifndef __SetTable_h_wanghan__
#define __SetTable_h_wanghan__

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

// angle range: -180 -> 180
class AngleSetTable1D
{
  vector<int > data;
  double bin;
public:
  void reinit (const string & file);
  int calIndex (const double & angle);
  int calIndicate (const double & angle);
  int calIndicate (const int & index);
  int size () const {return int(data.size());}
}
    ;


#endif
