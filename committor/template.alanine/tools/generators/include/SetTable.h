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
#include "StringSplit.h"

using namespace std;

// angle range: -180 -> 180
// class AngleSetTable1D
// {
//   vector<int > data;
//   double bin;
// public:
//   void reinit (const string & file);
//   int calIndex (const double & angle);
//   int calIndicate (const double & angle);
//   int calIndicate (const int & index);
//   int size () const {return int(data.size());}
// }
//     ;

class AngleSetTable2D
{
  int nx, ny;
  vector<int > data;
  double binx, biny;
public:
  void reinit (const string & file);
  void print  (const string & file);
  void plus   (const AngleSetTable2D & table);
  int calIndex (const double & phi,
		const double & psi);
  int calIndicate (const double & phi,
		   const double & psi);
  int calIndicate (const int & index);
  int sizePhi () const {return nx;}
  int sizePsi () const {return ny;}
}
    ;

#endif
