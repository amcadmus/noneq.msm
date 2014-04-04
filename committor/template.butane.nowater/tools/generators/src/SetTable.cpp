#include "SetTable.h"

void AngleSetTable1D::
reinit (const string & filename)
{
  FILE * fp = fopen (filename.c_str(), "r");
  if (fp == NULL){
    cerr << "cannot open file " << filename << endl;
    exit (1);
  }

  data.clear();
  int tmp;
  while (fscanf (fp, "%d", &tmp) == 1){
    data.push_back (tmp);
  }

  bin = 360. / (double(data.size()));
}

int AngleSetTable1D::
calIndex (const double & angle)
{
  int idx = (angle + 180) / bin;
  if (idx < 0) idx += int(data.size());
  if (idx >= int(data.size())) idx -= int(data.size());
  return idx;
}

int AngleSetTable1D::
calIndicate (const double & angle)
{
  return data[calIndex(angle)];
}

int AngleSetTable1D::
calIndicate (const int & idx)
{
  return data[idx];
}

