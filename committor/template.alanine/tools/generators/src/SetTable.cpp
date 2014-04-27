#include "SetTable.h"

// void AngleSetTable1D::
// reinit (const string & filename)
// {
//   FILE * fp = fopen (filename.c_str(), "r");
//   if (fp == NULL){
//     cerr << "cannot open file " << filename << endl;
//     exit (1);
//   }

//   data.clear();
//   int tmp;
//   while (fscanf (fp, "%d", &tmp) == 1){
//     data.push_back (tmp);
//   }

//   bin = 360. / (double(data.size()));
// }

// int AngleSetTable1D::
// calIndex (const double & angle)
// {
//   int idx = (angle + 180) / bin;
//   if (idx < 0) idx += int(data.size());
//   if (idx >= int(data.size())) idx -= int(data.size());
//   return idx;
// }

// int AngleSetTable1D::
// calIndicate (const double & angle)
// {
//   return data[calIndex(angle)];
// }

// int AngleSetTable1D::
// calIndicate (const int & idx)
// {
//   return data[idx];
// }


#define MaxLineLength 2048

int AngleSetTable2D::
calIndex (const double & phi,
	  const double & psi)
{
  int idxx = (phi + 180) / binx;
  int idxy = (psi + 180) / biny;
  if (idxx < 0)   idxx += nx;
  if (idxx >= nx) idxx -= nx;
  if (idxy < 0)   idxy += ny;
  if (idxy >= ny) idxy -= ny;
  
  return idxx * ny + idxy;
}

int AngleSetTable2D::
calIndicate (const double & phi,
	     const double & psi)
{
  return data[calIndex(phi, psi)];
}

int AngleSetTable2D::
calIndicate (const int & idx)
{
  return data[idx];
}


void AngleSetTable2D::
reinit (const string & filename)
{
  ifstream fpname0 (filename.c_str());  
  if (!fpname0){
    cerr << "cannot open file " << filename << endl;
    exit (1);
  }

  data.clear();
  nx = 0;
  ny = 0;
  
  char nameline[MaxLineLength];

  while (fpname0.getline(nameline, MaxLineLength)){
    if (nameline[0] == '#' || nameline[0] == '@') continue;
    vector<string > words;
    StringOperation::split (string(nameline), words);
    if (ny == 0){
      ny = words.size();
    }
    else {
      if (ny != int(words.size())){
	cerr << "lines of " << filename << " does not have the same number of indicators, format wrong, exit" << endl;
	exit(1);
      }
    }
    for (unsigned ii = 0; ii < words.size(); ++ii){
      data.push_back (atoi(words[ii].c_str()));
    }
    nx ++;
  }
  
  binx = 360. / (double(nx));
  biny = 360. / (double(ny));
}


void AngleSetTable2D::
print (const string & filename)
{
  FILE * fp = fopen (filename.c_str(), "w");
  for (int ii = 0; ii < nx; ++ii){
    for (int jj = 0; jj < ny; ++jj){
      int idx = ii * ny + jj;
      if (idx >= 0){
      fprintf (fp, "  %d", data[idx]);
      }
      else {
      fprintf (fp, " %d", data[idx]);
      }
    }
    fprintf (fp, "\n");
  }
  fclose (fp);
}

void AngleSetTable2D::
plus   (const AngleSetTable2D & table)
{
  if (nx != table.nx || ny != table.ny){
    cerr << "tables are not match, cannot add" << endl;
    return ;
  }

  for (int ii = 0; ii < nx; ++ii){
    for (int jj = 0; jj < ny; ++jj){
      int idx = ii * ny + jj;
      data[idx] += table.data[idx];
    }
  }
}
