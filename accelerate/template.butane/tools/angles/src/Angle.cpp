#include "Angle.h"

static void
normalize (vector<vector<ValueType> > & ala,
	   const VectorType & box)
{
  // if (ala.size() <= 1){
  //   return ;
  // }
  for (unsigned ii = 1; ii < ala.size() ; ++ii){
    if      (ala[ii][0] - ala[0][0] > 0.5 * box.x) ala[ii][0] -= box.x;
    else if (ala[ii][0] - ala[0][0] <-0.5 * box.x) ala[ii][0] += box.x;
    if      (ala[ii][1] - ala[0][1] > 0.5 * box.y) ala[ii][1] -= box.y;
    else if (ala[ii][1] - ala[0][1] <-0.5 * box.y) ala[ii][1] += box.y;
    if      (ala[ii][2] - ala[0][2] > 0.5 * box.z) ala[ii][2] -= box.z;
    else if (ala[ii][2] - ala[0][2] <-0.5 * box.z) ala[ii][2] += box.z;
  }
}

static vector<ValueType>
crossProd (const vector<ValueType> & diff0,
	   const vector<ValueType> & diff1)
{
  vector<ValueType> re (3);
  re[0] = (diff0[1] * diff1[2] - diff0[2] * diff1[1]);
  re[1] =-(diff0[0] * diff1[2] - diff0[2] * diff1[0]);
  re[2] = (diff0[0] * diff1[1] - diff0[1] * diff1[0]);  
  return re;
}

static ValueType
dotProd (const vector<ValueType> & diff0,
	 const vector<ValueType> & diff1)
{
  return (diff0[0] * diff1[0] +
	  diff0[1] * diff1[1] +
	  diff0[2] * diff1[2] );
}

static vector<ValueType>
calNormVec (vector<vector<ValueType> > & group)
{
  vector<ValueType> diff0 (3), diff1(3);

  for (unsigned dd = 0; dd < 3; ++dd){
    diff0[dd] = group[0][dd] - group[1][dd];
    diff1[dd] = group[2][dd] - group[1][dd];
  }

  return crossProd (diff0, diff1);
}

static ValueType
calAngle (const vector<ValueType> aa,
	  const vector<ValueType> bb)
{
  ValueType ep = 1e-8;
  ValueType la = aa[0] * aa[0] + aa[1] * aa[1] + aa[2] * aa[2];
  ValueType lb = bb[0] * bb[0] + bb[1] * bb[1] + bb[2] * bb[2];
  la = sqrt (la);
  lb = sqrt (lb);
  if (la < ep || lb < ep){
    return -1000;
  }
  ValueType cosv = dotProd(aa, bb) / la / lb;
  if (cosv > 1) {
    cerr << "cos value " << cosv << " >  1" << endl << endl;;
    cosv = 1;
  }
  else if (cosv < -1){
    cerr << "cos value " << cosv << " < -1" << endl << endl;    
    cosv = -1;
  }
  return acos (cosv) / M_PI * 180;
}

AngleCalculator::
AngleCalculator (const VectorType & box_)
    : box (box_)
{
}


void AngleCalculator::
calPhiPsi (const vector<vector<ValueType> > & ala,
	   ValueType & phi,
	   ValueType & psi)
{
  vector<vector<ValueType> > tmpala(ala);
  normalize (tmpala, box);

  // angle psi: 4, 6, 8, 14
  // angle phi: 6, 8, 14, 16
  if (ala.size() != 22){
    cerr << "size of alanine is not 22! may be a wrong molecule" << endl;
    exit (1);
  }
  vector<ValueType> aa(3), bb(3), cc(3), bond(3);
  vector<vector<ValueType> > group(3);
  
  group[0] = ala[4];
  group[1] = ala[6];
  group[2] = ala[8];
  aa = calNormVec (group);
  group[0] = ala[6];
  group[1] = ala[8];
  group[2] = ala[14];
  bb = calNormVec (group);
  
  psi = calAngle (aa, bb);
  if (psi < -999){
    cerr << "wrong psi detected!" << endl;
  }
  cc = crossProd (aa, bb);
  for (unsigned dd = 0; dd < 3; ++dd){
    bond[dd] = ala[8][dd] -  ala[6][dd];
  }
  if (dotProd (cc, bond) < 0) psi = -psi;


  group[0] = ala[6];
  group[1] = ala[8];
  group[2] = ala[14];
  aa = calNormVec (group);
  group[0] = ala[8];
  group[1] = ala[14];
  group[2] = ala[16];
  bb = calNormVec (group);
  
  phi = calAngle (aa, bb);
  cc = crossProd (aa, bb);
  for (unsigned dd = 0; dd < 3; ++dd){
    bond[dd] = ala[14][dd] -  ala[8][dd];
  }
  if (dotProd (cc, bond) < 0) phi = -phi;
  
}


