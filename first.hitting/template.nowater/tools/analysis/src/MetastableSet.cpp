#include "MetastableSet.h"

MetastableSet::
MetastableSet (const double & phb,
	       const double & phe)
    : phi_b (phb), phi_e (phe)
{
}

static bool
metastableSet_inrange (const double & bb,
		       const double & ee,
		       const double & vv)
{
  if (bb < ee){
    return (vv >= bb - 1e-8 && vv <= ee + 1e-8);
  }
  else {
    return ( (vv >= bb - 1e-8 && vv <= 180) || (vv >= -180 && vv <= ee + 1e-8) );
  }
}


bool MetastableSet::
inSet (const double & phi) const
{
  return ( metastableSet_inrange(phi_b, phi_e, phi) );
}

