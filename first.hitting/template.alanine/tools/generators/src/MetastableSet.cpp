#include "MetastableSet.h"


MetastableSet::
MetastableSet (const double & phb,
               const double & phe,
               const double & psb,
               const double & pse)
    : phi_b (phb), phi_e (phe), psi_b (psb), psi_e (pse)
{
}

static bool
metastableSet_inrange (const double & bb,
                       const double & ee,
                       const double & vv)
{
  if (bb < ee){
    return (vv >= bb && vv < ee);
  }
  else {
    return ( (vv >= bb && vv <= 180) || (vv >= -180 && vv < ee) );
  }
}


bool MetastableSet::
inSet (const double & phi,
       const double & psi) const
{
  return ( metastableSet_inrange(phi_b, phi_e, phi) &&
           metastableSet_inrange(psi_b, psi_e, psi) );
}

