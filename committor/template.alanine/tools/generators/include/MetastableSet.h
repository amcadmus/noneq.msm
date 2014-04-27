#ifndef __MetastableSet_h_wanghan__
#define __MetastableSet_h_wanghan__

class MetastableSet 
{
private:
  double phi_b, phi_e;
  double psi_b, psi_e;
public:
  MetastableSet (const double & phb,
                 const double & phe,
                 const double & psb,
                 const double & pse);
  bool inSet (const double & phi,
              const double & psi) const;
}
    ;

#endif

