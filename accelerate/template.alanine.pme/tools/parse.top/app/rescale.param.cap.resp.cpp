#include "GmxTop.h"
#include "GmxType.h"
#include <cmath>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;
using namespace GmxTop;

// return in unit s
// static double
// cal_period_bond (const double & mass,
// 		 const double & kk)
// {
//   if (kk != 0){
//     return 2. * M_PI * sqrt(mass / (kk * 1e24));
//   }
//   else {
//     return 0;
//   }
// }

static double
cal_period_bond2 (const double & mass0,
		  const double & mass1,
		  const double & kk)
{
  double mass = mass0 * mass1 / (mass0 + mass1);
  if (kk != 0){
    return 2. * M_PI * sqrt(mass / (kk * 1e24));
  }
  else {
    return 0;
  }
}

// return in unit s
// static double
// cal_period_angle (const double & mass,
// 		  const double & ll,
// 		  const double & kk)
// {
//   double II = mass * ll * ll;
//   if (kk != 0){
//     return 2. * M_PI * sqrt(II / (kk * 1e24));
//   }
//   else {
//     return 0;
//   }
// }

static double
cal_period_angle2 (const double & mass0,
		   const double & ll0,
		   const double & mass1,
		   const double & ll1,
		   const double & kk)
{
  double II0 = mass0 * ll0 * ll0;
  double II1 = mass1 * ll1 * ll1;
  double II = II0 * II1 / (II0 + II1);
  if (kk != 0){
    return 2. * M_PI * sqrt(II / (kk * 1e24));
  }
  else {
    return 0;
  }
}


int main (int argc, char **argv) {

  std::string ifile, ofile, opfile;
  double sbond, sangle, sdihedral, scmap;
  double capTbond, capTangle; 

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")      
      ("cap-period-bond", po::value<double > (&capTbond)->default_value (8.631997e-15), "cap period of the system, unit: s")
      ("cap-period-angle", po::value<double > (&capTangle)->default_value (6.955613e-15), "cap period of the system, unit: s")
      ("scal-bond,b", po::value<double > (&sbond)->default_value (1.0), "bond scale")
      ("scal-angle,a", po::value<double > (&sangle)->default_value (1.0), "angle scale")
      ("scal-dihedral,d", po::value<double > (&sdihedral)->default_value (1.0), "dihedral scale")
      ("scal-cmap,c", po::value<double > (&scmap)->default_value (1.0), "cmap scale")
      ("output,o", po::value<std::string > (&ofile)->default_value ("out.top"), "the output top")
      ("input,f",  po::value<std::string > (&ifile)->default_value ("topol.top"), "the input top");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  GmxTop::gmx_sys_top systop;
  GmxTop::parseTop (ifile.c_str(), systop);
  GmxTop::gmx_sys_types systype;
  GmxTop::parseType (ifile.c_str(), systype);

  // loop over all molecules
  for (unsigned mol_idx = 0; mol_idx < systop.moles.size(); ++mol_idx){
    // loop over all bonds of molecule moles[mole_idx]
    for (unsigned ii = 0; ii < systop.moles[mol_idx].bonds.size(); ++ii){
      vector<double > params = systop.moles[mol_idx].bonds[ii].params;
      int funct = systop.moles[mol_idx].bonds[ii].funct;
      int iiidx = systop.moles[mol_idx].bonds[ii].atom_idx_ii();
      int jjidx = systop.moles[mol_idx].bonds[ii].atom_idx_jj();
      gmx_atom iiatom = systop.moles[mol_idx].atoms[iiidx];
      gmx_atom jjatom = systop.moles[mol_idx].atoms[jjidx];
      if (params.size() == 0){
	gmx_bondtypes_item tmpbondtype;
	if (! matchBondType (iiatom.at_type, jjatom.at_type, systype, tmpbondtype)){
	  cerr << "cannot find bond between type: " << iiatom.at_type << " " << jjatom.at_type << endl;
	  return 1;
	}
	params = tmpbondtype.params;
      }
      if (funct != 1){
	cerr << "unexpected funct of the bond interaction " << funct << endl;
	return 1;
      }
      if (params.size() < 2){
	cerr << "wrong params size" << endl;
	return 1;
      }

      if (iiatom.mass == 0){
	gmx_atomtypes_item tmpatomtype;
	matchAtomType (iiatom.at_type, systype, tmpatomtype);
	iiatom.mass = tmpatomtype.mass;
      }
      if (jjatom.mass == 0){
	gmx_atomtypes_item tmpatomtype;
	matchAtomType (jjatom.at_type, systype, tmpatomtype);
	jjatom.mass = tmpatomtype.mass;
      }
      
      // double period0 = cal_period_bond (iiatom.mass, params[1]);
      // double period1 = cal_period_bond (jjatom.mass, params[1]);
      double period2 = cal_period_bond2 (iiatom.mass, jjatom.mass, params[1]);

      double periodScale = period2 / capTbond;
      periodScale = periodScale * periodScale;
      double realScale = 0.;
      double diffScale = 0.;
      if (periodScale < sbond) {
	realScale = periodScale;
	diffScale = sbond - periodScale;
      }
      else {
	realScale = sbond;
	diffScale = 0.;
      }
      if (systop.moles[mol_idx].name != string("SOL")){
	params.resize(4);
	params[2] = (params[0]);
	params[3] = (params[1] * diffScale);      
	params[1] *= realScale;
      }
      else {
	params.resize(4);
	params[2] = (params[0]);
	params[3] = (params[1] * 0.);	// do not calculate the response from the water bonds
	params[1] *= realScale;
      }
      systop.moles[mol_idx].bonds[ii].params = params;

      double period2new = cal_period_bond2 (iiatom.mass, jjatom.mass, params[1]);
      
      printf ("find bond between %s and %s, mass: %f %f, type %d, params: ",
	      iiatom.at_name.c_str(), jjatom.at_name.c_str(),
	      iiatom.mass, jjatom.mass,
	      funct);
      for (unsigned kk = 0; kk < params.size(); ++kk){
	printf (" %e ", params[kk]);
      }
      printf ("   periods: %f -> %f fs", period2 * 1e15, period2new * 1e15);
      cout << endl;
    }
    cout << endl;    


    for (unsigned ii = 0; ii < systop.moles[mol_idx].angles.size(); ++ii){
      vector<double > params = systop.moles[mol_idx].angles[ii].params;
      int funct = systop.moles[mol_idx].angles[ii].funct;
      int iiidx = systop.moles[mol_idx].angles[ii].atom_idx_ii();
      int jjidx = systop.moles[mol_idx].angles[ii].atom_idx_jj();
      int kkidx = systop.moles[mol_idx].angles[ii].atom_idx_kk();
      GmxTop::gmx_atom iiatom = systop.moles[mol_idx].atoms[iiidx];
      GmxTop::gmx_atom jjatom = systop.moles[mol_idx].atoms[jjidx];
      GmxTop::gmx_atom kkatom = systop.moles[mol_idx].atoms[kkidx];
      if (params.size() == 0){      
	GmxTop::gmx_angletypes_item tmpangletype;
	if (! matchAngleType (iiatom.at_type, jjatom.at_type, kkatom.at_type, systype, tmpangletype)){
	  cerr << "cannot find angle between type: " << iiatom.at_type
	       << " " << jjatom.at_type 
	       << " " << kkatom.at_type << endl;
	  return 1;
	}
	params = tmpangletype.params;
      }
      if (funct != 5 && funct != 1){
	cerr << "unexpected funct of the angle interaction " << funct << endl;
	return 1;
      }

      if (iiatom.mass == 0){
	gmx_atomtypes_item tmpatomtype;
	matchAtomType (iiatom.at_type, systype, tmpatomtype);
	iiatom.mass = tmpatomtype.mass;
      }
      if (jjatom.mass == 0){
	gmx_atomtypes_item tmpatomtype;
	matchAtomType (jjatom.at_type, systype, tmpatomtype);
	jjatom.mass = tmpatomtype.mass;
      }
      if (kkatom.mass == 0){
	gmx_atomtypes_item tmpatomtype;
	matchAtomType (kkatom.at_type, systype, tmpatomtype);
	kkatom.mass = tmpatomtype.mass;
      }

      gmx_bonds_item tmpbond;
      vector<double > params0, params1;
      int funct0, funct1;
      if (! matchBond (iiidx, jjidx, systop.moles[mol_idx], tmpbond)){
	cerr << "cannot find bond between type: " << iiatom.at_type << " " << jjatom.at_type << endl;
	return 1;
      }
      params0 = tmpbond.params;
      funct0 = tmpbond.funct;
      if (params0.size() == 0){
	gmx_bondtypes_item tmpbondtype;
	if (! matchBondType (iiatom.at_type, jjatom.at_type, systype, tmpbondtype)){
	  cerr << "cannot find bond between type: " << iiatom.at_type << " " << jjatom.at_type << endl;
	  return 1;
	}
	params0 = tmpbondtype.params;
      }
      if (! matchBond (kkidx, jjidx, systop.moles[mol_idx], tmpbond)){
      	cerr << "cannot find bond between type: " << kkatom.at_type << " " << jjatom.at_type << endl;
      	return 1;
      }
      params1 = tmpbond.params;
      funct1 = tmpbond.funct;
      if (params1.size() == 0){
	gmx_bondtypes_item tmpbondtype;
	if (! matchBondType (kkatom.at_type, jjatom.at_type, systype, tmpbondtype)){
	  cerr << "cannot find bond between type: " << kkatom.at_type << " " << jjatom.at_type << endl;
	  return 1;
	}
	params1 = tmpbondtype.params;
      }
      
      if (funct0 != 1){
	cerr << "unexpected funct of the bond interaction " << funct0 << endl;
      }
      if (funct1 != 1){
	cerr << "unexpected funct of the bond interaction " << funct1 << endl;
      }
      if (params0.size() < 2 || params1.size() < 2){
	cerr << "wrong params size" << endl;
	return 1;
      }

      // double period0 = cal_period_angle (iiatom.mass, params0[0], params[1]);
      // double period1 = cal_period_angle (kkatom.mass, params1[0], params[1]);
      double period2 = cal_period_angle2 (iiatom.mass, params0[0], kkatom.mass, params1[0], params[1]);
      double periodmin = period2;
      double period012 = 0.;
      if (funct == 5){	
	period012 = cal_period_bond2 (iiatom.mass, kkatom.mass, params[3]);
	if (period012 != 0){
	  if (period012 < periodmin) periodmin = period012;
	}
      }
      
      double periodScale = periodmin / capTangle;
      periodScale = periodScale * periodScale;
      double realScale = 0.;
      double diffScale = 0.;
      if (periodScale < sangle) {
	realScale = periodScale;
	diffScale = sangle - periodScale;
      }
      else {
	realScale = sangle;
	diffScale = 0.;
      }
      if (funct == 5){
	params[1] *= realScale;
	params[3] *= realScale;
      }
      else {
	params[1] *= realScale;
      }
      systop.moles[mol_idx].angles[ii].params = params;      

      double period2new = cal_period_angle2 (iiatom.mass, params0[0], kkatom.mass, params1[0], params[1]);

      // do not calculate the response from water!
      if (funct == 1 && systop.moles[mol_idx].name == string("SOL")){
	  systop.moles[mol_idx].angles[ii].params.resize(4);
	  systop.moles[mol_idx].angles[ii].params[2] = systop.moles[mol_idx].angles[ii].params[0];
	  systop.moles[mol_idx].angles[ii].params[3] = 0.;
      }
      if (funct == 5){
	// do not split the angle type 5, because will cause nb neighbor exclusion problem!
	if (systop.moles[mol_idx].name != string("SOL")){
	  systop.moles[mol_idx].angles[ii].params.resize(6);
	  systop.moles[mol_idx].angles[ii].params[4] = systop.moles[mol_idx].angles[ii].params[1] / realScale * diffScale;
	  systop.moles[mol_idx].angles[ii].params[5] = systop.moles[mol_idx].angles[ii].params[3] / realScale * diffScale;
	}
	// double period01a = cal_period_bond (iiatom.mass, params[3]);
	// double period01b = cal_period_bond (kkatom.mass, params[3]);
	double period012new = cal_period_bond2 (iiatom.mass, kkatom.mass, params[3]);
	printf ("find angle between %s %s and %s, type %d, params: ",
		iiatom.at_name.c_str(), jjatom.at_name.c_str(), kkatom.at_name.c_str(), funct);
	for (unsigned kk = 0; kk < params.size(); ++kk){
		printf (" %e ", params[kk]);
	}
	printf ("   periods %f -> %f    %f -> %f fs",
		period2 * 1e15, period2new * 1e15,
		period012 * 1e15, period012new * 1e15);
      }
      else {
	printf ("find angle between %s %s and %s, type %d, params: ",
		iiatom.at_name.c_str(), jjatom.at_name.c_str(), kkatom.at_name.c_str(), funct);
	for (unsigned kk = 0; kk < params.size(); ++kk){
		printf (" %e ", params[kk]);
	}
	printf ("   periods %f -> %f fs",
		period2 * 1e15, period2new * 1e15);
      }
      
      
      cout << endl;
    }
    cout << endl;
  }
  

  // if (systype.bondtypes.size() > 0){
  //   for (unsigned ii = 0; ii < systype.bondtypes.size(); ++ii){
  //     if (systype.bondtypes[ii].funct == 1 ||
  // 	  systype.bondtypes[ii].funct == 2 ||
  // 	  systype.bondtypes[ii].funct == 6 ||
  // 	  systype.bondtypes[ii].funct == 7 ||
  // 	  systype.bondtypes[ii].funct == 8 ||
  // 	  systype.bondtypes[ii].funct == 9 ){
  // 	systype.bondtypes[ii].params[1] *= sbond;
  //     }
  //     else {
  // 	cerr << "bond type is not 1,2,6,7,8,9, do nothing!" << endl;
  //     }
  //   }
  // }

  // if (systype.angletypes.size() > 0){
  //   for (unsigned ii = 0; ii < systype.angletypes.size(); ++ii){
  //     if (systype.angletypes[ii].funct == 1 ||
  // 	  systype.angletypes[ii].funct == 2 ||
  // 	  systype.angletypes[ii].funct == 8 ){
  // 	systype.angletypes[ii].params[1] *= sangle;
  //     }
  //     else if (systype.angletypes[ii].funct == 5 ){
  // 	systype.angletypes[ii].params[1] *= sangle;
  // 	systype.angletypes[ii].params[3] *= sanglebond;
  //     }
  //     else {
  // 	cerr << "angle type is not 1,2,5,8, do nothing!" << endl;
  //     }
  //   }
  // }
  if (systype.dihedraltypes.size() > 0){
    for (unsigned ii = 0; ii < systype.dihedraltypes.size(); ++ii){
      if (systype.dihedraltypes[ii].funct == 1 ||
  	  systype.dihedraltypes[ii].funct == 2 ||
  	  systype.dihedraltypes[ii].funct == 4 ||
  	  systype.dihedraltypes[ii].funct == 8 ||
  	  systype.dihedraltypes[ii].funct == 9 ){
  	systype.dihedraltypes[ii].params[1] *= sdihedral;
      }
      else {
  	cerr << "dihedral type is not 1,2,4,8,9, do nothing!" << endl;
      }
    }
  }
  if (systype.cmaptypes.size() > 0){
    for (unsigned ii = 0; ii < systype.cmaptypes.size(); ++ii){
      if (systype.cmaptypes[ii].funct == 1 ){
  	for (unsigned jj = 0; jj < systype.cmaptypes[ii].params.size(); ++jj){
  	  systype.cmaptypes[ii].params[jj] *= scmap;
  	}
      }
      else {
  	cerr << "cmap type is not 1, do nothing!" << endl;
      }
    }
  }


  systype.bondtypes.clear();
  systype.angletypes.clear();
  
  for (unsigned ii = 0; ii < systop.moles.size(); ++ii){
    // for (unsigned jj = 0; jj < systop.moles[ii].bonds.size(); ++jj){
    //   if (systop.moles[ii].bonds[jj].funct == 1 ||
    // 	  systop.moles[ii].bonds[jj].funct == 2 ||
    // 	  systop.moles[ii].bonds[jj].funct == 6 ||
    // 	  systop.moles[ii].bonds[jj].funct == 7 ||
    // 	  systop.moles[ii].bonds[jj].funct == 8 ||
    // 	  systop.moles[ii].bonds[jj].funct == 9 ){
    // 	if (systop.moles[ii].bonds[jj].params.size() >=2 ){
    // 	  systop.moles[ii].bonds[jj].params[1] *= sbond;
    // 	}
    //   }
    // }
    // for (unsigned jj = 0; jj < systop.moles[ii].angles.size(); ++jj){
    //   if (systop.moles[ii].angles[jj].funct == 1 ||
    // 	  systop.moles[ii].angles[jj].funct == 2 ||
    // 	  systop.moles[ii].angles[jj].funct == 8 ){
    // 	if (systop.moles[ii].angles[jj].params.size() >=2 ){
    // 	  systop.moles[ii].angles[jj].params[1] *= sangle;
    // 	}
    //   }
    // }
    for (unsigned jj = 0; jj < systop.moles[ii].dihedrals.size(); ++jj){
      if (systop.moles[ii].dihedrals[jj].funct == 1 ||
  	  systop.moles[ii].dihedrals[jj].funct == 2 ||
  	  systop.moles[ii].dihedrals[jj].funct == 4 ||
  	  systop.moles[ii].dihedrals[jj].funct == 8 ||
  	  systop.moles[ii].dihedrals[jj].funct == 9 ){
  	if (systop.moles[ii].dihedrals[jj].params.size() >=2 ){
  	  systop.moles[ii].dihedrals[jj].params[1] *= sdihedral;
  	}
      }
    }
  }  
  
  FILE * fp = fopen (ofile.c_str(), "w");
  systype.print (fp);
  systop.print (fp);

  return 0;
}
