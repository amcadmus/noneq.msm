#include "GmxTop.h"
#include "GmxType.h"
#include <cmath>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;
using namespace GmxTop;

// return in unit s
static double
cal_period_bond (const double & mass,
		 const double & kk)
{
  if (kk != 0){
    return 2. * M_PI * sqrt(mass / (kk * 1e24));
  }
  else {
    return 0;
  }
}

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
static double
cal_period_angle (const double & mass,
		  const double & ll,
		  const double & kk)
{
  double II = mass * ll * ll;
  if (kk != 0){
    return 2. * M_PI * sqrt(II / (kk * 1e24));
  }
  else {
    return 0;
  }
}

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

int regulate_idx (const int idx,
		  const int nn)
{
  int myidx (idx);
  if (myidx < 0) myidx += nn;
  else if (myidx >= nn) myidx -= nn;
  return myidx;
}


void read_angle (const string & filename_o,
		 const string & filename_r,
		 gmx_cmaptypes_item & cmap)
{
  vector<double > data_o, data_r;
  vector<double > data_p;
  vector<bool > data_valid;
  FILE * fp;
  double xx, yy, vv, ee;
  double kbT = 2.5;
  double eps = 1e-12;
  // double max = -kbT * log(eps);

  printf ("size of cmap param is %d\n", cmap.params.size());
  
  fp = fopen (filename_o.c_str(), "r");
  if (fp == NULL){
    cerr << "cannot open file " << filename_o << endl;
    exit (1);
  }
  while (4 == fscanf(fp, "%lf%lf%lf%lf", &xx, &yy, &vv, &ee)){
    data_o.push_back (vv);
  }
  fclose (fp);
  fp = fopen (filename_r.c_str(), "r");
  if (fp == NULL){
    cerr << "cannot open file " << filename_r << endl;
    exit (1);
  }
  while (4 == fscanf(fp, "%lf%lf%lf%lf", &xx, &yy, &vv, &ee)){
    data_r.push_back (vv);
  }
  fclose (fp);

  if (int(data_o.size()) != cmap.ngrid0 * cmap.ngrid1 ||
      int(data_r.size()) != cmap.ngrid0 * cmap.ngrid1 ){
    cerr << "inconsistent size of the grids" << endl;
    exit (1);
  }

  data_p.resize (cmap.ngrid0 * cmap.ngrid1 );
  data_valid.resize (cmap.ngrid0 * cmap.ngrid1 );
  
  for (int ii = 0; ii < cmap.ngrid0; ++ii){
    for (int jj = 0; jj < cmap.ngrid1; ++jj){
      int idx = ii * cmap.ngrid1 + jj;
      if (data_o[idx] < eps || data_r[idx] < eps){
	data_valid[idx] = false;
      }
      else {
	data_valid[idx] = true;
      }
    }
  }

  FILE * fp1 = fopen ("dpangle.out", "w");
  double hx = 360. / double(cmap.ngrid0);
  for (int ii = 0; ii < cmap.ngrid0; ++ii){
    for (int jj = 0; jj < cmap.ngrid1; ++jj){
      int idx = ii * cmap.ngrid1 + jj;
      if (data_valid[idx]){
	double pot_o, pot_r;
	pot_o = -kbT * log (data_o[idx]);
	pot_r = -kbT * log (data_r[idx]);
	data_p[idx] = pot_r - pot_o;
      }
      else {
	data_p[idx] = 0.;
      }
      fprintf (fp1, "%f %f %f\n", ii * hx, jj * hx, data_p[idx]);
    }
    fprintf (fp1, "\n");
  }
  fclose(fp1);

  // double hi = 360. / cmap.ngrid0;
  // double hj = 360. / cmap.ngrid1;
  for (int ii = 0; ii < cmap.ngrid0; ++ii){
    // int reg_iip1 = regulate_idx(ii+1, cmap.ngrid0);
    // int reg_iim1 = regulate_idx(ii-1, cmap.ngrid0);    
    for (int jj = 0; jj < cmap.ngrid1; ++jj){
      // int reg_jjp1 = regulate_idx(jj+1, cmap.ngrid1);
      // int reg_jjm1 = regulate_idx(jj-1, cmap.ngrid1);
      int idx = ii * cmap.ngrid1 + jj;
      // int idxjp1 = ii * cmap.ngrid1 + reg_jjp1;
      // int idxjm1 = ii * cmap.ngrid1 + reg_jjm1;
      // int idxip1 = reg_iip1 * cmap.ngrid1 + jj;
      // int idxim1 = reg_iim1 * cmap.ngrid1 + jj;
      // int idxip1jp1 = reg_iip1 * cmap.ngrid1 + reg_jjp1;
      // int idxim1jm1 = reg_iim1 * cmap.ngrid1 + reg_jjm1;
      // int idxip1jm1 = reg_iip1 * cmap.ngrid1 + reg_jjm1;
      // int idxim1jp1 = reg_iim1 * cmap.ngrid1 + reg_jjp1;
      // double dvdi = (cmap.params[idxip1*4] - cmap.params[idxim1*4]) / hi * 0.5;
      // double dvdj = (cmap.params[idxjp1*4] - cmap.params[idxjm1*4]) / hj * 0.5;
      // double d2vdidj = (cmap.params[idxip1jp1*4] + cmap.params[idxim1jm1*4] -
      // 			cmap.params[idxip1jm1*4] + cmap.params[idxim1jp1*4] ) / (hi * hj * 4);
      // printf ("%f %f %f %f   \t   %f %f %f\n",
      // 	      cmap.params[idx*4+0],
      // 	      cmap.params[idx*4+1], 
      // 	      cmap.params[idx*4+2], 
      // 	      cmap.params[idx*4+3], 
      // 	      dvdi,
      // 	      dvdj,
      // 	      d2vdidj);
      cmap.params[idx] += data_p[idx];
    }
  }
}


int main (int argc, char **argv) {

  std::string ifile, ofile, opfile, anglefile1, anglefile2;
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
      ("avg-angle-1", po::value<string > (&anglefile1)->default_value ("avg.angle.out.1"), "angle distrib file 1")
      ("avg-angle-2", po::value<string > (&anglefile2)->default_value ("avg.angle.out.8"), "angle distrib file 2")
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
      double realScale = 0;
      if (periodScale < sbond) {
	realScale = periodScale;
      }
      else {
	realScale = sbond;
      }
      params[1] *= realScale;
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
      double realScale = 0;
      if (periodScale < sangle) {
	realScale = periodScale;
      }
      else {
	realScale = sangle;
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
      
      if (funct == 5){
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

    for (unsigned ii = 0; ii < systop.moles[mol_idx].cmap.size(); ++ii){
      int iiidx = systop.moles[mol_idx].cmap[ii].atom_idx_ii();
      int jjidx = systop.moles[mol_idx].cmap[ii].atom_idx_jj();
      int kkidx = systop.moles[mol_idx].cmap[ii].atom_idx_kk();
      int llidx = systop.moles[mol_idx].cmap[ii].atom_idx_ll();
      int mmidx = systop.moles[mol_idx].cmap[ii].atom_idx_mm();
      gmx_atom iiatom = systop.moles[mol_idx].atoms[iiidx];
      gmx_atom jjatom = systop.moles[mol_idx].atoms[jjidx];
      gmx_atom kkatom = systop.moles[mol_idx].atoms[kkidx];
      gmx_atom llatom = systop.moles[mol_idx].atoms[llidx];
      gmx_atom mmatom = systop.moles[mol_idx].atoms[mmidx];
      gmx_cmaptypes_item tmpcmaptype;
      int tmpidx;
      if (! matchCMAPType (iiatom.at_type,
			   jjatom.at_type,
			   kkatom.at_type,
			   llatom.at_type,
			   mmatom.at_type,
			   systype,
			   tmpcmaptype,
			   tmpidx)){
	cerr << "cannot find cmap type " << endl;
	return 1;
      }
      read_angle (anglefile1, anglefile2, systype.cmaptypes[tmpidx]);
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
