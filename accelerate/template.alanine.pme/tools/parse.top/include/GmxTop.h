#ifndef __GmxTop_h_wanghan__
#define __GmxTop_h_wanghan__

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include "GmxType.h"

#define MAX_LINE_LENGTH 2048
#define MAX_NAME_LENGTH 32

using namespace std;

namespace GmxTop {
  struct gmx_atom
  {
    int			id;
    string		at_type;
    int			res_nr;
    string		res_name;
    string		at_name;
    int			cgnr;
    double		charge;
    double		mass;
    gmx_atom		();
    void		clear ();    
    void		print (FILE * fp) const;
  }
      ;

  struct gmx_pairs_item
  {
    int			ii;
    int			jj;
    int			funct;
    vector<double >	params;
    void		print (FILE * fp) const;
  }
      ;

  struct gmx_exclusions_item
  {
    int			ii;
    vector<int >	params;
    void		print (FILE * fp) const;
  }
      ;

  struct gmx_bonds_item
  {
    int			ii;
    int			jj;
    int			funct;
    vector<double >	params;
    void		print (FILE * fp) const;
    int			atom_idx_ii () const {return ii-1;}
    int			atom_idx_jj () const {return jj-1;}
  }
      ;

  struct gmx_angles_item
  {
    int			ii;
    int			jj;
    int			kk;
    int			funct;
    vector<double >	params;
    void		print (FILE * fp) const;
    int			atom_idx_ii () const {return ii-1;}
    int			atom_idx_jj () const {return jj-1;}
    int			atom_idx_kk () const {return kk-1;}
  }
      ;

  struct gmx_dihedrals_item
  {
    int			ii;
    int			jj;
    int			kk;
    int			ll;
    int			funct;
    vector<double >	params;
    void		print (FILE * fp) const;
    int			atom_idx_ii () const {return ii-1;}
    int			atom_idx_jj () const {return jj-1;}
    int			atom_idx_kk () const {return kk-1;}
    int			atom_idx_ll () const {return ll-1;}
  }
      ;
  
  struct gmx_cmap_item
  {
    int			ii;
    int			jj;
    int			kk;
    int			ll;
    int			mm;
    int			funct;
    int			ngrid0;
    int			ngrid1;
    vector<double >	params;
    gmx_cmap_item	();
    void		print (FILE * fp) const;
    int			atom_idx_ii () const {return ii-1;}
    int			atom_idx_jj () const {return jj-1;}
    int			atom_idx_kk () const {return kk-1;}
    int			atom_idx_ll () const {return ll-1;}
    int			atom_idx_mm () const {return mm-1;}
  }
      ;
  
  struct gmx_mol
  {
    string			name;
    int				nexcl;
    vector<gmx_atom>		atoms;
    vector<gmx_pairs_item>	pairs;
    vector<gmx_exclusions_item>	exclusions;
    vector<gmx_bonds_item>	bonds;
    vector<gmx_angles_item>	angles;
    vector<gmx_dihedrals_item>	dihedrals;
    vector<gmx_cmap_item>	cmap;
    void		print (FILE * fp) const;
  }
      ;

  struct gmx_sys_top
  {
    string		sys_name;
    vector<gmx_mol>	moles;
    vector<int>		numMol;
    void		print (FILE * fp) const;
  }
      ;

  void parseTop (const string & fname,
		 gmx_sys_top & top);

  bool matchAtomType (const string & type,
		      const gmx_sys_types & systypes,
		      gmx_atomtypes_item & atomtype);
  bool matchBond (const int & iitype,
		  const int & jjtype,
		  const gmx_mol & systypes,
		  gmx_bonds_item & bond);
  bool matchBondType (const string & iitype,
		      const string & jjtype,
		      const gmx_sys_types & systypes,
		      gmx_bondtypes_item & bond_type);
  bool matchAngleType (const string & iitype,
		       const string & jjtype,
		       const string & kktype,
		       const gmx_sys_types & systypes,
		       gmx_angletypes_item & angletype);
  bool matchCMAPType (const string & iitype,
		      const string & jjtype,
		      const string & kktype,
		      const string & lltype,
		      const string & mmtype,
		      const gmx_sys_types & systypes,
		      gmx_cmaptypes_item & cmaptype,
		      int & idx);
};

#endif
