#ifndef __GmxTop_h_wanghan__
#define __GmxTop_h_wanghan__

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>


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
    gmx_atom ();
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

  struct gmx_bonds_item
  {
    int			ii;
    int			jj;
    int			funct;
    vector<double >	params;
    void		print (FILE * fp) const;
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
  }
      ;

  
  struct gmx_mol
  {
    string			name;
    int				nexcl;
    vector<gmx_atom>		atoms;
    vector<gmx_pairs_item>	pairs;
    vector<gmx_bonds_item>	bonds;
    vector<gmx_angles_item>	angles;
    vector<gmx_dihedrals_item>	dihedrals;
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
};

#endif
