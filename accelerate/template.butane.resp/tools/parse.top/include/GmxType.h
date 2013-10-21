#ifndef __GmxType_h_wanghan__
#define __GmxType_h_wanghan__


#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;


namespace GmxTop {
  struct gmx_defaults_item
  {
    int		nbfunc;
    int		comb_rule;
    bool	gen_pair;
    double	fudgeLJ;
    double	fudgeQQ;
    void	print (FILE * fp) const;
  }
      ;

  struct gmx_atomtypes_item
  {
    string	name;
    int		atom_num;
    double	mass;
    double	charge;
    string	ptype;
    double	c6;
    double	c12;
    void	print (FILE * fp) const;
  }
      ;

  struct gmx_nonbond_params_item
  {
    string		name0;
    string		name1;
    int			funct;
    vector<double >	params;
    void		print (FILE * fp) const;
  }
      ;
  
  struct gmx_pairtypes_item
  {
    string		name0;
    string		name1;
    int			funct;
    vector<double >	params;
    void		print (FILE * fp) const;
  }
      ;
  
  struct gmx_bondtypes_item
  {
    string		name0;
    string		name1;
    int			funct;
    vector<double >	params;
    void		print (FILE * fp) const;
  }
      ;

  struct gmx_angletypes_item
  {
    string		name0;
    string		name1;
    string		name2;
    int			funct;
    vector<double >	params;
    void		print (FILE * fp) const;
  }
      ;

  struct gmx_dihedraltypes_item
  {
    string		name0;
    string		name1;
    string		name2;
    string		name3;
    int			funct;
    vector<double >	params;
    void		print (FILE * fp) const;
  }
      ;

  struct gmx_cmaptypes_item
  {
    string		name0;
    string		name1;
    string		name2;
    string		name3;
    string		name4;
    int			funct;
    int			ngrid0;
    int			ngrid1;    
    vector<double >	params;
    void		print (FILE * fp) const;
  }
      ;

  struct gmx_sys_types
  {
    gmx_defaults_item			defaults;
    vector<gmx_atomtypes_item>		atomtypes;
    vector<gmx_pairtypes_item>		pairtypes;
    vector<gmx_nonbond_params_item>	nonbond_params;
    vector<gmx_bondtypes_item>		bondtypes;
    vector<gmx_angletypes_item>		angletypes;
    vector<gmx_dihedraltypes_item>	dihedraltypes;
    vector<gmx_cmaptypes_item>		cmaptypes;
    void		print (FILE * fp) const;
  }
      ;
  
  void parseType (const string & fname,
		  gmx_sys_types & type);

  void convertType_2_1 (const gmx_sys_types	& old_types,
			gmx_sys_types		& new_types);
}


#endif
