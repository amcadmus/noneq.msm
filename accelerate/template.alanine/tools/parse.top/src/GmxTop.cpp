#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
#include <cmath>
#include <string.h>

#include "GmxTop.h"
#include "GmxType.h"
#include "StringSplit.h"

using namespace GmxTop;

static inline void
die_wrong_format (const string & file,
		  const int & line)
{
  cerr << "wrong format error happing at file "
       << file << " : "
       << line << endl;
  exit (1);
}

static void
normalizeLine (char *line)
{
  int ii = 0;
  while (line[ii] != '\0'){
    if (line[ii] == ';'){
      line[ii] = '\0';
      return;
    }
    ii++;
  }
}

static bool
ifKeyWord (const string & in_,
	   string & key)
{
  char tmpin [MAX_LINE_LENGTH];
  strncpy (tmpin, in_.c_str(), MAX_LINE_LENGTH);
  normalizeLine (tmpin);
  string in (tmpin);
  
  if (!(in[0] == '[')){
    return false;
  }
  key.clear ();
  
  for (unsigned ii = 0; ii < in.size(); ++ii){
    if (in[ii] != '[' && in[ii] != ']' && in[ii] != ' '){
      key.push_back (in[ii]);
    }
  }
  return true;
}


static bool
notTrivalLine (const string & line)
{
  if (line.size() > 0 && line[0] != ';' && line[0] != '*'){
    vector<string > words;
    StringOperation::split (line, words);
    if (words.size() > 0){
      return true;
    }
  }
  return false;
}

static void 
readBlocks (ifstream & file,
	    vector<string> & keys,
	    vector<vector<string > > & lines)
{
  bool inBlock = false;
  vector<string > blockLines;
  string tmpKey;
  char line[MAX_LINE_LENGTH];  
  
  while (! file.eof() ){
    file.getline (line, MAX_LINE_LENGTH);
    
    // cout << line <<endl;
    
    if (!inBlock){
      if (ifKeyWord(line, tmpKey)){
	normalizeLine (line);
	inBlock = true;
	// keys.push_back ("");
	keys.push_back (tmpKey);
	// lines.push_back (blockLines);
	blockLines.clear ();
      }
      else if (notTrivalLine(line)) {
	normalizeLine (line);
	blockLines.push_back(line);	
      }
    }
    else {
      if (ifKeyWord(line, tmpKey)){
	normalizeLine (line);
	keys.push_back (tmpKey);
	lines.push_back (blockLines);
	blockLines.clear ();	
      }
      else if (notTrivalLine(line)){
	normalizeLine (line);
	blockLines.push_back(line);
      }
    }
  }
  lines.push_back (blockLines);
}


// void
// realMol (ifstream & file,
// 	 gmx_mol & mol,
// 	 string & key)
// {
//   while (!file.eof()){
//     file.getline (line, MAX_LINE_LENGTH);
//     vector<string > words;
//     if (line.size() > 0 && line[0] != ';'){
//       StringOperation::split (line, words);
//       if (words.size() == 2){
// 	mol.name = words[0];
//       }
//     }    
//   }
// }

GmxTop::gmx_atom::
gmx_atom ()
    :  id(0), res_nr(0), cgnr(1), charge (0.0), mass (0.0)
{
}

void GmxTop::gmx_atom::
clear ()
{
  id = 0;
  at_type = "";
  res_nr = 0;
  res_name = "";
  at_name = "";
  cgnr = 1;
  charge = 0.;
  mass = 0.;
}


void GmxTop::gmx_atom::
print (FILE * fp) const
{
  fprintf (fp, "%d\t%s\t%d\t%s\t%s\t%d\t%f\t%f\n",
	   id, at_type.c_str(), res_nr, res_name.c_str(), at_name.c_str(), cgnr, charge, mass);
}

void GmxTop::gmx_pairs_item::
print (FILE * fp) const
{
  fprintf (fp, "%d\t%d\t%d", ii, jj, funct);
  for (unsigned ii = 0; ii < params.size(); ++ii){
    fprintf (fp, "\t%f", params[ii]);
  }
  fprintf (fp, "\n");
}

void GmxTop::gmx_bonds_item::
print (FILE * fp) const
{
  fprintf (fp, "%d\t%d\t%d", ii, jj, funct);
  for (unsigned ii = 0; ii < params.size(); ++ii){
    fprintf (fp, "\t%f", params[ii]);
  }
  fprintf (fp, "\n");
}

void GmxTop::gmx_angles_item::
print (FILE * fp) const
{
  fprintf (fp, "%d\t%d\t%d\t%d", ii, jj, kk, funct);
  for (unsigned ii = 0; ii < params.size(); ++ii){
    fprintf (fp, "\t%f", params[ii]);
  }
  fprintf (fp, "\n");
}

void GmxTop::gmx_dihedrals_item::
print (FILE * fp) const
{
  fprintf (fp, "%d\t%d\t%d\t%d\t%d", ii, jj, kk, ll, funct);
  for (unsigned ii = 0; ii < params.size(); ++ii){
    fprintf (fp, "\t%f", params[ii]);
  }
  fprintf (fp, "\n");
}


void GmxTop::gmx_mol::
print (FILE * fp) const
{
  fprintf (fp, "[ moleculetype ]\n");
  fprintf (fp, "%s\t%d\n", name.c_str(), nexcl);
  fprintf (fp, "\n");
  fprintf (fp, "[ atoms ]\n");
  for (unsigned ii = 0; ii < atoms.size(); ++ii){
    atoms[ii].print (fp);
  }
  fprintf (fp, "\n");
  if (pairs.size() > 0) {
    fprintf (fp, "[ pairs ]\n");
    for (unsigned ii = 0; ii < pairs.size(); ++ii){
      pairs[ii].print (fp);
    }
    fprintf (fp, "\n");
  }
  if (bonds.size() > 0) {
    fprintf (fp, "[ bonds ]\n");
    for (unsigned ii = 0; ii < bonds.size(); ++ii){
      bonds[ii].print (fp);
    }
    fprintf (fp, "\n");
  }
  if (angles.size() > 0) {
    fprintf (fp, "[ angles ]\n");
    for (unsigned ii = 0; ii < angles.size(); ++ii){
      angles[ii].print (fp);
    }
    fprintf (fp, "\n");
  }
  if (dihedrals.size() > 0) {
    fprintf (fp, "[ dihedrals ]\n");
    for (unsigned ii = 0; ii < dihedrals.size(); ++ii){
      dihedrals[ii].print (fp);
    }
    fprintf (fp, "\n");
  }  
}

void GmxTop::gmx_sys_top::
print (FILE * fp) const
{
  for (unsigned ii = 0; ii < moles.size(); ++ii){
    moles[ii].print (fp);
  }
  fprintf (fp, "[ system ]\n");
  fprintf (fp, "%s\n", sys_name.c_str());
  fprintf (fp, "\n");
  fprintf (fp, "[ molecules ]\n");
  for (unsigned ii = 0; ii < moles.size(); ++ii){
    fprintf (fp, "%s\t%d\n",
	     moles[ii].name.c_str(), numMol[ii]);
  }
}

void GmxTop::gmx_sys_types::
print (FILE * fp) const
{
  fprintf (fp, "[ defaults ]\n");
  defaults.print (fp);
  fprintf (fp, "\n");
  if (atomtypes.size() > 0){
    fprintf (fp, "[ atomtypes ]\n");
    for (unsigned ii = 0; ii < atomtypes.size(); ++ii){
      atomtypes[ii].print (fp);
    }
    fprintf (fp, "\n");
  }
  if (nonbond_params.size() > 0){
    fprintf (fp, "[ nonbond_params ]\n");
    for (unsigned ii = 0; ii < nonbond_params.size(); ++ii){
      nonbond_params[ii].print (fp);
    }
    fprintf (fp, "\n");
  }
  if (pairtypes.size() > 0){
    fprintf (fp, "[ pairtypes ]\n");
    for (unsigned ii = 0; ii < pairtypes.size(); ++ii){
      pairtypes[ii].print (fp);
    }
    fprintf (fp, "\n");
  }
  if (bondtypes.size() > 0){
    fprintf (fp, "[ bondtypes ]\n");
    for (unsigned ii = 0; ii < bondtypes.size(); ++ii){
      bondtypes[ii].print (fp);
    }
    fprintf (fp, "\n");
  }
  if (angletypes.size() > 0){  
    fprintf (fp, "[ angletypes ]\n");
    for (unsigned ii = 0; ii < angletypes.size(); ++ii){
      angletypes[ii].print (fp);
    }
    fprintf (fp, "\n");
  }
  if (dihedraltypes.size() > 0){
    fprintf (fp, "[ dihedraltypes ]\n");
    for (unsigned ii = 0; ii < dihedraltypes.size(); ++ii){
      dihedraltypes[ii].print (fp);
    }
    fprintf (fp, "\n");
  }
}


// void
void GmxTop::
parseTop (const string & fname,
	  gmx_sys_top & top)
{
  ifstream file (fname.c_str());
  if (! file.is_open()){
    cerr << "cannot open file " << fname << endl;
    exit (1);
  }

  vector<string > keys;
  vector<vector<string > > lines;
  readBlocks (file, keys, lines);
  vector<string > words;

  // cout << "n. keys " << keys.size() << endl;
  // cout << "n. lines " << lines.size() << endl;
  // for (unsigned ii = 0; ii < keys.size(); ++ii){
  //   cout << "key " << ii << ": " << keys[ii] << endl;
  //   for (unsigned jj = 0; jj < lines[ii].size(); ++jj){
  //     cout << lines[ii][jj] << endl;
  //   }
  // }
  
  for (unsigned ii = 0; ii < keys.size(); ++ii){
    if (keys[ii] == "system"){
      if (lines[ii].size() < 1) die_wrong_format (__FILE__, __LINE__);
      top.sys_name = lines[ii][0];
    }
  }
      
  for (unsigned ii = 0; ii < keys.size(); ++ii){
    if (keys[ii] == "moleculetype"){
      gmx_mol tmpmol;
      StringOperation::split (lines[ii][0], words);
      if (words.size() < 2) die_wrong_format (__FILE__, __LINE__);
      tmpmol.name = words[0];
      tmpmol.nexcl = atoi(words[1].c_str());
      unsigned jj = ii+1; 
      for (;jj < keys.size(); ++jj){
	if (keys[jj] == "moleculetype"){
	  break;
	}
      }
      for (unsigned kk = ii+1; kk < jj; ++kk){
	if (keys[kk] == "atoms"){
	  gmx_atom tmpatom;
	  for (unsigned ll = 0; ll < lines[kk].size(); ++ll){
	    tmpatom.clear ();
	    StringOperation::split (lines[kk][ll], words);
	    if (words.size() >= 1){
	      tmpatom.id = atoi(words[0].c_str());
	    }
	    if (words.size() >= 2){
	      tmpatom.at_type = words[1];
	    }
	    if (words.size() >= 3){
	      tmpatom.res_nr = atoi(words[2].c_str());
	    }
	    if (words.size() >= 4){
	      tmpatom.res_name = words[3];
	    }
	    if (words.size() >= 5){
	      tmpatom.at_name = words[4];
	    }
	    if (words.size() >= 6){
	      tmpatom.cgnr = atoi (words[5].c_str());
	    }
	    if (words.size() >= 7){
	      tmpatom.charge = atof (words[6].c_str());
	    }
	    if (words.size() >= 8){
	      tmpatom.mass = atof (words[7].c_str());
	    }
	    tmpmol.atoms.push_back (tmpatom);
	  }
	}
	if (keys[kk] == "pairs"){
	  for (unsigned ll = 0; ll < lines[kk].size(); ++ll){
	    StringOperation::split (lines[kk][ll], words);
	    gmx_pairs_item tmp;
	    if (words.size () < 3) die_wrong_format (__FILE__, __LINE__);
	    tmp.ii = atoi(words[0].c_str());
	    tmp.jj = atoi(words[1].c_str());
	    tmp.funct = atoi(words[2].c_str());
	    for (unsigned mm = 3; mm < words.size(); ++mm){
	      tmp.params.push_back (atof(words[mm].c_str()));
	    }
	    tmpmol.pairs.push_back(tmp);
	  } 
	}
	if (keys[kk] == "bonds"){
	  for (unsigned ll = 0; ll < lines[kk].size(); ++ll){
	    StringOperation::split (lines[kk][ll], words);
	    gmx_bonds_item tmp;
	    if (words.size () < 3) die_wrong_format (__FILE__, __LINE__);
	    tmp.ii = atoi(words[0].c_str());
	    tmp.jj = atoi(words[1].c_str());
	    tmp.funct = atoi(words[2].c_str());
	    for (unsigned mm = 3; mm < words.size(); ++mm){
	      tmp.params.push_back (atof(words[mm].c_str()));
	    }
	    tmpmol.bonds.push_back(tmp);
	  } 
	}
	if (keys[kk] == "angles"){
	  for (unsigned ll = 0; ll < lines[kk].size(); ++ll){
	    StringOperation::split (lines[kk][ll], words);
	    gmx_angles_item tmp;
	    if (words.size () < 4) die_wrong_format (__FILE__, __LINE__);
	    tmp.ii = atoi(words[0].c_str());
	    tmp.jj = atoi(words[1].c_str());
	    tmp.kk = atoi(words[2].c_str());
	    tmp.funct = atoi(words[3].c_str());
	    for (unsigned mm = 4; mm < words.size(); ++mm){
	      tmp.params.push_back (atof(words[mm].c_str()));
	    }
	    tmpmol.angles.push_back(tmp);
	  } 
	}
	if (keys[kk] == "dihedrals"){
	  for (unsigned ll = 0; ll < lines[kk].size(); ++ll){
	    StringOperation::split (lines[kk][ll], words);
	    gmx_dihedrals_item tmp;
	    if (words.size () < 5) die_wrong_format (__FILE__, __LINE__);
	    tmp.ii = atoi(words[0].c_str());
	    tmp.jj = atoi(words[1].c_str());
	    tmp.kk = atoi(words[2].c_str());
	    tmp.ll = atoi(words[3].c_str());
	    tmp.funct = atoi(words[4].c_str());
	    for (unsigned mm = 5; mm < words.size(); ++mm){
	      tmp.params.push_back (atof(words[mm].c_str()));
	    }
	    tmpmol.dihedrals.push_back(tmp);
	  } 
	}
      }
      top.moles.push_back (tmpmol);
      top.numMol.push_back (0);
      ii = jj - 1;
    }
  }

  for (unsigned ii = 0; ii < keys.size(); ++ii){
    if (keys[ii] == "molecules"){
      // if (lines[ii].size() != top.moles.size()){
      // 	cerr << "num. mole does not match the num. of lines in sec. [ molecules ], stop." << endl;
      // 	exit(1);
      // }
      for (unsigned jj = 0; jj < lines[ii].size(); ++jj){
	StringOperation::split (lines[ii][jj], words);
	if (words.size () < 2) die_wrong_format (__FILE__, __LINE__);
	for (unsigned kk = 0; kk < top.moles.size(); ++kk){
	  if (top.moles[kk].name == string(words[0])){
	    top.numMol[kk] = atoi(words[1].c_str());
	    break;
	  }
	}
      }
    }
  }   
}


void GmxTop::
parseType (const string & fname,
	   gmx_sys_types & type)
{
  ifstream file (fname.c_str());
  if (! file.is_open()){
    cerr << "cannot open file " << fname << endl;
    exit (1);
  }
  type.atomtypes.clear();
  type.bondtypes.clear();
  type.angletypes.clear();
  type.dihedraltypes.clear();

  vector<string > keys;
  vector<vector<string > > lines;
  readBlocks (file, keys, lines);
  vector<string > words;
  
  for (unsigned ii = 0; ii < keys.size(); ++ii){
    if (keys[ii] == "defaults"){
      StringOperation::split (lines[ii][0], words);
      if (words.size() < 5) die_wrong_format (__FILE__, __LINE__);
      type.defaults.nbfunc = atoi(words[0].c_str());
      type.defaults.comb_rule = atoi(words[1].c_str());
      if (words[2] == "yes"){
	type.defaults.gen_pair = true;
      }
      else {
	type.defaults.gen_pair = false;
      }
      type.defaults.fudgeLJ = atof(words[3].c_str());
      type.defaults.fudgeQQ = atof(words[4].c_str());
      break;
    }
  }

  for (unsigned ii = 0; ii < keys.size(); ++ii){
    if (keys[ii] == "atomtypes"){
      for (unsigned jj = 0; jj < lines[ii].size(); ++jj){
	StringOperation::split (lines[ii][jj], words);
	if (words.size() < 7) die_wrong_format (__FILE__, __LINE__);
	gmx_atomtypes_item tmp;
	tmp.name = words[0];
	tmp.atom_num = atoi(words[1].c_str());
	tmp.mass = atoi(words[2].c_str());
	tmp.charge = atoi(words[3].c_str());
	tmp.ptype = string(words[4].c_str());
	tmp.c6 = atof(words[5].c_str());
	tmp.c12 = atof(words[6].c_str());
	type.atomtypes.push_back (tmp);
      }
    }
  }

  for (unsigned ii = 0; ii < keys.size(); ++ii){
    if (keys[ii] == "bondtypes"){
      for (unsigned jj = 0; jj < lines[ii].size(); ++jj){
	StringOperation::split (lines[ii][jj], words);
	if (words.size() < 3) die_wrong_format (__FILE__, __LINE__);
	gmx_bondtypes_item tmp;
	tmp.name0 = words[0];
	tmp.name1 = words[1];
	tmp.funct = atoi(words[2].c_str());
	for (unsigned ii = 3; ii < words.size(); ++ii){
	  tmp.params.push_back (atof(words[ii].c_str()));
	}
	type.bondtypes.push_back (tmp);
      }
    }
  }
  
  for (unsigned ii = 0; ii < keys.size(); ++ii){
    if (keys[ii] == "pairtypes"){
      for (unsigned jj = 0; jj < lines[ii].size(); ++jj){
	StringOperation::split (lines[ii][jj], words);
	if (words.size() < 3) die_wrong_format (__FILE__, __LINE__);
	gmx_pairtypes_item tmp;
	tmp.name0 = words[0];
	tmp.name1 = words[1];
	tmp.funct = atoi(words[2].c_str());
	for (unsigned ii = 3; ii < words.size(); ++ii){
	  tmp.params.push_back (atof(words[ii].c_str()));
	}
	type.pairtypes.push_back (tmp);
      }
    }
  }

  for (unsigned ii = 0; ii < keys.size(); ++ii){
    if (keys[ii] == "nonbond_params"){
      for (unsigned jj = 0; jj < lines[ii].size(); ++jj){
	StringOperation::split (lines[ii][jj], words);
	if (words.size() < 3) die_wrong_format (__FILE__, __LINE__);
	gmx_nonbond_params_item tmp;
	tmp.name0 = words[0];
	tmp.name1 = words[1];
	tmp.funct = atoi(words[2].c_str());
	for (unsigned ii = 3; ii < words.size(); ++ii){
	  tmp.params.push_back (atof(words[ii].c_str()));
	}
	type.nonbond_params.push_back (tmp);
      }
    }
  }
  
  for (unsigned ii = 0; ii < keys.size(); ++ii){
    if (keys[ii] == "angletypes"){
      for (unsigned jj = 0; jj < lines[ii].size(); ++jj){
	StringOperation::split (lines[ii][jj], words);
	if (words.size() < 4) die_wrong_format (__FILE__, __LINE__);
	gmx_angletypes_item tmp;
	tmp.name0 = words[0];
	tmp.name1 = words[1];
	tmp.name2 = words[2];
	tmp.funct = atoi(words[3].c_str());
	for (unsigned ii = 4; ii < words.size(); ++ii){
	  tmp.params.push_back (atof(words[ii].c_str()));
	}
	type.angletypes.push_back (tmp);
      }
    }
  }

  for (unsigned ii = 0; ii < keys.size(); ++ii){
    if (keys[ii] == "dihedraltypes"){
      for (unsigned jj = 0; jj < lines[ii].size(); ++jj){
	StringOperation::split (lines[ii][jj], words);
	if (words.size() < 5) die_wrong_format (__FILE__, __LINE__);
	gmx_dihedraltypes_item tmp;
	tmp.name0 = words[0];
	tmp.name1 = words[1];
	tmp.name2 = words[2];
	tmp.name3 = words[3];
	tmp.funct = atoi(words[4].c_str());
	for (unsigned ii = 5; ii < words.size(); ++ii){
	  tmp.params.push_back (atof(words[ii].c_str()));
	}
	type.dihedraltypes.push_back (tmp);
      }
    }
  }
}

static void
param_2_to_1 (const double & eps,
	      const double & sigma,
	      double & c6,
	      double & c12)
{
  c6 = sigma * sigma * sigma;
  c6 = c6 * c6;
  c12 = 4. * eps * c6 * c6;
  c6  = 4. * eps * c6;
}

void GmxTop::
convertType_2_1 (const gmx_sys_types	& old_types,
		 gmx_sys_types		& new_types)
{
  if (old_types.defaults.comb_rule != 2){
    cout << "old comb rule is not 2, return without doing anything" <<endl;
    return ;
  }
  
  new_types = old_types;

  new_types.defaults.comb_rule = 1;
  new_types.nonbond_params.clear();
  
  for (unsigned ii = 0; ii < new_types.atomtypes.size(); ++ii){
    double eps0 = old_types.atomtypes[ii].c12;
    double sigma0 = old_types.atomtypes[ii].c6;
    double c6_0, c12_0;
    param_2_to_1 (eps0, sigma0, c6_0, c12_0);
    new_types.atomtypes[ii].c6 = c6_0;
    new_types.atomtypes[ii].c12 = c12_0;
    for (unsigned jj = ii; jj < new_types.atomtypes.size(); ++jj){
      gmx_nonbond_params_item tmp;
      tmp.name0 = old_types.atomtypes[ii].name;
      tmp.name1 = old_types.atomtypes[jj].name;
      double eps1 = old_types.atomtypes[jj].c12;
      double sigma1 = old_types.atomtypes[jj].c6;
      double eps01 = sqrt(eps0 * eps1);
      double sigma01 = 0.5 * (sigma0 + sigma1);
      double c6_01, c12_01;
      param_2_to_1 (eps01, sigma01, c6_01, c12_01);
      tmp.funct = 1;
      tmp.params.resize (2);
      tmp.params[0] = c6_01;
      tmp.params[1] = c12_01;
      new_types.nonbond_params.push_back(tmp);
    }
  }    
}



