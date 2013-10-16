#include "GmxType.h"

void GmxTop::gmx_defaults_item::
print (FILE * fp) const 
{
  if (gen_pair){
    fprintf (fp, "%d\t%d\t yes\t%f\t%f\n",
	     nbfunc, comb_rule, fudgeLJ, fudgeQQ);
  }
  else {
    fprintf (fp, "%d\t%d\t no \t%f\t%f\n",
	     nbfunc, comb_rule, fudgeLJ, fudgeQQ);
  }
}

void GmxTop::gmx_atomtypes_item::
print (FILE * fp) const
{
  fprintf (fp, "%s\t%d\t%f\t%f\t%s\t%.10e\t%.10e\n",
	   name.c_str(), atom_num, mass, charge, ptype.c_str(), c6, c12);
}

void GmxTop::gmx_nonbond_params_item::
print (FILE * fp) const
{
  fprintf (fp, "%s\t%s\t%d",
	   name0.c_str(), name1.c_str(), funct);
  for (unsigned ii = 0; ii < params.size(); ++ii){
    fprintf (fp, "\t%.10e", params[ii]);
  }
  fprintf (fp, "\n");
}

void GmxTop::gmx_pairtypes_item::
print (FILE * fp) const
{
  fprintf (fp, "%s\t%s\t%d",
	   name0.c_str(), name1.c_str(), funct);
  for (unsigned ii = 0; ii < params.size(); ++ii){
    fprintf (fp, "\t%.10e", params[ii]);
  }
  fprintf (fp, "\n");
}

void GmxTop::gmx_bondtypes_item::
print (FILE * fp) const
{
  fprintf (fp, "%s\t%s\t%d",
	   name0.c_str(), name1.c_str(), funct);
  for (unsigned ii = 0; ii < params.size(); ++ii){
    fprintf (fp, "\t%.10e", params[ii]);
  }
  fprintf (fp, "\n");
}
  
void GmxTop::gmx_angletypes_item::
print (FILE * fp) const
{
  fprintf (fp, "%s\t%s\t%s\t%d",
	   name0.c_str(), name1.c_str(), name2.c_str(), funct);
  for (unsigned ii = 0; ii < params.size(); ++ii){
    fprintf (fp, "\t%.10e", params[ii]);
  }
  fprintf (fp, "\n");
}

void GmxTop::gmx_dihedraltypes_item::
print (FILE * fp) const
{
  fprintf (fp, "%s\t%s\t%s\t%s\t%d",
	   name0.c_str(), name1.c_str(), name2.c_str(), name3.c_str(), funct);
  if (funct == 1 || funct == 4 || funct == 9){
    if (params.size() != 0){
      for (int myii = 0; myii < int(params.size()) - 1; ++myii){
	fprintf (fp, "\t%.10e", params[myii]);
      }
      fprintf (fp, "\t%d", int(params.back() + 0.5));
    }
  }
  else{
    for (int myii = 0; myii < int(params.size()); ++myii){
      fprintf (fp, "\t%.10e", params[myii]);
    }
  }
  fprintf (fp, "\n");
}

void GmxTop::gmx_cmaptypes_item::  
print (FILE * fp) const
{
  fprintf (fp, "%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d",
	   name0.c_str(), name1.c_str(), name2.c_str(), name3.c_str(), name4.c_str(),
	   funct,
	   ngrid0, ngrid1);
  fprintf (fp, "\\\n");
  for (unsigned ii = 0; ii < params.size(); ++ii){
    fprintf (fp, "%.10f ", params[ii]);
    if ((ii+1) % 10 == 0){
      fprintf (fp, "\\\n");
    }
  }
  fprintf (fp, "\n");
}

