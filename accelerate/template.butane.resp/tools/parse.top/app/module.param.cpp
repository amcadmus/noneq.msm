#include "GmxTop.h"
#include "GmxType.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;

int main (int argc, char **argv) {

  std::string ifile, ofile, opfile;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("output,o", po::value<std::string > (&ofile)->default_value ("out.top"), "the output top")
      ("output-param", po::value<std::string > (&opfile)->default_value ("top.param.out"), "the output param")
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

  FILE * fp = fopen (ofile.c_str(), "w");
  FILE * fp_param = fopen (opfile.c_str(), "w");
  
  fprintf (fp, "[ defaults ]\n");
  systype.defaults.print (fp);
  fprintf (fp, "\n");
  if (systype.atomtypes.size() > 0){
    fprintf (fp, "[ atomtypes ]\n");
    for (unsigned ii = 0; ii < systype.atomtypes.size(); ++ii){
      systype.atomtypes[ii].print (fp);
    }
    fprintf (fp, "\n");
  }
  if (systype.nonbond_params.size() > 0){
    fprintf (fp, "[ nonbond_params ]\n");
    for (unsigned ii = 0; ii < systype.nonbond_params.size(); ++ii){
      systype.nonbond_params[ii].print (fp);
    }
    fprintf (fp, "\n");
  }
  if (systype.pairtypes.size() > 0){
    fprintf (fp, "[ pairtypes ]\n");
    for (unsigned ii = 0; ii < systype.pairtypes.size(); ++ii){
      systype.pairtypes[ii].print (fp);
    }
    fprintf (fp, "\n");
  }

  if (systype.bondtypes.size() > 0){
    fprintf (fp, "[ bondtypes ]\n");
    for (unsigned ii = 0; ii < systype.bondtypes.size(); ++ii){
      fprintf (fp, "%s\t%s\t%d",
	       systype.bondtypes[ii].name0.c_str(),
	       systype.bondtypes[ii].name1.c_str(),
	       systype.bondtypes[ii].funct);
      if (systype.bondtypes[ii].funct == 1 ||
	  systype.bondtypes[ii].funct == 2 ||
	  systype.bondtypes[ii].funct == 6 ||
	  systype.bondtypes[ii].funct == 7 ||
	  systype.bondtypes[ii].funct == 8 ||
	  systype.bondtypes[ii].funct == 9 ){
	fprintf (fp, "\t%f", systype.bondtypes[ii].params[0]);
	char param_name [MAX_LINE_LENGTH];
	sprintf (param_name, "bond_k_%s_%s", 
		 systype.bondtypes[ii].name0.c_str(),
		 systype.bondtypes[ii].name1.c_str());
	fprintf (fp, "\t%s", param_name);
	fprintf (fp_param, "%s=%f\n", param_name, systype.bondtypes[ii].params[1]);
      }
      else {
	for (unsigned kk = 0; kk < systype.bondtypes[ii].params.size(); ++kk){
	  fprintf (fp, "\t%f", systype.bondtypes[ii].params[kk]);
	}
      }
      fprintf (fp, "\n");
    }
    fprintf (fp, "\n");
  }

  if (systype.angletypes.size() > 0){
    fprintf (fp, "[ angletypes ]\n");
    for (unsigned ii = 0; ii < systype.angletypes.size(); ++ii){
      fprintf (fp, "%s\t%s\t%s\t%d",
	       systype.angletypes[ii].name0.c_str(),
	       systype.angletypes[ii].name1.c_str(),
	       systype.angletypes[ii].name2.c_str(),
	       systype.angletypes[ii].funct);
      if (systype.angletypes[ii].funct == 1 ||
	  systype.angletypes[ii].funct == 2 ){
	fprintf (fp, "\t%f", systype.angletypes[ii].params[0]);
	char param_name [MAX_LINE_LENGTH];
	sprintf (param_name, "angle_k_%s_%s_%s", 
		 systype.angletypes[ii].name0.c_str(),
		 systype.angletypes[ii].name1.c_str(),
		 systype.angletypes[ii].name2.c_str());
	fprintf (fp, "\t%s", param_name);
	fprintf (fp_param, "%s=%f\n", param_name, systype.angletypes[ii].params[1]);
      }
      else {
	for (unsigned kk = 0; kk < systype.angletypes[ii].params.size(); ++kk){
	  fprintf (fp, "\t%f", systype.angletypes[ii].params[kk]);
	}
      }
      fprintf (fp, "\n");
    }
    fprintf (fp, "\n");
  }

  if (systype.dihedraltypes.size() > 0){
    fprintf (fp, "[ dihedraltypes ]\n");
    for (unsigned ii = 0; ii < systype.dihedraltypes.size(); ++ii){
      fprintf (fp, "%s\t%s\t%s\t%s\t%d",
	       systype.dihedraltypes[ii].name0.c_str(),
	       systype.dihedraltypes[ii].name1.c_str(),
	       systype.dihedraltypes[ii].name2.c_str(),
	       systype.dihedraltypes[ii].name3.c_str(),
	       systype.dihedraltypes[ii].funct);
      if (systype.dihedraltypes[ii].funct == 1 ||
	  systype.dihedraltypes[ii].funct == 4 ||
	  systype.dihedraltypes[ii].funct == 9 ){
	fprintf (fp, "\t%f", systype.dihedraltypes[ii].params[0]);
	char param_name [MAX_LINE_LENGTH];
	sprintf (param_name, "dihedral_k_%s_%s_%s_%s", 
		 systype.dihedraltypes[ii].name0.c_str(),
		 systype.dihedraltypes[ii].name1.c_str(),
		 systype.dihedraltypes[ii].name2.c_str(),
		 systype.dihedraltypes[ii].name3.c_str());
	fprintf (fp, "\t%s", param_name);
	fprintf (fp_param, "%s=%f\n", param_name, systype.dihedraltypes[ii].params[1]);
	fprintf (fp, "\t%d", int (systype.dihedraltypes[ii].params[2] + 0.5));
      }
      else {
	for (unsigned kk = 0; kk < systype.dihedraltypes[ii].params.size(); ++kk){
	  fprintf (fp, "\t%f", systype.dihedraltypes[ii].params[kk]);
	}
      }
      fprintf (fp, "\n");
    }
    fprintf (fp, "\n");
  }
  
  systop.print (fp);

  return 0;
}
