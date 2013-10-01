#include "GmxTop.h"
#include "GmxType.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;

int main (int argc, char **argv) {

  std::string ifile, ofile, opfile;
  double sbond, sangle, sdihedral;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")      
      ("scal-bond,b", po::value<double > (&sbond)->default_value (1.0), "bond scale")
      ("scal-angle,a", po::value<double > (&sangle)->default_value (1.0), "angle scale")
      ("scal-dihedral,d", po::value<double > (&sdihedral)->default_value (1.0), "dihedral scale")
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

  if (systype.bondtypes.size() > 0){
    for (unsigned ii = 0; ii < systype.bondtypes.size(); ++ii){
      if (systype.bondtypes[ii].funct == 1 ||
	  systype.bondtypes[ii].funct == 2 ||
	  systype.bondtypes[ii].funct == 6 ||
	  systype.bondtypes[ii].funct == 7 ||
	  systype.bondtypes[ii].funct == 8 ||
	  systype.bondtypes[ii].funct == 9 ){
	systype.bondtypes[ii].params[1] *= sbond;
      }
      else {
	cerr << "bond type is not 1,2,6,7,8,9, do nothing!" << endl;
      }
    }
  }

  if (systype.angletypes.size() > 0){
    for (unsigned ii = 0; ii < systype.angletypes.size(); ++ii){
      if (systype.angletypes[ii].funct == 1 ||
	  systype.angletypes[ii].funct == 2 ||
	  systype.angletypes[ii].funct == 8 ){
	systype.angletypes[ii].params[1] *= sangle;
      }
      else {
	cerr << "angle type is not 1,2,8, do nothing!" << endl;
      }
    }
  }
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
	cerr << "angle type is not 1,2,4,8,9, do nothing!" << endl;
      }
    }
  }

  
  FILE * fp = fopen (ofile.c_str(), "w");
  systype.print (fp);
  systop.print (fp);

  return 0;
}
