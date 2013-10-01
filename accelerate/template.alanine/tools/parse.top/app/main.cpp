#include "GmxTop.h"
#include "GmxType.h"


int main (int argc, char **argv) {
  GmxTop::gmx_sys_top systop;
  GmxTop::parseTop (argv[1], systop);
  GmxTop::gmx_sys_types systype;
  GmxTop::parseType (argv[1], systype);

  // cout << endl;
  // for (unsigned kk = 0; kk < systop.moles.size(); ++kk){
  //   cout << systop.moles[kk].name
  // 	 << "\t"
  // 	 << systop.numMol[kk]
  // 	 << endl;
  // }
  // cout << endl;
  
  // systype.print(stdout);
  // systop.print (stdout);

  GmxTop::gmx_sys_types newtype;
  GmxTop::convertType_2_1 (systype, newtype);
      
  newtype.print(stdout);
  systop.print (stdout);

  return 0;
}
