#include "TG3Application.h"
#include "TGeant3TGeo.h"
#include "TRint.h"
extern "C" 
{
   void uginit_();
   void uglast_();   
}
//______________________________________________________________________________
int main(int argc, char **argv)
{
   // Create an interactive ROOT application
   //TG3Application g3("G3","dummy");
   TG3Application g3;
   if (argc > 1 && !strcmp(argv[1],"TGeant3")) {
      new TGeant3("g",0);
      printf("+-----------------------------------------------------------+\n");
      printf("|                                                           |\n");
      printf("|       Running %s with TGeant3 <=========\n",argv[0]);
   } else {
      new TGeant3TGeo("g",0);
      printf("+-----------------------------------------------------------+\n");
      printf("|                                                           |\n");
      printf("|       Running %s with TGeant3TGeo <=========\n",argv[0]);
   }
   printf("|                                                           |\n");
   printf("+-----------------------------------------------------------+\n");
   
   uginit_();
   
   TRint theApp("Rint", &argc, argv);

   // and enter the event loop...
   theApp.Run(kTRUE);

   uglast_();

   return 0;
}
