/*
  File:      cloption.c

  Author: Xu KaiWen <xukaiwen@lsec.cc.ac.cn>
  Last modified: "2015-05-10 21:58:30 kaiwen"
*/

#include "cloption.h"
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

void GetCommandLineOptions(const int argc, char * argv[],
                           CommandLineOptions * cl){
   extern int opterr , optind;
   extern char* optarg;
   int opt;

   /* Note: turn off getopt() error checking */
   opterr = false;

   if (argc == 1){
      printf("Error, No command line option specified.\n");
      CommandLineOptionsHelp();
      exit(1);
   }

   while( ( opt = getopt( argc, argv, "h?cp" ) ) != -1 ){
      switch( opt ){
         case '?':         // getopt() does not recognize the argument
            printf("CommandLineOptions: ERROR: invalid command line option\n");
            CommandLineOptionsHelp();
            exit( 1 );
         case 'h':
            CommandLineOptionsHelp();
            exit( 0 );
         case 'c':{
            cl->if_convert = false;
            break;
         }
         case 'p':{
            cl->if_check = false;
            break;
         }
         default:
            printf("CommandLineOptions: ERROR: invalid command line option\n");
            CommandLineOptionsHelp();
            exit( 1 );
            break;
      }
   }
   if (argc > optind){
      strcpy(cl->filename, argv[optind]);
   } else {
      printf("CommandLineOptions: ERROR: no filename received.\n");
      exit(1);
   }
   return;
}

void CommandLineOptionsHelp(){
   printf("\n");
   printf("convert - convert 2d gmsh MSH2 file to easymesh-like .e .n files.\n");
   printf("\n");
   printf("Usage:\n");
   printf("   ./convert [OPTION] FILENAME\n");
   printf("Options:\n");
   printf("   -c, only check mesh validity. (Defualt is do both)\n");
   printf("   -p, only convert\n");
   printf("   -h, show help information\n");
   printf("Note:\n");
   printf("   1. GMSH mesh file is requared with .msh postfix, FILENAME.msh\n");
   printf("   2. Name of the output file is FILENAME.n for vertex\n");
   printf("      infomation and FILENAME.e for element related.\n");
   printf("   3. Format of .n file:\n");
   printf("      1) Fist line is the number of vertex.\n");
   printf("      2) Start from the second line is information of each vertex, in form:\n");
   printf("         [Point Index] [x coord.] [y coord.] [Point Mark]\n");
   printf("      3) Point Mark is the indicator of boundary type,\n");
   printf("         Point located on more than one edge is marked INTPTMK (-100)\n");
   printf("   4. Format of .e file:\n");
   printf("      1) Fist line is the number of element.\n");
   printf("      2) Start from the second line is information of each element, in form:\n");
   printf("         [Elemnt Index] [Vertex Index]X3 [Neighbour Index]X3 [Element Mark]\n");
   printf("      3) If one of the edges is overlaped with the boundary,\n");
   printf("         then there is no neighbour in that direction. \n");
   printf("         Instead, We use an negative value which is the\n");
   printf("         oppsite value of the boundary tag.\n");
}
