/*
  File:      cloption.c

  Author: Xu KaiWen <xukaiwen@lsec.cc.ac.cn>
  Last modified: "2015-06-03 23:48:21 kaiwen"
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

   while( ( opt = getopt( argc, argv, "h?cpn:" ) ) != -1 ){
      switch( opt ){
         case '?':         // getopt() does not recognize the argument
            printf("CommandLineOptions: ERROR: invalid command line option\n");
            if (optopt == 'n')
               fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            CommandLineOptionsHelp();
            exit( 1 );
         case 'h':
            CommandLineOptionsHelp();
            exit( 0 );
         case 'c':{
            cl->if_generate = false;
            break;
         }
         case 'n':{
            cl->n = atoi(optarg);
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
   printf("vect_trian_mesh - Generate 2D triangel unstructed mesh based on rectangles.\n");
   printf("\n");
   printf("Usage:\n");
   printf("   ./vect_trian_mesh [OPTION] FILENAME\n");
   printf("Options:\n");
   printf("   -c   -- only check mesh validity. (Defualt is do both)\n");
   printf("   -p   -- only generate\n");
   printf("   -n N -- Set number of boundary element on each edge. N must be even\n");
   printf("   -h   -- show help information\n");
   printf("Note:\n");
   printf("    See help info of \"convert\" for file standards.\n");
}
