/*
  File:      rect_trian_mesh.c
  Purpose:   生成基于对正则网格的三角形网格

  Author: Xu KaiWen <xukaiwen@lsec.cc.ac.cn>
  Last modified: "2016-04-29 01:44:00 kaiwen"
*/

#include "grid.h"
#include "cloption.h"
#include "timer.h"

int
main(int argc, char *argv[]){
   char filename[100],
      input_file_gmsh[100],
      out_file_node[100],
      out_file_element[100];

   double start, finish, elapsed;
   FILE *fp;
   CommandLineOptions cl;
   /* 命令行参数默认值 */
   cl.if_check = true;
   cl.if_generate = true;
   cl.n = 10;
   GetCommandLineOptions(argc, argv, &cl);
   strcpy(filename,cl.filename);
   strcpy(input_file_gmsh,filename);
   strcat(input_file_gmsh,".msh");
   strcpy(out_file_node,filename);
   strcat(out_file_node,".n");
   strcpy(out_file_element,filename);
   strcat(out_file_element,".e");

   if (cl.if_generate){
      printf("========================\n");
      printf("Start generating regular unstructed mesh...\n");
      /* printf("Input file is : %s\n",input_file_gmsh); */
      GRID *g = NewGrid(-1);
      GET_TIME(start);
      /* ReadMSHFile(input_file_gmsh, g); */
      /* RectTrianMeshGenerate(g, cl.n); */
      /* RectTrianMeshGenerateTypeA(g, cl.n); */
      /* RectTrianMeshGenerateTypeB(g, cl.n); */
      RectTrianMeshGenerateTypeC(g, cl.n);
      GET_TIME(finish);
      elapsed = finish - start;
      printf ("(%f s)\n",elapsed);
      printf("Save nodes information into file : %s\n",out_file_node);
      printf("Save element information into file : %s\n",out_file_element);
      ExportMesh(out_file_node,out_file_element,g);
   }
   if (cl.if_check){
      printf("======================\n");
      printf("Check mesh...\n");
      CheckMesh(out_file_node,out_file_element);
   }

   return 0;
} /* main */
