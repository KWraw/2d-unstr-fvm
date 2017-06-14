/*
  File:      convert.c
  Last modified: "2015-05-10 21:57:53 kaiwen"
  Author: Xu KaiWen <xukaiwen@lsec.cc.ac.cn>

  用法:

  处理 MSH (gmsh mesh format). 2D.

  具体用法参考 help 文档:
  ./convert -h

  输入:
  ./convert [OPTIONS] filename

  输出:
  - meshname.n
  [节点的个数]
  [节点编号] [节点位置] [节点的标号(用于确定是否在边界以及在那个边界)]
  - meshname.e
  [单元的个数]
  [单元编号] [单元的节点编号] [相邻单元(按节点顺序)的编号] [单元标记
  (确定是否在边界上)]
*/
#include "grid.h"
#include "timer.h"
#include "cloption.h"

int
main(int argc, char *argv[]){
   char filename[100],
      input_file_gmsh[100],
      out_file_node[100],
      out_file_element[100];

   double start, finish, elapsed;
   FILE *fp;

   /* 处理命令行参数 */
   /* if(argc<2){ */
   /*    printf("Error: no input meshfile! e.g. : ./convert double-mach-real\n"); */
   /*    exit(1); */
   /* } */
   CommandLineOptions cl;
   cl.if_check = true;
   cl.if_convert = true;
   GetCommandLineOptions(argc, argv, &cl);
   strcpy(filename,cl.filename);
   strcpy(input_file_gmsh,filename);
   strcat(input_file_gmsh,".msh");
   strcpy(out_file_node,filename);
   strcat(out_file_node,".n");
   strcpy(out_file_element,filename);
   strcat(out_file_element,".e");

   if (cl.if_convert){
      printf("========================\n");
      printf("Start converting mesh...\n");
      printf("Input file is : %s\n",input_file_gmsh);
      GRID *g = NewGrid(-1);
      GET_TIME(start);
      ReadMSHFile(input_file_gmsh, g);
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
