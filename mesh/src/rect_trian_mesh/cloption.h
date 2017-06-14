/*
  File:      cloption.h
  Purpose:   处理 convert 的命令行参数. 参数说明见 convert.c

  Author: Xu KaiWen <xukaiwen@lsec.cc.ac.cn>
  Last modified: "2015-05-27 22:15:09 kaiwen"
*/

#ifndef _CLOPTION_H_
#define _CLOPTION_H_

#include <stdio.h>
#include <string.h>
#include <stdbool.h>

typedef struct
{
   char filename[100];
   bool if_check;
   bool if_generate;
   /* 每条边上线单元的个数 */
   int n;
} CommandLineOptions;

void GetCommandLineOptions(const int argc, char * argv[], CommandLineOptions * cl);
void CommandLineOptionsHelp();


#endif /* _CLOPTION_H_ */
