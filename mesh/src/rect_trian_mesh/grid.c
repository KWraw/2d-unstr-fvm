/*
  File:      grid.c

  Author: Xu KaiWen <xukaiwen@lsec.cc.ac.cn>
  Last modified: "2016-04-29 01:13:39 kaiwen"
*/

#include "grid.h"
#include "timer.h"
#include "cloption.h"

/* static void */
/* PreconditMesh(GRID * g); */

GRID*
NewGrid(INT flag)
{
   GRID *g = (GRID*)malloc(sizeof(GRID));
   g->verts = NULL;
   g->vm = NULL;
   g->en = NULL;
   g->ee = NULL;
   g->em = NULL;

   if(flag==-1){
      return g;
   }

   printf("Grid Error!\n");
   exit(-1);
}

int
CheckMesh(char * nodefile, char * elemfile){
   GRID * g = NewGrid(-1);
   INT n,m,idx,k,l;
   FILE *fp;

   printf("Reading node file : %s ... ",nodefile);
   fp=fopen(nodefile,"r");
   if(!fp){
      printf("Nodes file is not found!\n");
      exit(-1);
   }

   fscanf(fp,"%d",&g->nvert);
   g->verts=(COORD*)malloc(sizeof(COORD)*g->nvert);
   g->vm   =(INT*)malloc(sizeof(INT)*g->nvert);
   for(n=0;n<g->nvert;n++){
      fscanf(fp,"%d: %le %le %d",&idx,&g->verts[n][0],&g->verts[n][1],&g->vm[n]);
   }

   fclose(fp);
   printf("Over!\n");

   printf("Reading element file : %s ... ",elemfile);
   fp=fopen(elemfile,"r");
   if(!fp){
      printf("Elements file is not found!\n");
      exit(-1);
   }
   fscanf(fp,"%d",&g->nelem);


   g->en   =(VERTS*)malloc(sizeof(VERTS)*g->nelem);
   g->ee   =(ELEMS*)malloc(sizeof(ELEMS)*g->nelem);
   g->em   =(INT*)malloc(sizeof(INT)*g->nelem);

   for(n=0;n<g->nelem;n++){
      fscanf(fp,"%d: %d %d %d %d %d %d %d",
             &idx,&g->en[n][0],&g->en[n][1],&g->en[n][2],
             &g->ee[n][0],&g->ee[n][1],&g->ee[n][2],
             &g->em[n]);
      g->ee[n][0] = g->ee[n][0]<0? -1 : g->ee[n][0];
      g->ee[n][1] = g->ee[n][1]<0? -1 : g->ee[n][1];
      g->ee[n][2] = g->ee[n][2]<0? -1 : g->ee[n][2];
   }
   printf("Over!\n");

   /*--------------------+
     | Checking node file |
     +--------------------*/
   printf("\n");
   printf("Mesh information:\n");
   FLOAT len,max_len=0,min_len=10000;
   FLOAT dx,dy;
   INT ii=0,jj=1,kk=2;
   for(n=0;n<g->nelem;n++){
      for(ii=0;ii<NEdge;++ii,jj=(jj+1)%NEdge,kk=(kk+1)%NEdge){
         dx = g->verts[g->en[n][jj]][0] - g->verts[g->en[n][kk]][0];
         dy = g->verts[g->en[n][jj]][1] - g->verts[g->en[n][kk]][1];
         len = sqrt(dx*dx + dy*dy);
         max_len = MAX(max_len, len);
         min_len = MIN(min_len, len);
      }
   }

   /*------------------+
     |  Mesh infomation |
     +------------------*/
   printf ("max_len = %lf\n",max_len);
   printf ("min_len = %lf\n",min_len);
   printf ("max_len/min_len = %lf\n",max_len/min_len);
   printf ("g->nvert = %d\n",g->nvert);
   printf ("g->nelem = %d\n",g->nelem);

   /*------------------------+
     |  Checking element file |
     +------------------------*/
   INT adjt;
   /* check */
   printf("Start checking...\n");
   printf(" -- Checking element neighbour consistency...");
   for(n=0;n<g->nelem;n++){
      for(m=0;m<NEdge;m++){
         adjt=g->ee[n][m];
         if(adjt<0)
            continue;
         l=0;
         for(k=0;k<NEdge;k++){
            if(g->ee[adjt][k]==n)
               l++;
         }
         if(l!=1){
            printf("Wrong!\n");
            return 1;
         }
      }
   }
   printf("Over!\n");

   printf(" -- Checking mesh index orientation...");
   ii=0,jj=1,kk=2;
   FLOAT dx1,dx2,dy1,dy2, delta;
   for(n=0;n<g->nelem;n++){
      dx1 = g->verts[g->en[n][jj]][0] - g->verts[g->en[n][ii]][0];
      dy1 = g->verts[g->en[n][jj]][1] - g->verts[g->en[n][ii]][1];
      dx2 = g->verts[g->en[n][kk]][0] - g->verts[g->en[n][jj]][0];
      dy2 = g->verts[g->en[n][kk]][1] - g->verts[g->en[n][jj]][1];
      delta = dx1*dy2 - dx2*dy1;
      if(delta < 0){
         printf("\nCheckMesh Error, wrong orientation element.\n");
         exit(1);
      }
   }
   printf("Over!\n");

   fclose(fp);
   printf("Finished!\n");
   return 0;
}

void
ExportMesh(char * nodefile, char * elementfile, GRID * g){
   INT n;
   FILE * fp;
   fp=fopen(nodefile,"w");
   printf("Write results to file %s ... \n",nodefile);
   /* 节点: 节点编号: 空间坐标 标志 */
   fprintf(fp,"%d\n",g->nvert);
   for ( n = 0; n < g->nvert; ++n){
      fprintf(fp,"%6d: %25.16e %25.16e %6d\n", n,g->verts[n][0],
              g->verts[n][1],g->vm[n]);
   }
   fclose(fp);

   fp=fopen(elementfile,"w");
   printf("Write results to file %s ... \n",elementfile);
   /* 面单元: 单元编号: 三个节点编号 相邻单元编号 标志 */
   fprintf(fp,"%d\n",g->nelem);
   for ( n = 0; n < g->nelem; ++n){
      fprintf(fp,"%6d: %6d %6d %6d %6d %6d %6d %6d\n", n,
              g->en[n][0],g->en[n][1],g->en[n][2],
              g->ee[n][0],
              g->ee[n][1],
              g->ee[n][2],
              g->em[n]);
   }
   fclose(fp);

   return;
}

#define EDGE_LEN 1.0L
static void
PreconditMesh(GRID * g);
void
RectTrianMeshGenerateTypeC(GRID * g, int n){
   /* check n */
   if(n%2){
      printf ("RectTrianMeshGenerate, ERROR! n must be a odd integer.\n");
      CommandLineOptionsHelp();
      exit(1);
   }

   INT n_node = (n+1) * (n+1);
   INT n_elem = n * n * 2;
   FLOAT   dd = EDGE_LEN / n;
   int i,j,k,e;
   /* float x,y; */

   /*----------+
     | 申请空间 |
     +----------*/
   g->nvert = n_node;
   g->verts = (COORD *)malloc(sizeof(COORD) * n_node);
   g->vm    = (INT *)malloc(sizeof(INT) * n_node);
   g->nelem = n_elem;
   g->en    = (VERTS *)malloc(sizeof(VERTS) * n_elem);
   g->ee    = (ELEMS *)malloc(sizeof(ELEMS) * n_elem);
   g->em    = (INT *)malloc(sizeof(INT) * n_elem);

   /*--------+
     | 初始化 |
     +--------*/
   for(k = 0; k<n_node; k++){
      g->verts[k][0] = 0.0;
      g->verts[k][1] = 0.0;
      g->vm[k] = -1;
   }
   for(e=0;e<n_elem;++e){
      g->em[e]=0;
      for(i=0;i<NEdge;i++){
         g->ee[e][i]=-1000;     /* 默认值, 观察是否有边没有被标识. */
         g->en[e][i]=-1000;
      }
   }

   /*--------------+
     | 计算节点信息 |
     +--------------*/
   for(k = 0; k<n_node; k++){
      /* x 方向标号 */
      i = k % (n+1);
      /* y */
      j = k / (n+1);
      /* 坐标 */
      g->verts[k][0] = i * dd;
      g->verts[k][1] = j * dd;
      /* 边界 mark, 注意初始化时赋值为 -1 */
      if((i*(n-i) + j*(n-j)) == 0){
         /* 位于四个角点之一 */
         g->vm[k] = INTPTMK;
      } else if (i == 0){g->vm[k] = 1;
      } else if (i == n){g->vm[k] = 3;
      } else if (j == 0){g->vm[k] = 2;
      } else if (j == n){g->vm[k] = 4;
      } else {g->vm[k] = 0;}
   }

   /*--------------+
     | 计算单元信息 |
     +--------------*/
   INT elem_type;
   INT base_vert;
   for(e = 0; e<n_elem; e++){
      /* 单元类型 */
      /* elem_type = e % 4; */
      elem_type = (e + (e / (2 * n)) * 2) % 4;
      /* 节点的编号                     */
      /* 单元所在正方形左下角节点的编号 */
      base_vert = e / 2;
      base_vert = base_vert + base_vert / n;
      switch(elem_type){
         case 0 :{
            g->en[e][0] = base_vert;
            g->en[e][1] = base_vert+1;
            g->en[e][2] = base_vert+n+1;
            break;
         }
         case 1 :{
            g->en[e][0] = base_vert+1;
            g->en[e][1] = base_vert+n+2;
            g->en[e][2] = base_vert+n+1;
            break;
         }
         case 2 :{
            g->en[e][0] = base_vert;
            g->en[e][1] = base_vert+n+2;
            g->en[e][2] = base_vert+n+1;
            break;
         }
         case 3 :{
            g->en[e][0] = base_vert;
            g->en[e][1] = base_vert+1;
            g->en[e][2] = base_vert+n+2;
            break;
         }
         default:{
            printf ("RectTrianMeshGenerate, ERROR! Unkown elem type.\n");
            exit(1);
         }
      }
   } /* each elem */

   /*---------------------------------------------------+
     | 通过单元与节点间的对应关系计算邻居单元与边界信息. |
     +---------------------------------------------------*/
   PreconditMesh(g);
   return;
}

void
RectTrianMeshGenerateTypeA(GRID * g, int n){
   /* check n */
   if(n%2){
      printf ("RectTrianMeshGenerate, ERROR! n must be a odd integer.\n");
      CommandLineOptionsHelp();
      exit(1);
   }

   INT n_node = (n+1) * (n+1);
   INT n_elem = n * n * 2;
   FLOAT   dd = EDGE_LEN / n;
   int i,j,k,e;
   /* float x,y; */

   /*----------+
    | 申请空间 |
    +----------*/
   g->nvert = n_node;
   g->verts = (COORD *)malloc(sizeof(COORD) * n_node);
   g->vm    = (INT *)malloc(sizeof(INT) * n_node);
   g->nelem = n_elem;
   g->en    = (VERTS *)malloc(sizeof(VERTS) * n_elem);
   g->ee    = (ELEMS *)malloc(sizeof(ELEMS) * n_elem);
   g->em    = (INT *)malloc(sizeof(INT) * n_elem);

   /*--------+
     | 初始化 |
     +--------*/
   for(k = 0; k<n_node; k++){
      g->verts[k][0] = 0.0;
      g->verts[k][1] = 0.0;
      g->vm[k] = -1;
   }
   for(e=0;e<n_elem;++e){
      g->em[e]=0;
      for(i=0;i<NEdge;i++){
         g->ee[e][i]=-1000;     /* 默认值, 观察是否有边没有被标识. */
         g->en[e][i]=-1000;
      }
   }

   /*--------------+
    | 计算节点信息 |
    +--------------*/
   for(k = 0 ; k < n_node ; k ++){
      i = k % (n+1);
      j = k / (n+1);
      g->verts[k][0] = i * dd;
      g->verts[k][1] = j * dd;
      if((i*(n-i) + j*(n-j)) == 0){
         /* 位于四个角点之一, 若不然则位于四条边 */
         g->vm[k] = INTPTMK;
      } else if (i == 0){g->vm[k] = 1;
      } else if (i == n){g->vm[k] = 3;
      } else if (j == 0){g->vm[k] = 2;
      } else if (j == n){g->vm[k] = 4;
      } else {g->vm[k] = 0;}
   }

   /*-------------------------------------------------------------+
    | 计算单元信息, 给出 en 即可, 其他的操作在 PreconditMesh 里面 |
    +-------------------------------------------------------------*/
   /* en: index of vertex */
   INT elem_type;
   INT base_vert;
   for(e = 0; e < n_elem; e++){
      elem_type = e % 2;
      base_vert = e / 2;
      base_vert = base_vert + base_vert / n;
      switch(elem_type){
         case 0 :{
            g->en[e][0] = base_vert;
            g->en[e][1] = base_vert + 1;
            g->en[e][2] = base_vert + n + 2;
            break;
         }
         case 1 :{
            g->en[e][0] = base_vert;
            g->en[e][1] = base_vert + n + 2;
            g->en[e][2] = base_vert + n + 1;
            break;
         }
         default:{
            printf ("RectTrianMeshGenerate, ERROR! Unkown elem type.\n");
            exit(1);
         }
      }
   }

   PreconditMesh(g);
   return;
}

void
RectTrianMeshGenerateTypeB(GRID * g, int n){
   /* check n */
   if(n%2){
      printf ("RectTrianMeshGenerate, ERROR! n must be a odd integer.\n");
      CommandLineOptionsHelp();
      exit(1);
   }

   INT n_node = (n+1) * (n+1);
   INT n_elem = n * n * 2;
   FLOAT   dd = EDGE_LEN / n;
   int i,j,k,e;
   /* float x,y; */

   /*----------+
    | 申请空间 |
    +----------*/
   g->nvert = n_node;
   g->verts = (COORD *)malloc(sizeof(COORD) * n_node);
   g->vm    = (INT *)malloc(sizeof(INT) * n_node);
   g->nelem = n_elem;
   g->en    = (VERTS *)malloc(sizeof(VERTS) * n_elem);
   g->ee    = (ELEMS *)malloc(sizeof(ELEMS) * n_elem);
   g->em    = (INT *)malloc(sizeof(INT) * n_elem);

   /*--------+
     | 初始化 |
     +--------*/
   for(k = 0; k<n_node; k++){
      g->verts[k][0] = 0.0;
      g->verts[k][1] = 0.0;
      g->vm[k] = -1;
   }
   for(e=0;e<n_elem;++e){
      g->em[e]=0;
      for(i=0;i<NEdge;i++){
         g->ee[e][i]=-1000;     /* 默认值, 观察是否有边没有被标识. */
         g->en[e][i]=-1000;
      }
   }

   /*--------------+
    | 计算节点信息 |
    +--------------*/
   for(k = 0 ; k < n_node ; k ++){
      i = k % (n+1);
      j = k / (n+1);
      g->verts[k][0] = i * dd;
      g->verts[k][1] = j * dd;
      if((i*(n-i) + j*(n-j)) == 0){
         /* 位于四个角点之一, 若不然则位于四条边 */
         g->vm[k] = INTPTMK;
      } else if (i == 0){g->vm[k] = 1;
      } else if (i == n){g->vm[k] = 3;
      } else if (j == 0){g->vm[k] = 2;
      } else if (j == n){g->vm[k] = 4;
      } else {g->vm[k] = 0;}
   }

   /*-------------------------------------------------------------+
    | 计算单元信息, 给出 en 即可, 其他的操作在 PreconditMesh 里面 |
    +-------------------------------------------------------------*/
   /* en: index of vertex */
   INT elem_type;
   INT base_vert;
   for(e = 0; e < n_elem; e++){
      elem_type = e % 2;
      base_vert = e / 2;
      base_vert = base_vert + base_vert / n;
      switch(elem_type){
         case 0 :{
            g->en[e][0] = base_vert;
            g->en[e][1] = base_vert + 1;
            g->en[e][2] = base_vert + n + 1;
            break;
         }
         case 1 :{
            g->en[e][0] = base_vert + 1;
            g->en[e][1] = base_vert + n + 2;
            g->en[e][2] = base_vert + n + 1;
            break;
         }
         default:{
            printf ("RectTrianMeshGenerate, ERROR! Unkown elem type.\n");
            exit(1);
         }
      }
   }

   PreconditMesh(g);
   return;
}




/* 找单元边界, 给单元标号 */
static void
PreconditMesh(GRID * g){
   INT mark,i,m,ii,jj,kk,l,e;
   INT mark_jj, mark_kk;
   ii=0,jj=1,kk=2;
   for (e =0; e <g->nelem; ++e){
      for (ii=0;ii<NEdge;++ii,jj=(jj+1)%NEdge,kk=(kk+1)%NEdge){
         /* /\\* 先判断 ii 对面的边是不是边界, 当两个点都位于边界这条边   *\\/ */
         /* /\\* 是边界. 边界上的点有非零的 mark 并且不同的 mark 值表示了 *\\/ */
         /* /\\* 不同的边界. 为了不至于引起混淆, 用边界 tag 值的负值做为  *\\/ */
         /* /\\* 标识. 这个标识占用相邻单元编号的位置.                    *\\/ */
         mark_jj = g->vm[g->en[e][jj]];
         mark_kk = g->vm[g->en[e][kk]];
         if((mark_jj!=0)&&(mark_kk!=0)){
            if((mark_jj == mark_kk) || (mark_jj == INTPTMK || mark_kk == INTPTMK)){
               /* 两个顶点位于同一条边上, 边位于边界. */
               g->ee[e][ii]=-MAX(g->vm[g->en[e][jj]],g->vm[g->en[e][kk]]);
               g->em[e]=-1;
               continue;
            }
            /* else{ */
               /* 两个顶点位于不同的边上, 此时 */
            /* } */
         }
         /* /\\* 不是边界单元, 找对面的单元 *\\/ */
         for(m=0;m<g->nelem;++m){
            l=0;
            if (m!=e){
               for ( i = 0; i < NEdge; ++i){
                  if (g->en[m][i]==g->en[e][jj]){
                     l++;
                  } else if(g->en[m][i]==g->en[e][kk]){
                     l++;
                  }
               }
            }
            if(l==2){           /* 说明是 ii 对面的相邻单元 */
               g->ee[e][ii]=m;
            }
         }
      }
   }
   /* INT mark, n,m,i,ii,jj,kk,l; */
   /* ii=0,jj=1,kk=2; */
   /* for (n = 1; n < g->nelem+1; ++n){ */
   /*    for (ii=0;ii<NEdge;++ii,jj=(jj+1)%NEdge,kk=(kk+1)%NEdge){ */
   /*       /\* 先判断 ii 对面的边是不是边界, 当两个点都位于边界这条边   *\/ */
   /*       /\* 是边界. 边界上的点有非零的 mark 并且不同的 mark 值表示了 *\/ */
   /*       /\* 不同的边界. 为了不至于引起混淆, 用边界 tag 值的负值做为  *\/ */
   /*       /\* 标识. 这个标识占用相邻单元编号的位置.                    *\/ */
   /*       if((g->vm[g->en[n][jj]]!=0)&&(g->vm[g->en[n][kk]]!=0)){ */
   /*          g->ee[n][ii]=-MAX(g->vm[g->en[n][jj]],g->vm[g->en[n][kk]]); */
   /*          g->em[n]=-1; */
   /*          continue; */
   /*       } */
   /*       /\* 不是边界单元, 找对面的单元 *\/ */
   /*       for(m=1;m<(g->nelem+1);++m){ */
   /*          l=0; */
   /*          if (m!=n){ */
   /*             for ( i = 0; i < NEdge; ++i){ */
   /*                if (g->en[m][i]==g->en[n][jj]){ */
   /*                   l++; */
   /*                } else if(g->en[m][i]==g->en[n][kk]){ */
   /*                   l++; */
   /*                } */
   /*             } */
   /*          } */
   /*          if(l==2){           /\* 说明是 ii 对面的相邻单元 *\/ */
   /*             g->ee[n][ii]=m; */
   /*          } */
   /*       } */
   /*    } */
   /* } */
   return;
}

/* void ReadMSHFile(char * filename, GRID * g){ */
/*    char charcache[100]; */
/*    INT i, n,vn,vn1,vn2,len,intcache; */
/*    INT intarray[10]; */
/*    FLOAT xx,yy,floatchache; */
/*    INT num_total_elements, num_point_elements, num_line_elements, */
/*       num_face_elements, l, lmark; */

/*    FILE *fp; */
/*    if(!filename){ */
/*       printf("Please input mesh file name!\n"); */
/*       exit(-1); */
/*    } */
/*    printf("Reading mesh file : %s ... \n",filename); */
/*    fp=fopen(filename,"r"); */
/*    if(!fp){ */
/*       printf("Mesh file is not found!\n"); */
/*       exit(-2); */
/*    } */

/*    while(1){ */
/*       fscanf(fp,"%s",charcache); */
/*       /\*-----------------------------------------------------------------------+ */
/*         |  读 nodes, 注意: GMSH 文件点的编号是从 1 开始的, 预处理网格保留这个设 | */
/*         |  置, 在输出的时候再修改为从 0 开始编号.                               | */
/*         +-----------------------------------------------------------------------*\/ */
/*       if (0==strcmp(charcache,"$Nodes")){ */
/*          fscanf(fp,"%d",&g->nvert); */
/*          g->verts=(COORD*)malloc(sizeof(COORD)*(g->nvert+1)); */
/*          g->vm =(INT*)malloc(sizeof(INT)*(g->nvert+1)); */
/*          for (n = 1; n < g->nvert+1; ++n){ */
/*             fscanf(fp,"%d %lf %lf %lf",&intcache,&xx,&yy,&floatchache); */
/*             vn=intcache; */
/*             g->verts[vn][0]=xx; */
/*             g->verts[vn][1]=yy; */
/*             g->vm[n]=0; */
/*          } */
/*       } */

/*       /\*------------------------------------------+ */
/*         | 读 element, 同样注意编号开始是 1 的问题. | */
/*         +------------------------------------------*\/ */
/*       if (0==strcmp(charcache,"$Elements")){ */
/*          /\* num_elements is the total number of all kinds of       *\/ */
/*          /\* element. Including vertics element, line element, face *\/ */
/*          /\* element.                                               *\/ */
/*          fscanf(fp,"%d",&num_total_elements); */

/*          /\*------------------------------------+ */
/*            | 处理点线单元信息, 记录面单元的个数 | */
/*            +------------------------------------*\/ */
/*          num_face_elements=0; */
/*          for (n=0;n<num_total_elements;++n){ */
/*             while(1){ */
/*                /\* 格式: 单元编号 类型 标记的个数(?) 第一个标记 第二 *\/ */
/*                /\* 个 节点编号(点线面单元分别有 1-3 个)              *\/ */
/*                fscanf(fp,"%d %d",&intarray[0],&intarray[1]); */

/*                if (15==intarray[1]){ /\* 点单元 *\/ */
/*                   len=6; */
/*                   for (i=2; i<len; ++i) */
/*                      fscanf(fp,"%d",&intarray[i]); */
/*                   /\*  读取点单元之后的操作. *\/ */
/*                   /\* 每个点单元的点一定是位于两个边界上(注意定义 geo */
/*                    * 文件的时候要定义 Physical Point), 给 vm 取 */
/*                    * 为 -100 (INTPTMK), 做为标识. *\/ */
/*                   vn=intarray[5]; */
/*                   g->vm[vn]=INTPTMK; */
/*                   break; */
/*                } else if (1==intarray[1]){ /\* 线单元 *\/ */
/*                   len=7; */
/*                   for (i=2; i<len; ++i) */
/*                      fscanf(fp,"%d",&intarray[i]); */
/*                   /\*                                                   *\/ */
/*                   /\* intarray[3] 是线单元所在边界的 Physical entity    *\/ */
/*                   /\* number, 一般取成正值. 但是边界点 mark 位一般取    *\/ */
/*                   /\* 负值.  /\* Sun May 10 17:34:31 2015 更改这里的策略,处 *\\/ *\/ */
/*                   /\* /\\* 理成正值会更方便一些.                             *\\/ *\/ */
/*                   /\* lmark=-intarray[3]; *\/ */
/*                   lmark=intarray[3]; */
/*                   vn1=intarray[5]; */
/*                   vn2=intarray[6]; */
/*                   if (INTPTMK!=g->vm[vn1]) /\* 如果不是拐点. *\/ */
/*                      g->vm[vn1]=lmark; */
/*                   if (INTPTMK!=g->vm[vn2]) */
/*                      g->vm[vn2]=lmark; */
/*                   break; */
/*                } else if (2==intarray[1]){ /\* 面单元 *\/ */
/*                   /\* 这部分只统计面单元的个数. *\/ */
/*                   len=8; */
/*                   for (i=2; i<len; ++i) */
/*                      fscanf(fp,"%d",&intarray[i]); */
/*                   num_face_elements++; */
/*                   break; */
/*                } else { */
/*                   printf("Error, unknown element type : %d\n",intarray[1]); */
/*                   exit(-3); */
/*                } */
/*             } */
/*          } /\* num_total_elements for loop. *\/ */
/*          break; */
/*       } /\* if $Elements *\/ */
/*    } /\* while(1) *\/ */

/*    /\*----------------+ */
/*      | 处理面单元信息 | */
/*      +----------------*\/ */
/*    g->nelem=num_face_elements; */
/*    g->en=(VERTS*)malloc(sizeof(VERTS)*(g->nelem+1)); */
/*    g->ee=(ELEMS*)malloc(sizeof(ELEMS)*(g->nelem+1)); */
/*    g->em=(INT*)malloc(sizeof(INT)*(g->nelem+1)); */

/*    /\* 初始化 *\/ */
/*    for(n=1;n<(g->nelem+1);++n){ */
/*       g->em[n]=0; */
/*       for(i=0;i<NEdge;i++){ */
/*          g->ee[n][i]=-1000;     /\* 默认值, 观察是否有边没有被标识. *\/ */
/*          g->en[n][i]=-1000; */
/*       } */
/*    } */

/*    fseek(fp,0,SEEK_SET); */
/*    INT nelem; */
/*    INT residue=num_total_elements-num_face_elements; */
/*    /\* 这一次循环就是为了处理面单元. *\/ */
/*    while(1){ */
/*       fscanf(fp,"%s",charcache); */
/*       if (0==strcmp(charcache,"$Elements")){ */
/*          for (n=0;n<num_total_elements;++n){ */
/*             while(1){ */
/*                fscanf(fp,"%d %d",&intarray[0],&intarray[1]); */
/*                if (15==intarray[1]){ /\* 点单元 *\/ */
/*                   len=6; */
/*                   for (i=2; i<len; ++i) */
/*                      fscanf(fp,"%d",&intarray[i]); */
/*                   break; */
/*                } else if (1==intarray[1]){ /\* 线单元 *\/ */
/*                   len=7; */
/*                   for (i=2; i<len; ++i) */
/*                      fscanf(fp,"%d",&intarray[i]); */
/*                   break; */
/*                } else if (2==intarray[1]){ /\* 面单元 *\/ */
/*                   len=8; */
/*                   for (i=2; i<len; ++i) */
/*                      fscanf(fp,"%d",&intarray[i]); */
/*                   nelem=intarray[0]-residue; */
/*                   /\* 单元顶点 *\/ */
/*                   g->en[nelem][0]=intarray[5]; */
/*                   g->en[nelem][1]=intarray[6]; */
/*                   g->en[nelem][2]=intarray[7]; */
/*                   break; */
/*                } else { */
/*                   printf("Error, unknown element type : %d\n",intarray[1]); */
/*                   exit(-3); */
/*                } */
/*             } */
/*          } /\* num_total_elements for loop. *\/ */
/*          break; */
/*       } */
/*    } */

/*    fclose(fp); */
/*    printf("Finished!\n"); */

/*    printf("Precondition mesh ...\n"); */
/*    PreconditMesh(g); */
/*    printf("Finished!\n"); */
/*    return; */
/* } */
