#ifndef _GRID_H_
#define _GRID_H_


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define Dim   2
#define NVert 3
#define NEdge 3

#define INTPTMK -100              /* 拐点的 mark, 这限制了边界的数量不能
                                   * 超 100  */

#define DEFAULT_NBR -1000       /* 边界编号的初始化值, 为了检测是否出
                                 * 现没有找到邻居单元的情形 */
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

typedef double      FLOAT;
typedef int           INT;

typedef FLOAT COORD[Dim];
typedef INT VERTS[NVert];
typedef INT ELEMS[NEdge];

typedef struct {
   COORD *verts;     //coordinate of verts
   INT   *vm;        /* vertical mark. */
   INT   nvert;      //number of vertex
   INT   nelem;      //number of element

   /* 单元上 */
   VERTS *en;        //index of vert;
   ELEMS *ee;        //index of neighbour element
   INT   *em;        //mark of element
} GRID;

GRID *
NewGrid(INT flag);

void
ReadMSHFile(char * filename, GRID * g);

void
ExportMesh(char * nodefile, char * elementfile, GRID * g);

int
CheckMesh(char * nodefile, char * elemfile);
#endif /* _GRID_H_ */
