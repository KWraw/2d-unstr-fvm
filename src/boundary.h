/*
  File:      boundary.h
  Purpose:   边界条件相关的处理 (暂时只处理在区域 [0,1]\times[0,1] 上的周期边界条件.)
             (可以通过新的继承类实现更多的边界条件类型.)
  Author: <kw.xu@hotmail.com>
  Last modified: "2016-05-13 10:39:57 kaiwen"
*/

#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "grid.h"


//////////////////////////////////////////////////////////////////////////////////////
// 注意 : 当边界 ghost 单元的类型为 PERIODIC, boundary_type 是无效的.
// 不需要指定边界条件的类型. 此时 nghbs 位的值(负值), 标识了边界位于哪条边.
// 四条边的标号标准如下:
//                        -----4-----
//                        |         |
//                        1         3
//                        |         |
//                        -----2-----
// 默认变长为 PERIODIC_BOUNDARY::edge_length = 1.0 . 长宽相等. (gmsh 生成的时候就要调好.)
//////////////////////////////////////////////////////////////////////////////////////

// enum BOUNDARY_TYPE
// {DIRICHILET=1, EXTRAPOLATE=2, CHARACTERISTIC=3, PERIODIC=4};
enum BOUNDARY_TYPE
{DIRICHILET=1, EXTRAPOLATE=2, CHARACTERISTIC=3, REFLECTION=4};

enum PERIODIC_EDGE
{PERIODIC_LEFT=1, PERIODIC_DOWN=2, PERIODIC_RIGHT=3, PERIODIC_UP=4};
// 两种生成 ghost cell 的方式.
// MIRROR, 镜像对称.
// PERIODIC, 周期的边界(只适用于矩形的网格)
//           在处理对流方程的时候就默认了计算区域是矩形, 并且网格本身允许进行周期拓展.
enum GHOST_CELL_TYPE
{MIRROR=1, PERIODIC=2};

const FLOAT EDGE_LENGTH = 1.0;

class BOUNDARY_ELEMENT;
class PERIODIC_BOUNDARY_ELEMENT;


// ===================================================================
class BOUNDARY_ELEMENT : public ELEMENT
{
public:
   BOUNDARY_ELEMENT(){};
   BOUNDARY_ELEMENT(GRID * _grid) : ELEMENT(_grid), original_elem(-1){
      is_ghost = true;
   }
   virtual ~BOUNDARY_ELEMENT(){};
   // 边界类型, 如果需要指定的话.
   BOUNDARY_TYPE bd_type;
   // 相邻的区域内单元的编号.
   INT adjust_elem;
   // ghost 单元突出顶点对面的那个顶点.
   INT internal_vert;
   // ghost cell 所对应的单元的编号.
   INT original_elem;
   // 凸出去那个点对应的原来的点, 的全局编号.
   INT original_vert;
   // 返回 ghost 单元的第 i 个顶点在原始单元上的顶点编号.
   virtual INT GetOriginalVertIndex(INT i) = 0;
};

class PERIODIC_BOUNDARY_ELEMENT : public BOUNDARY_ELEMENT
{
public:
   PERIODIC_BOUNDARY_ELEMENT(){};
   PERIODIC_BOUNDARY_ELEMENT(GRID * _grid)
      : BOUNDARY_ELEMENT(_grid){};
   ~PERIODIC_BOUNDARY_ELEMENT(){};
   // 返回 ghost 单元的第 i 个顶点在原始单元上的顶点编号.
   virtual INT GetOriginalVertIndex(INT i);
   // 位于边界的边的中点坐标.
   VECT mid_pt;
   // 位于哪条边. PERIODIC_EDGE 是一个 enum 类型.
   PERIODIC_EDGE which_edge;
};

// ===================================================================
class BOUNDARY
{
// private:
//    GRID * _grid;
public:
   GRID * _grid;
   // 如果不关联一个 grid 的话就没意义了, 所以只设定这样一种构造方式.
   BOUNDARY(GRID * grid) : _grid(grid) {};
   // 新生成的 elem 仍然可以通过 GRID 的 destructor 释放内存.
   // 所以这里不做任何操作.
   virtual ~BOUNDARY(){};
   // 边界线单元的数量.
   // 注意在生成 BOUNDARY 对象的时候不知道边界的个数.
   INT n_boundary;
   // 新生成的边界单元编号的集合. 用于对所有边界单元的直接遍历.
   std::vector<BOUNDARY_ELEMENT*> boundary_elem;
public:
   virtual void ShowInfo() = 0;
   virtual bool ExpendBoundary() = 0; // pure virtual memeber function.
   // 返回全局编号对应的 boundary 编号.
   virtual INT GetBoundaryIndex(INT nghbs_index) = 0;
   // 返回 boundary 编号对应的全局编号.
   virtual INT GetBoundaryIndexVersus(INT boundary_index) = 0;
   // 返回角点的 index, 这只在 periodic boundary 的时候有用.
   virtual INT GetConnerIndex(INT k) = 0;
// private:
   // INT n_dirichilet_bd, n_extrapolate_bd, n_characteristric_bd, n_periodic_bd;
};

// ===================================================================
class PERIODIC_BOUNDARY : public BOUNDARY{
public:
   PERIODIC_BOUNDARY(GRID * grid) : BOUNDARY(grid), edge_length(EDGE_LENGTH){};
   ~PERIODIC_BOUNDARY(){};
   PERIODIC_BOUNDARY_ELEMENT & operator[](int i){
      // return static_cast<PERIODIC_BOUNDARY_ELEMENT *> (*(boundary_elem[i]));
      // return *(static_cast<PERIODIC_BOUNDARY_ELEMENT *> (boundary_elem[i]));
      PERIODIC_BOUNDARY_ELEMENT * periodic_bd_ptr =
         dynamic_cast<PERIODIC_BOUNDARY_ELEMENT *> (boundary_elem[i]);
      // return *(dynamic_cast<PERIODIC_BOUNDARY_ELEMENT *> (boundary_elem[i]));
      if(periodic_bd_ptr)
         return *periodic_bd_ptr;
      else
      {
         std::cout<<
            "PERIODIC_BOUNDARY, ERROR! Can not convert BOUNDARY_ELEMENT* to PERIODIC_BOUNDARY_ELEMENT"<<std::endl;
         exit(1);
      }
   }
   virtual INT GetConnerIndex(INT k){return corner_index[k];}
private:
   std::vector<PERIODIC_BOUNDARY_ELEMENT*> _periodic_left;
   std::vector<PERIODIC_BOUNDARY_ELEMENT*> _periodic_right;
   std::vector<PERIODIC_BOUNDARY_ELEMENT*> _periodic_up;
   std::vector<PERIODIC_BOUNDARY_ELEMENT*> _periodic_down;
   INT _n_left_elems,_n_right_elems,_n_up_elems,_n_down_elems;
   const FLOAT edge_length;
   FLOAT corner_index[4];
public:
   virtual bool ExpendBoundary();
   virtual void ShowInfo();
   virtual INT GetBoundaryIndex(INT nghbs_index);
   virtual INT GetBoundaryIndexVersus(INT boundary_index);
   // bool myfunction (int i,int j) { return (i<j); }
   // 用于对各条边上的 ghost cell 按位置进行排序.
   // 只是对指针 vector 进行排序, 并没有调换对象.
   // boundary_elem 是没有任何改变的.
   friend bool _sortVertical(PERIODIC_BOUNDARY_ELEMENT* ,PERIODIC_BOUNDARY_ELEMENT*);
   friend bool _sortHorizontal(PERIODIC_BOUNDARY_ELEMENT* ,PERIODIC_BOUNDARY_ELEMENT*);
};


#endif /* _BOUNDARY_H_ */
