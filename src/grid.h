/*
  File:      grid.h
  Purpose:   管理计算网格

  Author: <kw.xu@hotmail.com>
  Last modified: "2016-05-13 10:36:55 kaiwen"
*/

#ifndef _GRID_H_
#define _GRID_H_

#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include "def.h"

class VECT;
class BOUNDARY;
class ELEMENT;

// ===================================================================
class GRID
{
public:
   GRID() : n_vert(0), n_elem(0){};
   // GRID(const string & elemfile, const string & vertfile);
   // GRID(const string & filename);
   ~GRID();

   // read mesh
   bool ReadGridFile(const std::string & elemfile, const std::string & vertfile);
   bool Precondition();
   // bool GridCheck();
   // void UpdateInfo();
   void ShowInfo();
public:
   INT n_vert;
   INT n_elem;
   // VECT * vert_coord;
   std::vector<VECT> vert_coord;
   // ELEMENT * elem;
   // 注意是指针!
   std::vector<ELEMENT*> elem;
   // BOUNDARY * bdry;
private:
   std::string _elemfile, _vertfile;
   friend class ELEMENT;
   friend class BOUNDARY;
};

// ===================================================================
class ELEMENT
{
public:
   ELEMENT(){};
   // 这种初始化方式只能使用 vector template 进行处理
   ELEMENT(GRID * _grid){
      isbd = false;
      is_ghost = false;
      grid = _grid;
   }
   virtual ~ELEMENT(){};
   INT verts[NVert];
   // INT nhgbs[NEdge];
   // 该单元的编号.
   INT index;
   // 用于标识单元上某条边在这条边相邻的单元中的编号.
   // 每个位置上合理的值只可能是 0, 1, 2.
   // 初始化为 -1.
   // 如果 nghbs[i] < 0
   // INT nghbr_edge_index[NEdge];
   INT nghbs[NEdge];
   FLOAT area;
   VECT barycenter;
   VECT orient[NEdge];
   FLOAT edge_len[NEdge];
// 是否位于边界, 这里指的不是 ghost cell.
   bool isbd;
   // ghost cell
   bool is_ghost;
private:
   GRID * grid;
public:
   // 计算(更新) element 的成员变量, 如 area, bary, oreint, etc.
   // bool UpdateElemInfo(const GRID * _grid);
   bool UpdateElemInfo();
   // 获取第 i 个顶点的位置
   const VECT & GetVertPosition(INT i){
      return grid->vert_coord[verts[i]];
   }
   // 获得第 i 条边的中点位置.
   const VECT & GetEdgeMidPosition(INT i){
      return 0.5 * (GetVertPosition((i+1)%NEdge) + GetVertPosition((i+2)%NEdge));
   }
   // 返回邻居单元的指针
   ELEMENT * GetNeighborElem(INT i){
      return grid->elem[nghbs[i]];
   }
   // 返回第 i 边在邻居单元上的编号.
   INT GetNghbrEdgeIndex(INT i);
   GRID * GetGrid(){
      return grid;
   }
   // void ShowElemInfo();
   // 正交坐标系到重心坐标的转换
   void XY2Lambda(BARY_COOD & bary_cood,const VECT &vc);
   // ABSS _XY2Lambda(const VECT &vc);
public:
   // 直角坐标与重心坐标的转换所用变量
   FLOAT omega[NVert];
   FLOAT eta[NVert];
   FLOAT xi[NVert];
   // std::vector<std::vector<FLOAT, NVert>, Dim> elem_jacobian;
   // std::vector<BARY_COOD, Dim> elem_jacobian;
   std::array<BARY_COOD, Dim> elem_jacobian;
};

#endif /* _GRID_H_ */
