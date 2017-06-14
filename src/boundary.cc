/*
  File:      boundary.cc
  Purpose:

  Author: <kw.xu@hotmail.com>
  Last modified: "2016-05-13 10:38:55 kaiwen"
*/

#include "boundary.h"
#include <cmath>
#include <vector>
#include <algorithm>

// ===================================================================
// BOUDNARY
// ===================================================================

// ===================================================================
// PERIODIC_BOUNDARY
// ===================================================================

// 值小的排在前面.
bool _sortVertical(PERIODIC_BOUNDARY_ELEMENT* bd_elem1,
                   PERIODIC_BOUNDARY_ELEMENT* bd_elem2){
   return (bd_elem1->mid_pt.y < bd_elem2->mid_pt.y);
}

bool _sortHorizontal(PERIODIC_BOUNDARY_ELEMENT* bd_elem1,
                     PERIODIC_BOUNDARY_ELEMENT* bd_elem2){
   return (bd_elem1->mid_pt.x < bd_elem2->mid_pt.x);
}

bool PERIODIC_BOUNDARY::ExpendBoundary(){
/*
    1. 从 GRID 对象中提取边界相关的信息.
    2. 生成 ghost 节点, 生成 ghost 单元.
    3. 更新 ghost 单元的信息.
    4. 更新原 GRID 对象的信息.

    算法:
    1. 遍历所有单元找到所有边界, 生成 boundary_elem, _periodic_xxxx.
       以及新的 _grid->elem 元素.
       此时可以获得的信息:
       1) 边界所在边.
       2) 边界相邻的单元编号.
       3) 新生成的单元的全局编号.
       4) 需要生成的 ghost 单元的个数.
       5) 对每个 ghost 单元来说, 除去需要新构造的节点之外的其他两个节点.
          默认的顺序如下,
                                2
                               / \
                              / O \
                             0--*--1
                            kk--*--jj
                             \  I /
                              \  /
                               ii
          其中 * 表示位于计算区域边界的那条边.
       6) 利用上一条信息, 可以计算出 * 的中点坐标.
    2. 利用前面所提到的中点信息, 对 _periodic_xxxx 进行排序. 找到对应的内
       部单元.
       1) left, right 对纵坐标进行排序.
       2) up, down 对横坐标进行排序.
       3) 到这一步可以对当前网格是否可以做周期的边界沿拓进行一个判断.
          判断的标准可以按下面几条进行,
          a. left -- right, up -- down 的单元个数是否相同.
          b. 相对应的边单元的中点位置是否对应.
       4) 到这一步, 已经建立了对边上从对应的两个边单元上沿拓出去的 ghost 单
          元的对应关系. 又因为通过前面第一步, ghost 单元中存储了与 ghost单
          元相邻的那个单元的编号. 也就找到了每个 ghost 单元对应的那个内部单
          元.
    3. 完善 ghost 单元.
       1) 从内部单元可以确定出那个位于单元内部的点. 平移 edge_length 作为新
          添 加的节点, push_back 到 vert_coord 后面. (注意与 elem 的处理不
          同, 此时 push_back 的是对象而不是指针.)
       2) 建立 ghost 单元与新节点的联系.
       3) update ghost 单元上的信息.
 */

   INT i,j;
   INT nghb;
   INT ii,jj,kk;
   n_boundary=0;
   // std::vector<INT> internal_vert; // 用于记录内单元的内节点.
   _n_left_elems=0;
   _n_up_elems=0;
   _n_down_elems=0;
   _n_right_elems=0;

   //,-----------------------------
   //| Step 1 -- find the boundary.
   //`-----------------------------

   // if(ghost_type==PERIODIC){
   PERIODIC_BOUNDARY_ELEMENT* periodic_bd; // 一个临时变量
   for(i=0;i<_grid->n_elem;i++){
      for(j=0;j<NEdge;j++){     // 考虑可能会出现一个单元有两条边位于区域边界的情形.
         nghb=_grid->elem[i]->nghbs[j];
         if(nghb<0){
            n_boundary++;
            // create boundary element.
            periodic_bd = new PERIODIC_BOUNDARY_ELEMENT(_grid);
            // boundary_elem 中的存储顺序为 ghost 单元的标准编号.
            boundary_elem.push_back(periodic_bd); // upcasting convert.
            // copy into the GRID elem array for convenience.
            _grid->elem.push_back(periodic_bd); // upcasting convert.
            // adjust internal element
            periodic_bd->adjust_elem=i;
            // 两个已知顶点的标号.
            ii=j;
            jj=(ii+1)%NEdge;
            kk=(jj+1)%NEdge;
            periodic_bd->verts[0] = _grid->elem[i]->verts[kk];
            periodic_bd->verts[1] = _grid->elem[i]->verts[jj];
            // 将内侧顶点的坐标进行缓存, 用于后面新节点的生成.
            // internal_vert 存储的顺序与 boundary_elem
            // 的记录顺序是相同的.
            // internal_vert.push_back(_grid->elem[i]->verts[ii]);
            // 换一种策略,
            // 在 BOUNDARY_ELEMEMT 中添加 internal_vert 成员变量
            // 在 PERIODIC_BOUNDARY_ELEMENT 中添加 original_vert
            periodic_bd->internal_vert=_grid->elem[i]->verts[ii];
            // 边中点的坐标.
            periodic_bd->mid_pt = 0.5 *
               (_grid->vert_coord[periodic_bd->verts[0]] +
                _grid->vert_coord[periodic_bd->verts[1]]);
            // on which domain edge.
            nghb=-nghb;
            switch (nghb){
               case PERIODIC_LEFT :
                  _n_left_elems++;
                  periodic_bd->which_edge=PERIODIC_LEFT;
                  _periodic_left.push_back(periodic_bd);
                  break;
               case PERIODIC_DOWN :
                  _n_down_elems++;
                  periodic_bd->which_edge=PERIODIC_DOWN;
                  _periodic_down.push_back(periodic_bd);
                  break;
               case PERIODIC_RIGHT :
                  _n_right_elems++;
                  periodic_bd->which_edge=PERIODIC_RIGHT;
                  _periodic_right.push_back(periodic_bd);
                  break;
               case PERIODIC_UP :
                  _n_up_elems++;
                  periodic_bd->which_edge=PERIODIC_UP;
                  _periodic_up.push_back(periodic_bd);
                  break;
               default:
                  std::cout<<"ExpendBoundary, Error! Unknown PERIODIC_EDGE type : "
                           <<nghb<<std::endl;
                  exit(1);
            }
            // update internal element's boundary edges neighbour index.
            nghb = _grid->n_elem + n_boundary - 1;
            _grid->elem[i]->nghbs[ii] = nghb;
            // set ghost cell index
            periodic_bd->index = nghb;
         } // if edge is boundary.
      } // each edge.
   } // each element.

   //,-----------------------------------------
   //| Step 2 -- find related internal element.
   //`-----------------------------------------

   // sort by mid point coordinate.
   std::sort (_periodic_left.begin(),_periodic_left.end(),_sortVertical);
   std::sort (_periodic_right.begin(),_periodic_right.end(),_sortVertical);
   std::sort (_periodic_down.begin(),_periodic_down.end(),_sortHorizontal);
   std::sort (_periodic_up.begin(),_periodic_up.end(),_sortHorizontal);

   //,--------------------------------------------------------------
   //| NOTION: Check whether periodic boundary type can implemented.
   //`--------------------------------------------------------------
   if(!((_periodic_left.size()==_periodic_right.size())&&
        (_periodic_down.size()==_periodic_up.size()))){
      std::cout<<"ExpendBoundary, Error! Cannot implement periodic boudnary type."
               <<std::endl
               <<"                       Number of boundary element not equal."
               <<std::endl;
      exit(1);
   }
   std::vector<PERIODIC_BOUNDARY_ELEMENT*>::iterator it1,it2;
   // 比较左右两条边.
   for(it1=_periodic_left.begin(), it2=_periodic_right.begin();
       it1 != _periodic_left.end(); it1++, it2++){
      if(std::fabs((*it1)->mid_pt.y - (*it2)->mid_pt.y) > TOL){
         std::cout<<"ExpendBoundary, Error! Cannot implement periodic boudnary type."
                  <<std::endl
                  <<"                       None-homologous vertical mid-points."
                  <<std::endl;
         exit(1);
      }
   }
   // 比较上下两条边.
   for(it1=_periodic_up.begin(), it2=_periodic_down.begin();
       it1 != _periodic_up.end(); it1++, it2++){
      if(std::fabs((*it1)->mid_pt.x - (*it2)->mid_pt.x) > TOL){
         std::cout<<"ExpendBoundary, Error! Cannot implement periodic boudnary type."
                  <<std::endl
                  <<"                       None-homologous horizontal mid-points."
                  <<std::endl;
         exit(1);
      }
   }

   // build the relationship, between internal element and ghost elements.
   // 1. index of originial elements.
   // 2. vertex.
   for(it1=_periodic_left.begin(), it2=_periodic_right.begin();
       it1 != _periodic_left.end(); it1++, it2++){
      (*it1)->original_elem = (*it2)->adjust_elem;
      (*it1)->original_vert = (*it2)->internal_vert;
      (*it2)->original_elem = (*it1)->adjust_elem;
      (*it2)->original_vert = (*it1)->internal_vert;
   }
   for(it1=_periodic_up.begin(), it2=_periodic_down.begin();
       it1 != _periodic_up.end(); it1++, it2++){
      (*it1)->original_elem = (*it2)->adjust_elem;
      (*it1)->original_vert = (*it2)->internal_vert;
      (*it2)->original_elem = (*it1)->adjust_elem;
      (*it2)->original_vert = (*it1)->internal_vert;
   }

   //,-------------------------------------
   //| Step 3, Complete the ghost elements.
   //`-------------------------------------
   std::vector<BOUNDARY_ELEMENT*>::iterator it_bd;
   // std::vector<VECT>::iterator it_vect;
   PERIODIC_BOUNDARY_ELEMENT* ptr_periodic_bd;
   // INT new_vert, original_vert;
   INT new_vert;
   FLOAT xx,yy;
   for(it_bd=boundary_elem.begin(), new_vert=_grid->n_vert;
       it_bd!=boundary_elem.end(); it_bd++, new_vert++){
      ptr_periodic_bd = static_cast<PERIODIC_BOUNDARY_ELEMENT *> (*it_bd); // downcast
      // original_vert=ptr_periodic_bd->original_vert;
      // location of original vertical.
      xx = _grid->vert_coord[ptr_periodic_bd->original_vert].x;
      yy = _grid->vert_coord[ptr_periodic_bd->original_vert].y;
      // 1. create new vertex.
      // 2. relate elements to new vertex.
      switch(ptr_periodic_bd->which_edge){
         case PERIODIC_LEFT :
            _grid->vert_coord.push_back(VECT((xx - edge_length),yy));
            break;
         case PERIODIC_DOWN :
            _grid->vert_coord.push_back(VECT(xx,(yy - edge_length)));
            break;
         case PERIODIC_RIGHT :
            _grid->vert_coord.push_back(VECT((xx + edge_length),yy));
            break;
         case PERIODIC_UP :
            _grid->vert_coord.push_back(VECT(xx,(yy + edge_length)));
            break;
         default:
            std::cout<<"ExpendBoundary, Error! Unknown PERIODIC_EDGE type."
                     <<std::endl;
            exit(1);
      }
      ptr_periodic_bd->verts[2]=new_vert;
      // update the new element.
      ptr_periodic_bd->UpdateElemInfo();
      // 如何处理 ghost cell 的 nghbs 项. TODO.
   } // each boundary element.

   // } // if ghost_type is PERIODIC
   // else{
   //    std::cout<<"BOUDNARY, Error! MIRROR-type boundary extention method not implement."
   //             <<endl;
   //    exit(1);
   // } // if ghost_type is MIRROR, todo

   //,---------------------
   //| 找到四个角点的 index
   //`---------------------

   // 以之前已经经过排序的 _periodic_XXXX 为基础.
   corner_index[0] = _periodic_left[0]->verts[0];
   corner_index[1] = _periodic_down[_n_down_elems - 1]->verts[0];
   corner_index[2] = _periodic_right[_n_right_elems - 1]->verts[0];
   corner_index[3] = _periodic_up[0]->verts[0];

   return true;
}

void PERIODIC_BOUNDARY::ShowInfo(){
   std::cout<<std::endl
            <<"=========================="<< std::endl
            <<"BOUNDARY.ShowInfo:"<< std::endl
            <<"   ghost type  -- "<< "PERIODIC" << std::endl
            <<"   edge_length -- "<< edge_length  << std::endl
            <<"   n_boundary  -- "<< n_boundary << std::endl
            <<"=========================="<< std::endl;
   // std::vector<PERIODIC_BOUNDARY_ELEMENT*>::iterator it;
   // std::cout<<"Each boundary:"<<std::endl;
   // std::cout<<"left:" << std::endl;
   // for (it = _periodic_left.begin(); it != _periodic_left.end(); ++it)
   //    std::cout<<(*it)->barycenter;
   // std::cout<<"right:" << std::endl;
   // for (it = _periodic_right.begin(); it != _periodic_right.end(); ++it)
   //    std::cout<<(*it)->barycenter;
   // std::cout<<"up:" << std::endl;
   // for (it = _periodic_up.begin(); it != _periodic_up.end(); ++it)
   //    std::cout<<(*it)->barycenter;
   // std::cout<<"down:" << std::endl;
   // for (it = _periodic_down.begin(); it != _periodic_down.end(); ++it)
   //    std::cout<<(*it)->barycenter;
   return;
}

INT PERIODIC_BOUNDARY::GetBoundaryIndex(INT nghbs_index){
   if(nghbs_index < _grid->n_elem){
      // Wrong Message
      std::cout<< "PERIODIC_BOUNDARY, Wrong. can not convert to boundary index."
               << " Exceeding the boundary."
               <<std::endl;
      exit(1);
   } else {
      // convert to boundary index. cf. PERIODIC_BOUNDARY::ExpendBoundary()
      // End of step.1
      return nghbs_index - _grid->n_elem;
   }
}

INT PERIODIC_BOUNDARY::GetBoundaryIndexVersus(INT boundary_index){
   if(boundary_index >= n_boundary){
      // Wrong Message
      std::cout<< "PERIODIC_BOUNDARY, Wrong. can not convert to Grid index."
               << " Exceeding the boundary."
               <<std::endl;
         exit(1);
   }
   return _grid->n_elem + boundary_index;
}

// ===================================================================
// bool GRID::UpdateBoundaryInfo(GHOST_CELL_TYPE ghost_type){
//    return true;
// }

// ===================================================================
// PERIODIC_BOUNDARY_ELEMENT
// ===================================================================
// ===================================================================
INT PERIODIC_BOUNDARY_ELEMENT::GetOriginalVertIndex(INT j){
   // original_vert.
   // INT index;
   INT original_vert_0_index;
   ELEMENT * elem_ptr = GetGrid()->elem[original_elem];

   // 先找到 0 号顶点所对应的顶点在 original 单元上的编号.
   for(int i = 0; i<NVert; i++){
      if(elem_ptr->verts[i] == original_vert){
         // original_vert_0_index = i;
         original_vert_0_index = (i+1)%NVert;
         break;
      }
   }
   return (original_vert_0_index + j)%NVert;
}
