/*
  File:      grid.cc
  Purpose:

  Author: <kw.xu@hotmail.com>
  Last modified: "2016-04-22 16:34:31 kaiwen"
*/

#include "grid.h"
#include "def.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <cmath>


using namespace std;

// ===================================================================
// class GRID
// ===================================================================

// GRID::GRID() : NVert(0), NElem(0), vert(NULL), elem(NULL), bdry(NULL){}
// GRID::GRID() : n_vert(0), n_elem(0){}
GRID::~GRID(){
//    // delete [] elem;
//    // delete [] vert_coord;
   std::vector<ELEMENT*>::iterator it;
   for(it = elem.begin(); it != elem.end(); it++)
      delete (*it);
}
// GRID::~GRID(){};
bool GRID::ReadGridFile(const string & elemfile, const string & vertfile){
   string temp_str;
   INT temp_int;
   INT i,j;
   ifstream f_in;

   // 存储读取文件信息
   _elemfile=elemfile;
   _vertfile=vertfile;

   // ===================================================================
   f_in.open(elemfile.c_str());
   if( !f_in ){
      cerr << "GRID, ERROR: unable to open \""
           << elemfile << "\"\n";
      return false;
   }
   cout<<"GRID, Start reading element file (.e)...";
   // read element number.
   f_in>>n_elem;
   // create the such number of element objects.
   // elem = new ELEMENT [NElem];
   // elem.reserve(n_elem); 这里不应该用 reserve, 应该用 resize.
   // 注意 elem 里面存储的是指针.
   elem.resize(n_elem);

   for(i=0; i<n_elem; i++){
      // elem.push_back(ELEMENT(this));
      elem[i] = new ELEMENT(this); // 注意到 elem 是指针的 array.
      elem[i]->index = i;
      f_in >> temp_str;
      for(j=0; j<NEdge; j++)
         f_in >> elem[i]->verts[j];
      for(j=0; j<NEdge; j++)
         f_in >> elem[i]->nghbs[j]; // 不做任何关于边界方面的处理,
      // 边界相关的处理独立放在 BOUNDARY 里面.
      f_in >> temp_int;
      if(temp_int == BDMARK)
         elem[i]->isbd = true;   // false in ELEMENT initializition.
   }
   cout<<"Over!"<<endl;
   f_in.close();

   // ===================================================================
   f_in.open(vertfile.c_str());
   if( !f_in ){
      cerr << "GRID, ERROR: unable to open \""
           << vertfile << "\"\n";
      return false;
   }
   cout<<"GRID, Start reading vertex file (.n)...";
   // read vertex number
   f_in>>n_vert;
   // vert_coord = new VECT [NVert];
   // _vert_coord.reserve(n_vert);
   vert_coord.resize(n_vert);
   // 注意 vertex 的 mark 位在程序里面没有用到.
   for(i=0; i<n_vert; i++){
      f_in >> temp_str;
      f_in >> vert_coord[i].x;
      f_in >> vert_coord[i].y;
      f_in >> temp_str;
   }
   cout<<"Over!"<<endl;
   cout<<"NElement: "<<n_elem <<endl;
   cout<<"NVert: "<<n_vert<<endl;
   f_in.close();
   return true;
}

// ===================================================================
bool GRID::Precondition(){
   // 更新单元相关的物理量
   std::cout<<"GRID, Preconditioning...";
   std::vector<ELEMENT*>::iterator it;
   for (it = elem.begin(); it != elem.end(); ++it)
      (*it)->UpdateElemInfo();
   std::cout<<"Over!"<<endl;
   return true;
}

// ===================================================================
void GRID::ShowInfo(){
   std::cout<<endl
            <<"=========================="<< endl
            <<"GRID.ShowInfo:"<< endl
            <<"   elemfile -- "<< _elemfile << endl
            <<"   vertfile -- "<< _vertfile << endl
            <<"   n_vert   -- "<< n_vert << endl
            <<"   n_elem   -- "<< n_elem << endl
            <<"=========================="<< endl;
   // std::cout<<"Vertex:"<<endl;
   // for(int i = 0; i< n_vert; i++){
   //    std::cout<<vert_coord[i];
   // }
   // std::vector<ELEMENT*>::iterator it;
   // for (it = elem.begin(); it != elem.end(); ++it)
   //    std::cout<<(*it)->barycenter;
   // return;
   // std::cout<<"neighbours:"<<endl;
   // std::vector<ELEMENT*>::iterator it;
   // for (it = elem.begin(); it != elem.end(); ++it)
   //    std::cout<<setw(10)<<(*it)->nghbs[0]
   //             <<setw(10)<<(*it)->nghbs[1]
   //             <<setw(10)<<(*it)->nghbs[2]
   //             <<endl;
   // std::cout<<(*it)->barycenter;
   return;

}

// ===================================================================
// class ELEMENT
// ===================================================================

bool ELEMENT::UpdateElemInfo(){
   // grid = _grid;
   VECT vc[NEdge];
   INT ii, jj, kk;
   for(ii=0,jj=1;ii<NEdge;++ii,jj=(jj+1)%NEdge)
      vc[(ii+2)%NEdge] = grid->vert_coord[verts[jj]] - grid->vert_coord[verts[ii]];

   // area
   // VECT vc1, vc2;
   // vc1.x = _grid->vert[verts[jj]].x - _grid->vert[verts[ii]].x;
   // vc1.y = _grid->vert[verts[jj]].y - _grid->vert[verts[ii]].y;
   // vc[0] = _grid->vert_coord[verts[jj]] - _grid->vert_coord[verts[ii]];
   // vc[1] = _grid->vert_coord[verts[kk]] - _grid->vert_coord[verts[jj]];
   area = OuterProduct(vc[0], vc[1]);
   area = 0.5 * std::fabs(area);

   // barycenter
   barycenter.x=0.0;
   barycenter.y=0.0;
   for(ii=0;ii<NVert;ii++)
      barycenter = barycenter + grid->vert_coord[verts[ii]];
   barycenter = (1.0/3.0) * barycenter;

   // orient
   // 因为 vc1, vc2, vc3 按逆时针, 所以外侧指的是每个向量的右侧.
   // 小心标号的对应.
   for(ii=0;ii<NEdge;ii++){
      // 顺时针旋转 \pi/2.
      orient[ii].x = vc[ii].y;
      orient[ii].y = - vc[ii].x;
      // 单位化.
      orient[ii] = orient[ii].Normalize();
   }

   FLOAT dx, dy;
   // edge_len
   for(ii = 0; ii < NEdge; ii++){
      jj = (ii + 1)%NEdge;
      kk = (ii + 2)%NEdge;
      dx = grid->vert_coord[verts[jj]].x - grid->vert_coord[verts[kk]].x;
      dy = grid->vert_coord[verts[jj]].y - grid->vert_coord[verts[kk]].y;
      edge_len[ii] = std::sqrt(dx * dx + dy * dy);
   }

   // 坐标变换相关的几何量.
   for(ii=0,jj=1,kk=2; ii<NEdge; ii++, jj=(jj+1)%NEdge, kk=(kk+1)%NEdge){
      omega[ii] =
         grid->vert_coord[verts[jj]].x * grid->vert_coord[verts[kk]].y -
         grid->vert_coord[verts[kk]].x * grid->vert_coord[verts[jj]].y;
      eta[ii] = grid->vert_coord[verts[jj]].y - grid->vert_coord[verts[kk]].y;
      xi[ii]  = grid->vert_coord[verts[jj]].x - grid->vert_coord[verts[kk]].x;
   }
   // Jacobian
   int e;
   e = 0;
   for(int l = 0; l < NVert; l++)
      elem_jacobian[e][l] = eta[l] / (2.0 * area);
   e = 1;
   for(int l = 0; l < NVert; l++)
      elem_jacobian[e][l] = - xi[l] / (2.0 * area);

   // nghbr_edge_index
   // CANCEAL
   return true;
}

INT ELEMENT::GetNghbrEdgeIndex(INT i){
   INT edge_tail = verts[(i+1)%NEdge];
   INT index;
   for(int j = 0; j<NEdge; j++){
      if(GetNeighborElem(i)->verts[j] == edge_tail){
         index = (j+1)%NEdge;
         break;
      }
   }
   return index;
}

void ELEMENT::XY2Lambda(BARY_COOD & bary_cood,const VECT &vc){
   for (int i=0; i<NVert; i++)
      bary_cood[i] = (omega[i] + eta[i] * vc.x - xi[i] * vc.y) / (2 * area);
   return;
}
