/*
  File:      recon.cc
  Purpose:

  Author: <kw.xu@hotmail.com>
  Last modified: "2016-05-15 10:35:06 kaiwen"
*/

#include "recon.h"
#include <cmath>
#include <limits>

// ===================================================================

void MUSCL_RECON::Reconstruction(ELEMENT * elem,
                                 FLOAT * recon_value){
   INT i;
   INT elem_index = elem->index;
   VECT elem_centre = elem->barycenter;
   VECT edge_mid;
   FLOAT cell_average = _solution[elem_index];

   // pre, Pre 应该在 elem 循环的外层做. (是对整个流场的操作)
   // Pre_Treatment();

   // predict
   _predicter.init_grad(elem, _init_grad);

   // limite
   _limiter.apply_limiter(elem, _init_grad, _limit_grad);

   // 计算单元边界中点值.
   for(i = 0 ; i < NEdge; i ++)
   {
      edge_mid = elem->GetEdgeMidPosition(i);
      recon_value[i] = InerProduct(_limit_grad, edge_mid - elem_centre) + cell_average;
   }

   return;
}

void MUSCL_RECON::Pre_Treatment()
{
   _limiter.pretreatment();
}


// ===================================================================

void GRAD_PREDICT_LS::init_grad(ELEMENT * elem, VECT & grad)
{
   FLOAT l_11, l_22, l_12, lf_1, lf_2;
   FLOAT multiplier;
   VECT elem_centre = elem->barycenter;
   VECT nbr_elem_centre;
   ELEMENT * nbr_elem;
   FLOAT L_1[NEdge],L_2[NEdge],f[NEdge];
   FLOAT elem_average = (*_solution)[elem->index];

   for(int i = 0 ; i < NEdge; i ++)
   {
      nbr_elem = (elem->GetGrid())->elem[elem->nghbs[i]];
      nbr_elem_centre = nbr_elem->barycenter;
      L_1[i] = nbr_elem_centre.x - elem_centre.x;
      L_2[i] = nbr_elem_centre.y - elem_centre.y;
      // Be careful about the boundary treatment.
      f[i]   = (*_solution)[elem->nghbs[i]] - elem_average;
   }

   l_11 = l_22 = l_12 = 0;
   lf_1 = lf_2 = 0;
   for(int i = 0 ; i < NEdge; i ++)
   {
      l_11 += L_1[i] * L_1[i];
      l_22 += L_2[i] * L_2[i];
      l_12 += L_1[i] * L_2[i];
      lf_1 += L_1[i] * f[i];
      lf_2 += L_2[i] * f[i];
   }

   multiplier = 1 / (l_11 * l_22 - l_12 * l_12);

   grad.x = multiplier * (l_22 * lf_1 - l_12 * lf_2);
   grad.y = multiplier * (l_11 * lf_2 - l_12 * lf_1);

   return;
}

// ===================================================================

void GRAD_LIMITER_MLP::pretreatment()
{
   INT i, k, l;
   INT vert_index;
   FLOAT cell_average;
   ELEMENT * current_elem;
   GRID * grid = _solution->GetGrid();
   BOUNDARY * boundary = _solution->GetBoundary();

   for(i = 0; i<grid->n_vert; i++){
      _max_min_vt[i][0] = std::numeric_limits<float>::min();
      _max_min_vt[i][1] = std::numeric_limits<float>::max();
   }

   //,---------------------------
   //| 计算顶点处允许的最大最小值
   //`---------------------------

   // 对区域内单元的循环.
   // 更精确的确定断点值的范围 <2015-07-03 Fri 16:55>
   // <2015-07-03 Fri 17:02> 修改回原来的方案.
   for(l = 0; l<grid->n_elem; l++){
      // u0 = _dofs[l]->GetAverage();
      cell_average = (*_solution)[l];
      current_elem = grid->elem[l];
      for(k = 0; k<NVert; k++){
         // vert_index = _dofs[l]->GetVertexIndex(k);
         vert_index = current_elem->verts[k];
         _max_min_vt[vert_index][0] = std::max(_max_min_vt[vert_index][0],
                                               cell_average);
         _max_min_vt[vert_index][1] = std::min(_max_min_vt[vert_index][1],
                                               cell_average);
      }
   }

   INT boundary_vert_index, original_vert_index;
   BOUNDARY_ELEMENT * boundary_elem;
   ELEMENT * original_elem;
   // 对边界上的点进行处理. 将对应的两个节点上的最值进行整合.
   for(l = 0 ; l < boundary->n_boundary ; l ++)
   {
      // 定位对边上的两个对应点.
      boundary_elem = boundary->boundary_elem[l];
      boundary_vert_index = boundary_elem->verts[0];
      original_elem = grid->elem[boundary_elem->original_elem];
      original_vert_index = original_elem->verts[boundary_elem->GetOriginalVertIndex(0)];

      // 最大值取较大.
      // _max_min_vt[boundary_vert_index][0] > _max_min_vt[original_vert_index][0] ?
      //    _max_min_vt[original_vert_index][0] = _max_min_vt[boundary_vert_index][0] :
      //    _max_min_vt[boundary_vert_index][0] = _max_min_vt[original_vert_index][0];
      if(_max_min_vt[boundary_vert_index][0] > _max_min_vt[original_vert_index][0])
      {
         _max_min_vt[original_vert_index][0] = _max_min_vt[boundary_vert_index][0];
      }
      else
      {
         _max_min_vt[boundary_vert_index][0] = _max_min_vt[original_vert_index][0];
      }
      // 最小值取较小.
      // _max_min_vt[boundary_vert_index][1] > _max_min_vt[original_vert_index][1] ?
      //    _max_min_vt[boundary_vert_index][1] = _max_min_vt[original_vert_index][1] :
      //    _max_min_vt[original_vert_index][1] = _max_min_vt[boundary_vert_index][1];
      if(_max_min_vt[boundary_vert_index][1] > _max_min_vt[original_vert_index][1])
      {
         _max_min_vt[boundary_vert_index][1] = _max_min_vt[original_vert_index][1];
      }
      else
      {
         _max_min_vt[original_vert_index][1] = _max_min_vt[boundary_vert_index][1];
      }
   }
   // 四个角点上的处理要非常小心.
   FLOAT conner_max = std::numeric_limits<float>::min();
   FLOAT conner_min = std::numeric_limits<float>::max();
   for(k = 0 ; k < 4 ; k ++)
   {
      conner_max = conner_max > _max_min_vt[boundary->GetConnerIndex(k)][0] ?
         conner_max : _max_min_vt[boundary->GetConnerIndex(k)][0];
      conner_min = conner_min < _max_min_vt[boundary->GetConnerIndex(k)][1] ?
         conner_min : _max_min_vt[boundary->GetConnerIndex(k)][1];
   }

   for(k = 0 ; k < 4 ; k ++)
   {
      _max_min_vt[boundary->GetConnerIndex(k)][0] = conner_max;
      _max_min_vt[boundary->GetConnerIndex(k)][1] = conner_min;
   }

   return;
}



void GRAD_LIMITER_MLP::apply_limiter(ELEMENT * elem,
                                     const VECT & init_grad,
                                     VECT & limited_grad)
{
   INT k;
   GRID * grid = _solution->GetGrid();
   FLOAT cell_average = (*_solution)[elem->index];
   FLOAT alpha = 1.0;
   FLOAT duv;
   VECT cell_centre = elem->barycenter;
   VECT vertex_pos;
   INT vert_index;

   for(k = 0 ; k < NVert ; k ++)
   {
      vert_index = elem->verts[k];
      vertex_pos = grid->vert_coord[vert_index];
      duv = InerProduct(init_grad, vertex_pos - cell_centre);
      if(std::fabs(duv)>TOL) // 防止出现 duv \aprox 0 的情形. TOL defined in def.h
      {
         alpha = std::min(alpha,
                          std::max(0.0,
                                   std::max((_max_min_vt[vert_index][0]-
                                             cell_average)/duv,
                                            (_max_min_vt[vert_index][1]-
                                             cell_average)/duv)));
      }
      else
      {
         alpha = 1.0;
      }
   }

   limited_grad = alpha * init_grad;

   return;
}


// ===================================================================
// GRAD_LIMITER_MLP_MIDPOINT

void GRAD_LIMITER_MLP_MIDPOINT::pretreatment()
{
   // 方便起见, 直接 copy 的上面 original MLP 的部分.
   INT i, k, l;
   INT vert_index;
   FLOAT cell_average;
   ELEMENT * current_elem;
   GRID * grid = _solution->GetGrid();
   BOUNDARY * boundary = _solution->GetBoundary();

   for(i = 0; i<grid->n_vert; i++){
      _max_min_vt[i][0] = std::numeric_limits<float>::min();
      _max_min_vt[i][1] = std::numeric_limits<float>::max();
   }

   //,---------------------------
   //| 计算顶点处允许的最大最小值
   //`---------------------------

   // 对区域内单元的循环.
   // 更精确的确定断点值的范围 <2015-07-03 Fri 16:55>
   // <2015-07-03 Fri 17:02> 修改回原来的方案.
   for(l = 0; l<grid->n_elem; l++){
      // u0 = _dofs[l]->GetAverage();
      cell_average = (*_solution)[l];
      current_elem = grid->elem[l];
      for(k = 0; k<NVert; k++){
         // vert_index = _dofs[l]->GetVertexIndex(k);
         vert_index = current_elem->verts[k];
         _max_min_vt[vert_index][0] = std::max(_max_min_vt[vert_index][0],
                                               cell_average);
         _max_min_vt[vert_index][1] = std::min(_max_min_vt[vert_index][1],
                                               cell_average);
      }
   }

   INT boundary_vert_index, original_vert_index;
   BOUNDARY_ELEMENT * boundary_elem;
   ELEMENT * original_elem;
   // 对边界上的点进行处理. 将对应的两个节点上的最值进行整合.
   for(l = 0 ; l < boundary->n_boundary ; l ++)
   {
      // 定位对边上的两个对应点.
      boundary_elem = boundary->boundary_elem[l];
      boundary_vert_index = boundary_elem->verts[0];
      original_elem = grid->elem[boundary_elem->original_elem];
      original_vert_index = original_elem->verts[boundary_elem->GetOriginalVertIndex(0)];

      // 最大值取较大.
      // _max_min_vt[boundary_vert_index][0] > _max_min_vt[original_vert_index][0] ?
      //    _max_min_vt[original_vert_index][0] = _max_min_vt[boundary_vert_index][0] :
      //    _max_min_vt[boundary_vert_index][0] = _max_min_vt[original_vert_index][0];
      if(_max_min_vt[boundary_vert_index][0] > _max_min_vt[original_vert_index][0])
      {
         _max_min_vt[original_vert_index][0] = _max_min_vt[boundary_vert_index][0];
      }
      else
      {
         _max_min_vt[boundary_vert_index][0] = _max_min_vt[original_vert_index][0];
      }
      // 最小值取较小.
      // _max_min_vt[boundary_vert_index][1] > _max_min_vt[original_vert_index][1] ?
      //    _max_min_vt[boundary_vert_index][1] = _max_min_vt[original_vert_index][1] :
      //    _max_min_vt[original_vert_index][1] = _max_min_vt[boundary_vert_index][1];
      if(_max_min_vt[boundary_vert_index][1] > _max_min_vt[original_vert_index][1])
      {
         _max_min_vt[boundary_vert_index][1] = _max_min_vt[original_vert_index][1];
      }
      else
      {
         _max_min_vt[original_vert_index][1] = _max_min_vt[boundary_vert_index][1];
      }
   }
   // 四个角点上的处理要非常小心.
   FLOAT conner_max = std::numeric_limits<float>::min();
   FLOAT conner_min = std::numeric_limits<float>::max();
   for(k = 0 ; k < 4 ; k ++)
   {
      conner_max = conner_max > _max_min_vt[boundary->GetConnerIndex(k)][0] ?
         conner_max : _max_min_vt[boundary->GetConnerIndex(k)][0];
      conner_min = conner_min < _max_min_vt[boundary->GetConnerIndex(k)][1] ?
         conner_min : _max_min_vt[boundary->GetConnerIndex(k)][1];
   }

   for(k = 0 ; k < 4 ; k ++)
   {
      _max_min_vt[boundary->GetConnerIndex(k)][0] = conner_max;
      _max_min_vt[boundary->GetConnerIndex(k)][1] = conner_min;
   }

   return;
}

void GRAD_LIMITER_MLP_MIDPOINT::apply_limiter(ELEMENT * elem, const VECT & init_grad,
                                              VECT & limited_grad)
{
   INT k;
   GRID * grid = _solution->GetGrid();
   FLOAT cell_average = (*_solution)[elem->index];
   FLOAT alpha = 1.0;
   FLOAT duv;
   VECT cell_centre = elem->barycenter;
   // VECT vertex_pos;
   VECT edge_mid_pos;
   INT vert_index;
   FLOAT local_max, local_min;
   INT ii,jj,kk;

   for(k = 0 ; k < NEdge ; k ++)
   {
      ii = k;
      jj = (ii+1)%NEdge;
      kk = (ii+2)%NEdge;

      local_max = std::max(
         _max_min_vt[elem->verts[jj]][0],
         _max_min_vt[elem->verts[kk]][0]);
      local_min = std::min(
         _max_min_vt[elem->verts[jj]][1],
         _max_min_vt[elem->verts[kk]][1]);

      // local_max = 0.5 *
      //    (_max_min_vt[elem->verts[jj]][0],
      //     _max_min_vt[elem->verts[kk]][0]);
      // local_max = 0.5 *
      //    (_max_min_vt[elem->verts[jj]][1],
      //     _max_min_vt[elem->verts[kk]][1]);

      // vert_index = elem->verts[k];
      // vertex_pos = grid->vert_coord[vert_index];
      edge_mid_pos = elem->GetEdgeMidPosition(k);
      duv = InerProduct(init_grad, edge_mid_pos - cell_centre);
      if(std::fabs(duv)>TOL) // 防止出现 duv \aprox 0 的情形. TOL defined in def.h
      {
         alpha = std::min(alpha,
                          std::max(0.0,
                                   std::max((local_max - cell_average)/duv,
                                            (local_min - cell_average)/duv)));
      }
      else
      {
         alpha = 1.0;
      }
   }

   limited_grad = alpha * init_grad;

   return;
}

GRAD_LIMITER_MLP_LIST::GRAD_LIMITER_MLP_LIST(SOLUTION * solution) :
   GRAD_LIMITER(solution){
   GRID * grid = _solution->GetGrid();
   BOUNDARY * boundary = _solution->GetBoundary();
   INT j,k,l;
   INT vert_index;
   ELEMENT * elem;

   vertex_neighbour_elem.resize(grid->n_vert);

   for(j = 0 ; j < grid->n_elem ; j ++){
      elem = grid->elem[j];
      for(k = 0 ; k < NVert ; k ++){
         vert_index = elem->verts[k];
         vertex_neighbour_elem[vert_index].push_front(elem->index);
      }
   }

   // 对边界的特殊处理.
   // 1. 对边上的对应点的 list 整合到一起.
   INT boundary_vert_index, original_vert_index;
   BOUNDARY_ELEMENT * boundary_elem;
   ELEMENT * original_elem;
   for(l = 0 ; l < boundary->n_boundary ; l ++)
   {
      // 定位对边上的两个对应点.
      boundary_elem = boundary->boundary_elem[l];
      boundary_vert_index = boundary_elem->verts[0];
      original_elem = grid->elem[boundary_elem->original_elem];
      original_vert_index =
         original_elem->verts[boundary_elem->GetOriginalVertIndex(0)];

      vertex_neighbour_elem[boundary_vert_index].insert_after(
         vertex_neighbour_elem[boundary_vert_index].begin(),
         vertex_neighbour_elem[original_vert_index].begin(),
         vertex_neighbour_elem[original_vert_index].end()
         );
      vertex_neighbour_elem[original_vert_index] =
         vertex_neighbour_elem[boundary_vert_index];
   }

   // 2. 对四个角点进行特殊处理.
   for(k = 1 ; k < 4 ; k ++){
      vertex_neighbour_elem[boundary_vert_index].insert_after(
         vertex_neighbour_elem[boundary->GetConnerIndex(0)].begin(),
         vertex_neighbour_elem[boundary->GetConnerIndex(k)].begin(),
         vertex_neighbour_elem[boundary->GetConnerIndex(k)].end()
         );
   }
   for(k = 1 ; k < 4 ; k ++){
      vertex_neighbour_elem[boundary->GetConnerIndex(k)] =
         vertex_neighbour_elem[boundary->GetConnerIndex(0)];
   }
}


void GRAD_LIMITER_MLP_LIST::apply_limiter(ELEMENT * elem,
                                          const VECT & init_grad,
                                          VECT & limited_grad)
{
   INT k;
   GRID * grid = _solution->GetGrid();
   FLOAT cell_average = (*_solution)[elem->index];
   FLOAT alpha = 1.0;
   FLOAT duv;
   VECT cell_centre = elem->barycenter;
   VECT vertex_pos;
   INT vert_index;

   FLOAT vertex_max, vertex_min;

   for(k = 0 ; k < NVert ; k ++)
   {
      vert_index = elem->verts[k];
      vertex_pos = grid->vert_coord[vert_index];
      duv = InerProduct(init_grad, vertex_pos - cell_centre);

      vertex_max = std::numeric_limits<float>::min();
      vertex_min = std::numeric_limits<float>::max();
      for(auto it = vertex_neighbour_elem[vert_index].begin() ;
          it != vertex_neighbour_elem[vert_index].end() ;
          it ++
         ) {
         vertex_max = std::max(vertex_max, (*_solution)[*it]);
         vertex_min = std::min(vertex_min, (*_solution)[*it]);
      }

      if(std::fabs(duv)>TOL) // 防止出现 duv \aprox 0 的情形. TOL defined in def.h
      {
         alpha = std::min(alpha,
                          std::max(0.0,
                                   std::max((vertex_max - cell_average)/duv,
                                            (vertex_min - cell_average)/duv)));
      }
      else
      {
         alpha = 1.0;
      }
   }

   limited_grad = alpha * init_grad;

   return;
}


// ===================================================================
// Barth-Jespersen

void GRAD_LIMITER_BARTH_JESPERSEN::apply_limiter(ELEMENT * elem, const VECT & init_grad,
                                                 VECT & limited_grad)
{
   INT k;
   // GRID * grid = _solution->GetGrid();
   FLOAT local_max, local_min;
   FLOAT nbr_value[NEdge];
   VECT cell_centre = elem->barycenter;
   FLOAT cell_average = (*_solution)[elem->index];
   FLOAT alpha = 1.0;
   FLOAT duv;
   VECT edge_mid_pos;

   local_max = cell_average;
   local_min = cell_average;
   for(k = 0 ; k < NEdge ; k ++)
   {
      nbr_value[k] = (*_solution)[elem->nghbs[k]];
      local_max = std::max(local_max, nbr_value[k]);
      local_min = std::min(local_min, nbr_value[k]);
   }

   for(k = 0 ; k < NVert ; k ++)
   {
      edge_mid_pos = elem->GetEdgeMidPosition(k);
      duv = InerProduct(init_grad, edge_mid_pos - cell_centre);
      if(std::fabs(duv)>TOL) // 防止出现 duv \aprox 0 的情形. TOL defined in def.h
      {
         alpha = std::min(alpha,
                          std::max(0.0,
                                   std::max((local_max - cell_average)/duv,
                                            (local_min - cell_average)/duv)));
      }
      else
      {
         alpha = 1.0;
      }
   }

   limited_grad = alpha * init_grad;

   return;
}

// ===================================================================
// LCD

void GRAD_LIMITER_LCD::apply_limiter(ELEMENT * elem, const VECT & init_grad,
                                     VECT & limited_grad)
{
   INT k;
   FLOAT cell_average = (*_solution)[elem->index];
   FLOAT local_max, local_min;
   FLOAT alpha = 1.0;
   FLOAT duv;
   // FLOAT nbr_value[NEdge];
   VECT cell_centre = elem->barycenter;
   VECT edge_mid_pos;

   for(k = 0 ; k < NVert ; k ++)
   {
      local_max = std::max(cell_average, (*_solution)[elem->nghbs[k]]);
      local_min = std::min(cell_average, (*_solution)[elem->nghbs[k]]);
      edge_mid_pos = elem->GetEdgeMidPosition(k);
      duv = InerProduct(init_grad, edge_mid_pos - cell_centre);
      if(std::fabs(duv)>TOL) // 防止出现 duv \aprox 0 的情形. TOL defined in def.h
      {
         alpha = std::min(alpha,
                          std::max(0.0,
                                   std::max((local_max - cell_average)/duv,
                                            (local_min - cell_average)/duv)));
      }
      else
      {
         alpha = 1.0;
      }
   }

   limited_grad = alpha * init_grad;

   return;
}

// ===================================================================
// Vertex version of barth.

void GRAD_LIMITER_VERTEX_BARTH::apply_limiter(ELEMENT * elem, const VECT & init_grad,
                                              VECT & limited_grad)
{
   INT k;
   GRID * grid = _solution->GetGrid();
   FLOAT nbr_value[NEdge];
   FLOAT local_max, local_min;
   FLOAT cell_average = (*_solution)[elem->index];
   FLOAT alpha = 1.0;
   FLOAT duv;
   VECT cell_centre = elem->barycenter;
   VECT vertex_pos;
   INT vert_index;

   local_max = cell_average;
   local_min = cell_average;
   for(k = 0 ; k < NEdge ; k ++)
   {
      nbr_value[k] = (*_solution)[elem->nghbs[k]];
      local_max = std::max(local_max, nbr_value[k]);
      local_min = std::min(local_min, nbr_value[k]);
   }

   for(k = 0 ; k < NVert ; k ++)
   {
      vert_index = elem->verts[k];
      vertex_pos = grid->vert_coord[vert_index];
      duv = InerProduct(init_grad, vertex_pos - cell_centre);
      if(std::fabs(duv)>TOL) // 防止出现 duv \aprox 0 的情形. TOL defined in def.h
      {
         alpha = std::min(alpha,
                          std::max(0.0,
                                   std::max((local_max - cell_average)/duv,
                                            (local_min - cell_average)/duv)));
      }
      else
      {
         alpha = 1.0;
      }
   }

   limited_grad = alpha * init_grad;

   return;
}

// ===================================================================
// Cockburn & Shu

void GRAD_LIMITER_SHU_COCKBURN::apply_limiter(ELEMENT * elem, const VECT & init_grad,
                                              VECT & limited_grad)
{
   FLOAT nu = 1.0;
   FLOAT L_1[2], L_2[2], f[2], r_i[2], alpha[2];
   bool is_stencil;
   INT i, j, l, k, ii, jj;
   INT nb_ii,nb_jj;
   VECT nb_cen[2];
   FLOAT nb_value[2], cell_value;
   VECT mid_pos;
   VECT cell_cen;
   FLOAT du[NEdge], dup[NEdge];
   FLOAT pos, neg;
   FLOAT theta_pos, theta_neg;
   FLOAT res;
   // INT vert_index[2];
   INT cell_index = elem->index;
   GRID * grid = _solution->GetGrid();

   // for(l = 0 ; l < _grid->n_elem ; l ++){
   //,-------
   //| Step.1 calculate \Delta_i
   //`-------
   // cell_value = _dofs[l]->GetAverage();
   cell_value = (*_solution)[cell_index];
   // cell_cen = _dofs[l]->GetCellCenter();
   cell_cen = elem->barycenter;
   for(k = 0 ; k < NEdge ; k++){
      // vert_index[0] = _dofs[l]->GetVertexIndex((k+1)%NEdge);
      // vert_index[1] = _dofs[l]->GetVertexIndex((k+2)%NEdge);
      // mid_pos = 0.5 * (_dofs[l]->GetVertexIndex((k+1)%NEdge) +
      //                  _dofs[l]->GetVertexIndex((k+2)%NEdge));
      // mid_pos = 0.5 * (_grid->vert_coord[vert_index[0]] +
      //                  _grid->vert_coord[vert_index[1]]);
      mid_pos = elem->GetEdgeMidPosition(k);
      r_i[0] = mid_pos.x - cell_cen.x;
      r_i[1] = mid_pos.y - cell_cen.y;
      // 到边中点处的增量.
      // du[k] = (*_dofs[l])[k] - _dofs[l]->GetAverage();
      du[k] = InerProduct(init_grad, mid_pos - cell_cen);
      // Step 1.1 find stencil
      is_stencil = false;
      for(ii = 0; ii < NEdge ; ii ++){
         jj = (ii + 1) % NEdge;
         // nb_ii = _dofs[l]->GetNeighborElemIndex(ii);
         // nb_jj = _dofs[l]->GetNeighborElemIndex(jj);
         nb_ii = elem->nghbs[ii];
         nb_jj = elem->nghbs[jj];
         // nb_cen[0] = _dofs[nb_ii]->GetCellCenter();
         // nb_cen[1] = _dofs[nb_jj]->GetCellCenter();
         nb_cen[0] = grid->elem[nb_ii]->barycenter;
         nb_cen[1] = grid->elem[nb_jj]->barycenter;
         // nb_value[0] = _dofs[nb_ii]->GetAverage();
         // nb_value[1] = _dofs[nb_jj]->GetAverage();
         nb_value[0] = (*_solution)[nb_ii];
         nb_value[1] = (*_solution)[nb_jj];
         L_1[0] = nb_cen[0].x - cell_cen.x;
         L_1[1] = nb_cen[1].x - cell_cen.x;
         L_2[0] = nb_cen[0].y - cell_cen.y;
         L_2[1] = nb_cen[1].y - cell_cen.y;
         f[0] = nb_value[0] - cell_value;
         f[1] = nb_value[1] - cell_value;

         alpha[0] =
            (L_2[1] * r_i[0] - L_1[1] * r_i[1]) /
            (L_1[0] * L_2[1] - L_2[0] * L_1[1]);
         alpha[1] =
            (- L_2[0] * r_i[0] + L_1[0] * r_i[1]) /
            (L_1[0] * L_2[1] - L_2[0] * L_1[1]);

         // if(alpha[0] > 0 && alpha[1] > 0){
         if(alpha[0] > -TOL && alpha[1] > -TOL){
            is_stencil = true;
            break;
         }
      }
      if(!is_stencil){
         // 没有找到合适的模板, 计算失败.
         std::cout << "_ShuCockburnLimiter Error, can't find the appropriate stencil"
                   << std::endl;
         exit(1);
      }
      // Step 1.2 cal dup (通过之前的模板预测的 m_i 处的值)
      dup[k] = 0;
      for(j = 0 ; j < 2 ; j ++){
         dup[k] += alpha[j] * f[j];
      }

      // Step 1.3 cal \Delta_i using TVB type scaler limiter.
      if(du[k] * dup[k] < 0){
         du[k] = 0;
      } else if(du[k] > 0){
         du[k] = std::min(du[k], nu * dup[k]);
      } else {
         // < 0
         du[k] = std::max(du[k], nu * dup[k]);
      }
   } // each edge

   //,------------------
   //| Step.2 modify the increamental in satisfaction of conservation.
   //`------------------

   // 算法分三种情况
   // 1. 三个增量的和为 0
   //    随意选两个点计算 grad 即可.
   // 2. 需要被 modify 的只有一个点,
   //    即 Shu 1998 pp.211 最后一个式子只对三个点中的一个起作用.
   //    a. 首先计算另外两个点与当前单元中点得到梯度.
   //    b. 将这个梯度值根据剩下的一个点进行 modify.
   // 3. 需要被 modify 的有两个点.
   //    pp.211 最后一个式子对两个点起作用.
   //    用这两个点以及单元重心可以即可得到斜率.
   // 这样的限制按理来讲应该比只计算单一的 alpha 值耗散要低一些.
   // (考虑针对 MLP 进行一下尝试)

   res = 0;
   for(k = 0; k < NEdge ; k ++){
      res += du[k];
   }
   // if(res == 0){
   //    // 选择 0, 1 两个边中点.
   //    ii = 0;
   //    jj = 1;
   //    is_single_modified = false;
   // } else {
   if(res != 0){
      // step.1 计算 modified increment.
      pos = 0;
      neg = 0;
      for(k = 0; k < NEdge ; k ++){
         pos += std::max(0.0, du[k]);
         neg += std::max(0.0, -du[k]);
      }
      theta_pos = std::min(1.0, neg/pos);
      theta_neg = std::min(1.0, pos/neg);
      for(k = 0 ; k < NEdge ; k ++){
         du[k] =
            theta_pos * std::max(0.0, du[k]) -
            theta_neg * std::max(0.0, -du[k]);
         // (*_dofs[l])[k] = cell_value + du[k];
         // (*_dofs[l])[k] = cell_value + du[k];
      }
   } // end of if

   ii = 0;
   jj = 1;
   f[0] = du[ii];
   f[1] = du[jj];
   mid_pos = elem->GetEdgeMidPosition(ii);
   L_1[0] = mid_pos.x - cell_cen.x;
   L_2[0] = mid_pos.y - cell_cen.y;
   mid_pos = elem->GetEdgeMidPosition(jj);
   L_1[1] = mid_pos.x - cell_cen.x;
   L_2[1] = mid_pos.y - cell_cen.y;

   limited_grad.x =
      (L_2[1] * f[0] - L_2[0] * f[1]) /
      (L_1[0] * L_2[1] - L_2[0] * L_1[1]);
   limited_grad.y =
      (- L_1[1] * f[0] + L_1[0] * f[1]) /
      (L_1[0] * L_2[1] - L_2[0] * L_1[1]);

   return;
}


// ===================================================================
// MLP_LIMITE_ON_VALUE


void GRAD_LIMITER_MLP_LIMITE_ON_VALUE::pretreatment()
{
   // 方便起见, 直接 copy 的上面 original MLP 的部分.
   INT i, k, l;
   INT vert_index;
   FLOAT cell_average;
   ELEMENT * current_elem;
   GRID * grid = _solution->GetGrid();
   BOUNDARY * boundary = _solution->GetBoundary();

   for(i = 0; i<grid->n_vert; i++){
      _max_min_vt[i][0] = std::numeric_limits<float>::min();
      _max_min_vt[i][1] = std::numeric_limits<float>::max();
   }

   //,---------------------------
   //| 计算顶点处允许的最大最小值
   //`---------------------------

   // 对区域内单元的循环.
   // 更精确的确定断点值的范围 <2015-07-03 Fri 16:55>
   // <2015-07-03 Fri 17:02> 修改回原来的方案.
   for(l = 0; l<grid->n_elem; l++){
      // u0 = _dofs[l]->GetAverage();
      cell_average = (*_solution)[l];
      current_elem = grid->elem[l];
      for(k = 0; k<NVert; k++){
         // vert_index = _dofs[l]->GetVertexIndex(k);
         vert_index = current_elem->verts[k];
         _max_min_vt[vert_index][0] = std::max(_max_min_vt[vert_index][0],
                                               cell_average);
         _max_min_vt[vert_index][1] = std::min(_max_min_vt[vert_index][1],
                                               cell_average);
      }
   }

   INT boundary_vert_index, original_vert_index;
   BOUNDARY_ELEMENT * boundary_elem;
   ELEMENT * original_elem;
   // 对边界上的点进行处理. 将对应的两个节点上的最值进行整合.
   for(l = 0 ; l < boundary->n_boundary ; l ++)
   {
      // 定位对边上的两个对应点.
      boundary_elem = boundary->boundary_elem[l];
      boundary_vert_index = boundary_elem->verts[0];
      original_elem = grid->elem[boundary_elem->original_elem];
      original_vert_index = original_elem->verts[boundary_elem->GetOriginalVertIndex(0)];

      // 最大值取较大.
      // _max_min_vt[boundary_vert_index][0] > _max_min_vt[original_vert_index][0] ?
      //    _max_min_vt[original_vert_index][0] = _max_min_vt[boundary_vert_index][0] :
      //    _max_min_vt[boundary_vert_index][0] = _max_min_vt[original_vert_index][0];
      if(_max_min_vt[boundary_vert_index][0] > _max_min_vt[original_vert_index][0])
      {
         _max_min_vt[original_vert_index][0] = _max_min_vt[boundary_vert_index][0];
      }
      else
      {
         _max_min_vt[boundary_vert_index][0] = _max_min_vt[original_vert_index][0];
      }
      // 最小值取较小.
      // _max_min_vt[boundary_vert_index][1] > _max_min_vt[original_vert_index][1] ?
      //    _max_min_vt[boundary_vert_index][1] = _max_min_vt[original_vert_index][1] :
      //    _max_min_vt[original_vert_index][1] = _max_min_vt[boundary_vert_index][1];
      if(_max_min_vt[boundary_vert_index][1] > _max_min_vt[original_vert_index][1])
      {
         _max_min_vt[boundary_vert_index][1] = _max_min_vt[original_vert_index][1];
      }
      else
      {
         _max_min_vt[original_vert_index][1] = _max_min_vt[boundary_vert_index][1];
      }
   }
   // 四个角点上的处理要非常小心.
   FLOAT conner_max = std::numeric_limits<float>::min();
   FLOAT conner_min = std::numeric_limits<float>::max();
   for(k = 0 ; k < 4 ; k ++)
   {
      conner_max = conner_max > _max_min_vt[boundary->GetConnerIndex(k)][0] ?
         conner_max : _max_min_vt[boundary->GetConnerIndex(k)][0];
      conner_min = conner_min < _max_min_vt[boundary->GetConnerIndex(k)][1] ?
         conner_min : _max_min_vt[boundary->GetConnerIndex(k)][1];
   }

   for(k = 0 ; k < 4 ; k ++)
   {
      _max_min_vt[boundary->GetConnerIndex(k)][0] = conner_max;
      _max_min_vt[boundary->GetConnerIndex(k)][1] = conner_min;
   }

   return;
}


void GRAD_LIMITER_MLP_LIMITE_ON_VALUE::apply_limiter(ELEMENT * elem,
                                                     const VECT & init_grad,
                                                     VECT & limited_grad)
{
   INT k;
   GRID * grid = _solution->GetGrid();
   FLOAT cell_average = (*_solution)[elem->index];
   FLOAT alpha = 1.0;
   FLOAT duv[NVert];
   VECT cell_centre = elem->barycenter;
   VECT vertex_pos;
   INT vert_index;

   for(k = 0 ; k < NVert ; k ++)
   {
      vert_index = elem->verts[k];
      vertex_pos = grid->vert_coord[vert_index];
      duv[k] = InerProduct(init_grad, vertex_pos - cell_centre);
      if(std::fabs(duv[k])>TOL) // 防止出现 duv \aprox 0 的情形. TOL defined in def.h
      {
         alpha = std::min(alpha,
                          std::max(0.0,
                                   std::max((_max_min_vt[vert_index][0]-
                                             cell_average)/duv[k],
                                            (_max_min_vt[vert_index][1]-
                                             cell_average)/duv[k])));
      }
      else
      {
         alpha = 1.0;
      }
      duv[k] = alpha * duv[k];
   }

   // limited_grad = alpha * init_grad;

   FLOAT pos, neg;
   FLOAT theta_pos, theta_neg;
   FLOAT res;
   INT ii, jj;
   FLOAT f[2], L_1[2], L_2[2];
   // VECT mid_pos;
   VECT cell_cen = elem->barycenter;

   res = 0;
   for(k = 0; k < NEdge ; k ++){
      res += duv[k];
   }
   if(res != 0){
      // step.1 计算 modified increment.
      pos = 0;
      neg = 0;
      for(k = 0; k < NEdge ; k ++){
         pos += std::max(0.0, duv[k]);
         neg += std::max(0.0, -duv[k]);
      }
      theta_pos = std::min(1.0, neg/pos);
      theta_neg = std::min(1.0, pos/neg);
      for(k = 0 ; k < NEdge ; k ++){
         duv[k] =
            theta_pos * std::max(0.0, duv[k]) -
            theta_neg * std::max(0.0, -duv[k]);
      }
   } // end of if

   ii = 0;
   jj = 1;
   f[0] = duv[ii];
   f[1] = duv[jj];
   // mid_pos = elem->GetEdgeMidPosition(ii);
   // L_1[0] = mid_pos.x - cell_cen.x;
   // L_2[0] = mid_pos.y - cell_cen.y;
   // mid_pos = elem->GetEdgeMidPosition(jj);
   // L_1[1] = mid_pos.x - cell_cen.x;
   // L_2[1] = mid_pos.y - cell_cen.y;
   vert_index = elem->verts[ii];
   vertex_pos = grid->vert_coord[vert_index];
   L_1[0] = vertex_pos.x - cell_cen.x;
   L_2[0] = vertex_pos.y - cell_cen.y;
   vert_index = elem->verts[jj];
   vertex_pos = grid->vert_coord[vert_index];
   L_1[1] = vertex_pos.x - cell_cen.x;
   L_2[1] = vertex_pos.y - cell_cen.y;

   limited_grad.x =
      (L_2[1] * f[0] - L_2[0] * f[1]) /
      (L_1[0] * L_2[1] - L_2[0] * L_1[1]);
   limited_grad.y =
      (- L_1[1] * f[0] + L_1[0] * f[1]) /
      (L_1[0] * L_2[1] - L_2[0] * L_1[1]);

   return;
}

// ===================================================================

void GRAD_LIMITER_ORIGINAL_BARTH::apply_limiter(ELEMENT * elem, const VECT & init_grad,
                                                VECT & limited_grad){

   GRID * grid = elem->GetGrid();
   FLOAT c_max, c_min;
   FLOAT nb_max, nb_min;
   FLOAT duv_max, duv_min;
   ELEMENT * nb_elem;
   INT c_index = elem->index;
   INT nb_index;
   VECT cell_centre = elem->barycenter;
   VECT edge_mid_pos;
   FLOAT cell_average = (*_solution)[elem->index];
   FLOAT alpha = 1.0;
   FLOAT duv;
   INT j;

   get_min_max(elem, c_max, c_min);

   for(j = 0 ; j < NEdge ; j ++){

      nb_index = elem->nghbs[j];
      nb_elem = grid->elem[nb_index];
      get_min_max(nb_elem, nb_max, nb_min);

      duv_max = std::min(nb_max, c_max);
      duv_min = std::max(nb_min, c_min);

      edge_mid_pos = elem->GetEdgeMidPosition(j);
      duv = InerProduct(init_grad, edge_mid_pos - cell_centre);
      if(std::fabs(duv)>TOL) // 防止出现 duv \aprox 0 的情形. TOL defined in def.h
      {
         alpha = std::min(alpha,
                          std::max(0.0,
                                   std::max((duv_max - cell_average)/duv,
                                            (duv_min - cell_average)/duv)));
      }
      else
      {
         alpha = 1.0;
      }
   }

   limited_grad = alpha * init_grad;

   return;
}

void GRAD_LIMITER_ORIGINAL_BARTH::get_min_max(ELEMENT * elem, FLOAT & max, FLOAT & min){
   INT c_index = elem->index;
   INT nb_index;
   INT j;
   bool is_bd = elem->is_ghost;
   BOUNDARY_ELEMENT * bd_elem;
   ELEMENT * original_elem;

   if(is_bd){
      bd_elem = static_cast <BOUNDARY_ELEMENT*> (elem);
      c_index = bd_elem->original_elem;
      original_elem = elem->GetGrid()->elem[c_index];
   } else {
      c_index = elem->index;
   }

   max = (*_solution)[c_index];
   min = (*_solution)[c_index];

   if(is_bd){
      for(j = 0 ; j < NEdge ; j ++){
         nb_index = original_elem->nghbs[j];
         max = std::max(max , (*_solution)[nb_index]);
         min = std::min(min , (*_solution)[nb_index]);
      }
   } else {
      for(j = 0 ; j < NEdge ; j ++){
         nb_index = elem->nghbs[j];
         max = std::max(max , (*_solution)[nb_index]);
         min = std::min(min , (*_solution)[nb_index]);
      }
   }

   return;
}

// ===================================================================

void GRAD_LIMITER_MLP_ORIGINAL::apply_limiter(ELEMENT * elem, const VECT & init_grad,
                                              VECT & limited_grad){
   INT j,k;
   GRID * grid = _solution->GetGrid();
   FLOAT cell_average = (*_solution)[elem->index];
   VECT cell_centre = elem->barycenter;
   INT cell_index = elem->index;
   INT next_index;
   ELEMENT * next_elem;
   ELEMENT * pre_elem;
   BOUNDARY_ELEMENT * next_elem_bd;
   FLOAT alpha = 1.0;
   FLOAT duv;
   VECT vertex_pos;
   INT vert_index;
   FLOAT vertex_min, vertex_max;
   INT pre_edge_index, next_edge_index;

   for(j = 0 ; j < NVert ; j ++){
      vertex_min = cell_average;
      vertex_max = cell_average;
      next_elem = elem;
      next_edge_index = (j + 1) % NVert;
      next_index = elem->nghbs[next_edge_index];
      pre_elem = next_elem;
      next_elem = grid->elem[next_index];
      pre_edge_index = pre_elem->GetNghbrEdgeIndex(next_edge_index);
      if(next_elem->is_ghost){
         next_elem_bd = static_cast <BOUNDARY_ELEMENT*> (next_elem);
         next_index = next_elem_bd->original_elem;
         next_elem = grid->elem[next_index];
         pre_edge_index = next_elem_bd->GetOriginalVertIndex(pre_edge_index);
      }

      while(next_index != cell_index){
         // 收集最大最小值信息
         vertex_min = std::min(vertex_min, (*_solution)[next_index]);
         vertex_max = std::max(vertex_max, (*_solution)[next_index]);

         // 确定下一个遍历单元
         next_edge_index = (pre_edge_index + 2) % NVert;
         next_index = next_elem->nghbs[next_edge_index];
         pre_elem = next_elem;
         next_elem = grid->elem[next_index];
         pre_edge_index = pre_elem->GetNghbrEdgeIndex(next_edge_index);
         if(next_elem->is_ghost){
            next_elem_bd = static_cast <BOUNDARY_ELEMENT*> (next_elem);
            next_index = next_elem_bd->original_elem;
            next_elem = grid->elem[next_index];
            pre_edge_index = next_elem_bd->GetOriginalVertIndex(pre_edge_index);
         }
      } // end of while

      // cal alpha.
      vert_index = elem->verts[j];
      vertex_pos = grid->vert_coord[vert_index];
      duv = InerProduct(init_grad, vertex_pos - cell_centre);
      if(std::fabs(duv)>TOL) // 防止出现 duv \aprox 0 的情形. TOL defined in def.h
      {
         alpha = std::min(alpha,
                          std::max(0.0,
                                   std::max((vertex_max - cell_average)/duv,
                                            (vertex_min - cell_average)/duv)));
      }

   } // each vertex

   // limite the gradiant
   limited_grad = alpha * init_grad;

   return;
}
