/*
  File:      solver.cc
  Purpose:

  Author: <kw.xu@hotmail.com>
  Last modified: "2016-05-15 10:53:10 kaiwen"
*/

#include "solver.h"
#include <fstream>
#include <cmath>
#include <limits>
#include <iomanip>
#include <time.h>

// ===================================================================

SOLVER::SOLVER(SOLUTION & solution,
               MUSCL_RECON & recon,
               RK_TYPE rk_type,
               FLOAT cfl,
               INIT_FUNC_2D * init_func_2d,
               FLUX_FUNC_2D * flux_function,
               FIELD_FUNC_2D * field_function,
               EXACT_FUNC_2D * exact_function) :
   _solution(solution), _recon(recon),
   _rk(rk_type), _cfl(cfl),
   _init_function(init_func_2d), _flux_function(flux_function),
   _field_function(field_function), _exact_function(exact_function)
{
   _grid = _solution.GetGrid();
   _boundary = _solution.GetBoundary();
   _recon_value.resize(_grid->n_elem);
   _bd_recon_value.resize(_boundary->n_boundary);

   _solution_1 = 0;
   _solution_2 = 0;

   switch(_rk){
      case RK_TYPE_2:{
         _solution_1 = new SOLUTION(_grid,_boundary);
         break;
      }
      default:{
         std::cerr<<"SOLVER::Error! Unimplemented RK Type : "<<_rk<<std::endl;
         exit(1);
      }
   }
}

SOLVER::~SOLVER()
{
   if(_solution_1 != 0) delete _solution_1;
   if(_solution_2 != 0) delete _solution_2;
}

void SOLVER::ShowSolverInfo()
{
   std::cout<<"==================================================="<< std::endl
            <<" 2d-unstr-fvm-adv -- FVM solver for 2D Advection equations."<< std::endl;
   _grid->ShowInfo();
   _boundary->ShowInfo();
   // _dof_type->ShowInfo();
   // _quad->ShowInfo();
   // _rk.ShowInfo();
   // 关于用于自定义的 function
   std::cout<<std::endl
            <<"=========================="<< std::endl
            <<"User defined functions:"<< std::endl
            <<"   init function  -- "<< _init_function->GetName() << std::endl
            <<"   flux function  -- "<< _flux_function->GetName() << std::endl
            <<"   field function -- "<< _field_function->GetName()<<std::endl
            <<"=========================="<< std::endl;
   std::cout<<"==================================================="<< std::endl;
   return;
}


void SOLVER::AdvanceSolution(FLOAT t_end, INT max_time_step)
{
   INT time_step = 0;
   FLOAT dt;
   bool if_end = false;
   FLOAT time_start, time_end;

   //,--------------
   //| Step.1 附初值
   //`--------------

   _setInitial();

   //,----------------
   //| Step.2 时间推进
   //`----------------

   if(time_step > max_time_step || t_end < 0){
      if_end = true;
   }

   while(!if_end)
   {
      time_step++;
      if(time_step == max_time_step) if_end = true;

      time_start = std::clock();

      dt = _getGridDT();
      if(dt + _solution.GetCurrentTime() > t_end){
         if_end = true;
         dt = t_end - _solution.GetCurrentTime();
      }

      switch(_rk){
         case RK_TYPE_2:
            _rk2(dt);
            break;
         default:
            std::cerr<<"SOLVER, Error. Unimplemented time step type." << std::endl;
            exit(1);
      }

      time_end = std::clock();
      step_time_elapse = ((FLOAT)(time_end - time_start)) / CLOCKS_PER_SEC;

      // show step information
      std::cout<<std::fixed<<std::setprecision(1);
      std::cout<<"Progress : "
               <<"["<<std::setw(4)<<100.0 * _solution.GetCurrentTime() / t_end<<"\%]"
               <<", dt -- "
               <<std::scientific<<std::setw(7)<<dt
               <<", n_time_step : " << time_step;
      std::cout<<std::setprecision(3);
      std::cout<<", recon_time : " << recon_time << " (s)";
      std::cout<<", total step time : " << step_time_elapse;
      // <<std::endl;
      std::cout<<std::flush<<'\r';

      _solution.IncrementDT(dt);
   }
   std::cout<<std::endl;
}

FLOAT SOLVER::GetL1Error()
{
   FLOAT l1_error = 0;
   INT l;
   FLOAT c_time = _solution.GetCurrentTime();
   ELEMENT * elem;

   for(l = 0 ; l < _grid->n_elem ; l ++)
   {
      elem = _grid->elem[l];
      l1_error += elem->area *
         std::abs(_solution[l] - _getExactCellAverage(elem, c_time));
   }

   return l1_error;
}

FLOAT SOLVER::GetL8Error()
{
   FLOAT l8_error = 0;
   INT l;
   FLOAT c_time = _solution.GetCurrentTime();

   for(l = 0 ; l < _grid->n_elem ; l ++)
   {
      l8_error =
         std::max(l8_error,
                  std::abs(_solution[l] -
                           _getExactCellAverage(_grid->elem[l], c_time)));
   }

   return l8_error;
}

void SOLVER::GetPeak(FLOAT & max_peak, FLOAT & min_peak)
{

   INT l;

   max_peak = std::numeric_limits<float>::min();
   min_peak = std::numeric_limits<float>::max();

   for(l = 0 ; l < _grid->n_elem ; l ++)
   {
      max_peak = std::max(max_peak, _solution[l]);
      min_peak = std::min(min_peak, _solution[l]);
   }

   return;
}

void SOLVER::_setInitial()
{
   FLOAT weight = .3333333333333333333333333333333333L;
   INT n_face_quad = 3;
   FLOAT f_val;
   // FLOAT area;
   FLOAT res;
   VECT pos;

   for(int l = 0; l < _grid->n_elem; l++){
      // area = _grid->elem[l]->area;
      res = 0;
      for(int j = 0; j<n_face_quad; j++){
         pos = _grid->elem[l]->GetEdgeMidPosition(j);
         f_val = (*_init_function)(pos);
         // res += weight * f_val * area;
         res += weight * f_val;
      } // each quad point.

      // if(1){
      //    // only for square wave test
      //    if(res > 0.5) {res = 1.0;}
      //    else {res = 0;}
      // }
      _solution[l] = res;
   }

   return;
}

FLOAT SOLVER::_getGridDT(){
   FLOAT grid_dt = std::numeric_limits<float>::max();
   FLOAT elem_dt;
   FLOAT temp;
   VECT orient;
   VECT pos;
   ELEMENT * c_elem;
   INT l, k;

   for(l = 0 ; l < _grid->n_elem ; l ++)
   {
      elem_dt = std::numeric_limits<float>::max();
      c_elem = _grid->elem[l];
      for(k = 0 ; k < NEdge ; k ++)
      {
         orient = c_elem->orient[k];
         pos = c_elem->GetEdgeMidPosition(k);
         temp = (c_elem->area)/(c_elem->edge_len[k] *
                                InerProduct(orient, (*_field_function)(pos)));
         temp = std::abs(temp);
         elem_dt = std::min(elem_dt,temp);
      }
      grid_dt = std::min(grid_dt, elem_dt);
   }

   return grid_dt;
}

void SOLVER::_calRecon()
{
   INT l, k;
   ELEMENT * c_elem;
   FLOAT recon_value[NEdge];
   BOUNDARY_ELEMENT * bd_elem;
   INT original_edge_index;
   INT original_elem_index;
   clock_t time_start, time_end;

   // Step.1 预处理
   time_start = std::clock();
   _recon.Pre_Treatment();

   // Step.2 对内部单元单元边界中点进行重构.
   for(l = 0 ; l < _grid->n_elem ; l ++)
   {
      c_elem = _grid->elem[l];
      _recon.Reconstruction(c_elem, recon_value);
      for(k = 0 ; k < NEdge ; k ++) {_recon_value[l][k] = recon_value[k];}
   }

   time_end = std::clock();
   recon_time = ((FLOAT)(time_end - time_start)) / CLOCKS_PER_SEC;
   // recon_time = (time_end - time_start);

   // Step.3 处理边界边中点的外侧重构.
   for(l = 0 ; l < _boundary->n_boundary; l ++)
   {
      bd_elem = _boundary->boundary_elem[l];
      original_edge_index = bd_elem->GetOriginalVertIndex(2);
      original_elem_index = bd_elem->original_elem;
      _bd_recon_value[l] = _recon_value[original_elem_index][original_edge_index];
   }

   return;
}

FLOAT SOLVER::_calElemSpatDiscrete(ELEMENT * elem)
{
   INT k;
   FLOAT elem_spatial;
   FLOAT recon_ext, recon_int;
   INT boundary_index;
   ELEMENT * nbr_elem;
   VECT orient, pos;
   FLOAT num_flux, edge_len;

   elem_spatial = 0;
   for(k = 0; k < NEdge; k ++)
   {
      // find ext / int recon value.
      recon_int = _recon_value[elem->index][k];
      if(elem->nghbs[k] < _grid->n_elem)
      {                         // internal edge
         nbr_elem = _grid->elem[elem->nghbs[k]];
         recon_ext = _recon_value[nbr_elem->index][elem->GetNghbrEdgeIndex(k)];
      }
      else
      {                         // external edge
         boundary_index = _boundary->GetBoundaryIndex(elem->nghbs[k]);
         recon_ext = _bd_recon_value[boundary_index];
      }
      // add on
      orient = elem->orient[k];
      pos = elem->GetEdgeMidPosition(k);
      edge_len = elem->edge_len[k];
      // numrical flux
      if(InerProduct(orient, (*_field_function)(pos))>0)
      {
         num_flux = InerProduct(orient, (*_flux_function)(pos, recon_int));
      }
      else
      {
         num_flux = InerProduct(orient, (*_flux_function)(pos, recon_ext));
      }
      elem_spatial += num_flux * edge_len;
   } // each edge
   elem_spatial = - elem_spatial / elem->area;

   return elem_spatial;
}

void SOLVER::_updateBoundaryValue()
{
   INT l;
   INT elem_index, original_index;
   BOUNDARY_ELEMENT * bd_elem;

   for(l = 0; l < _boundary->n_boundary; l ++)
   {
      bd_elem = _boundary->boundary_elem[l];
      elem_index = _boundary->GetBoundaryIndexVersus(l);
      original_index = bd_elem->original_elem;
      _solution[elem_index] = _solution[original_index];
   }

   return;
}

// void SOLVER::_updataBoundaryRecon()
// {
//    INT l;
//    BOUNDARY_ELEMENT * bd_elem;
//    INT original_index;

//    for(l = 0 ; l < _boundary->n_boundary ; l ++)
//    {
//       bd_elem = _boundary->boundary_elem[l];
//       original_index = bd_elem->GetOriginalVertIndex(2);
//       _bd_recon_value[l] = _recon_value[bd_elem->original_elem][original_index];
//    }

//    return;
// }

void SOLVER::_rk2(FLOAT dt)
{
   INT l;
   ELEMENT * c_elem;
   INT elem_index;

   //,-------
   //| Step.1
   //`-------
   (*_solution_1) = _solution;

   //,-------
   //| Step.2
   //`-------
   // 计算完重构值之后单元平均值即可更新.
   _calRecon();
   // 更新边界的 recon.
   // _updataBoundaryRecon();
   // 利用重构值进行时间推进.
   for(l = 0 ; l < _grid->n_elem ; l ++)
   {
      c_elem = _grid->elem[l];
      elem_index = c_elem->index;
      _solution[elem_index] += dt * _calElemSpatDiscrete(c_elem);
   }
   // 更新边界
   _updateBoundaryValue();

   //,-------
   //| Step.3
   //`-------
   _calRecon();
   // _updataBoundaryRecon();
   for(l = 0 ; l < _grid->n_elem ; l ++)
   {
      c_elem = _grid->elem[l];
      elem_index = c_elem->index;
      _solution[elem_index] += dt * _calElemSpatDiscrete(c_elem);
      _solution[elem_index] =
         0.5 * _solution[elem_index] +
         0.5 * (*_solution_1)[elem_index];
   }

   return;
}

FLOAT SOLVER::_getExactCellAverage(ELEMENT * elem, FLOAT t)
{
   FLOAT exact_average;
   FLOAT weight = .3333333333333333333333333333333333L;
   INT n_face_quad = 3;
   FLOAT f_val;
   // FLOAT area;
   FLOAT res;
   VECT pos;
   // std::vector<FLOAT> rhs(n_base, 0.0);
   // for(int l = 0; l < _grid->n_elem; l++){
   // area = elem->area;
   // for(int i = 0; i<n_base; i++)
   //    rhs[i] = 0.0;
   res = 0;
   for(int j = 0; j<n_face_quad; j++){
      // pos = _u_h[l].GetFaceQuadPosition(j);
      pos = elem->GetEdgeMidPosition(j);
      f_val = (*_exact_function)(pos, t);
      // weight = quad->face_quad_pt[j]->weight;
      // for(int i = 0; i<n_base; i++)
      //    rhs[i] +=weight * f_val * _u_h[l].GetFaceQuadBaseValue(i, j) * area;
      // res += weight * f_val * area
      res += f_val;
   } // each quad point.
   // _solution[l] = res;
   exact_average = weight * res;
   // }

   return exact_average;
}
