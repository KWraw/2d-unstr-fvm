/*
  File:      solver.h
  Purpose:   二维非结构网格有限体积法求解器

  Author: <kw.xu@hotmail.com>
  Last modified: "2016-05-15 10:41:12 kaiwen"
*/

#ifndef _SOLVER_H_
#define _SOLVER_H_

#include <vector>
#include <string>
#include "def.h"
#include "grid.h"
#include "boundary.h"
#include "test_suite.h"
#include "recon.h"
#include "solution.h"

// 1. 定义 solution. (单独用一个头文件管理)
// 2. 定义 solver.

// ===================================================================

class SOLVER{
public:
   // 构造函数
   // INPUT
   // solution -- 数值解.
   // recon -- 重构方式.
   // rk_type -- 时间推进方式.
   // cfl -- 控制时间推进步的常数.
   // init_func_2d -- 初值
   // flux_function -- 真实通量函数 (以位置 pos 和值的大小为变量)
   // field_function -- 速度向量场 (以位置 pos 为变量)
   // exact_function -- 精确解.
   SOLVER(SOLUTION & solution,
          MUSCL_RECON & recon,
          RK_TYPE rk_type,
          FLOAT cfl,
          INIT_FUNC_2D * init_func_2d,
          FLUX_FUNC_2D * flux_function,
          FIELD_FUNC_2D * field_function,
          EXACT_FUNC_2D * exact_function
      );

   // 析构函数
   ~SOLVER();

   // 输出解法器信息.
   void ShowSolverInfo();

   // 时间推进, 推进到指定时间重点或指定的时间步.
   // INPUT
   // t_end -- 终止时间.
   // max_time_step -- 最大时间步数.
   void AdvanceSolution(FLOAT t_end, INT max_time_step);

   // 误差分析, L1
   FLOAT GetL1Error();

   // 误差分析, L^\infty
   FLOAT GetL8Error();

   // 输出最大值和最小值
   void GetPeak(FLOAT & max, FLOAT & min);

private:
   // 为 SOLUTION 附初值
   void _setInitial();

   // 获取时间推进步长
   FLOAT _getGridDT();

   // 计算重构值 (该函数在全场范围进行计算)
   // 1. 先做一下预处理 (MUSCL_RECON::Pre_Treatment())
   // 2. 依单元进行梯度预测.
   // 3. 依单元进行梯度限制.
   // 4. 计算单元边中点处的重构值.
   // (2 ~ 4 步集成在 MUSCL_RECON::Reconstruction 函数里面)
   // 5. 计算位于边界单元上的中点处外侧的重构值.
   void _calRecon();

   // 计算单元上的空间离散算子 (计算半离散形式的右端部分)
   // 数值通量的处理包含在这里部分. 使用迎风通量. (如有必要可以考虑其他形式的通量)
   // 1. 计算数值通量.
   // 2. 计算完整的右端项.
   FLOAT _calElemSpatDiscrete(ELEMENT * elem);

   // 更新 ghost cell 上的平均值.
   void _updateBoundaryValue();

   // 更新边界单元边界边上的外侧重构值. 直接放在重构步 _calRecon() 里面.
   // void _updataBoundaryRecon();

   // 时间推进 (暂时只实现了二阶的 RK, 但是为了扩展性仍然保留 rk_type)
   // 将 solution 推进一个时间步.
   void _rk2(FLOAT dt);

   // 获取精确的单元平均值 (需要用到相应精度的数值积分)
   // INPUT
   // elem -- 需要计算平均值的单元.
   // t -- 时间.
   FLOAT _getExactCellAverage(ELEMENT * elem, FLOAT t);

private:
   SOLUTION & _solution;
   MUSCL_RECON & _recon;
   GRID * _grid;
   BOUNDARY * _boundary;
   // 内部单元边界中点处的重构值
   std::vector<std::array<FLOAT, NEdge>> _recon_value;
   // 边界单元的重构值
   std::vector<FLOAT> _bd_recon_value;
   // 高阶时间推进的中间时间层. 根据输入的 GRID 大小以及时间推进方式在
   // SOLUTION 初始化时进行内存分配.
   SOLUTION * _solution_1;
   SOLUTION * _solution_2;
   INIT_FUNC_2D * _init_function;
   FLUX_FUNC_2D * _flux_function;
   FIELD_FUNC_2D * _field_function;
   EXACT_FUNC_2D * _exact_function;
   RK_TYPE _rk;
   FLOAT _cfl;
   FLOAT recon_time;
   FLOAT step_time_elapse;
};

#endif /* _SOLVER_H_ */
