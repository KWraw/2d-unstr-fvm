/*
  File:      recon.h
  Purpose:   定义二维非结构网格上的 MUSCL 重构过程.
             包含梯度预测和限制器两个步骤

  Author: <kw.xu@hotmail.com>
  Last modified: "2016-05-14 23:06:06 kaiwen"
*/

#ifndef _RECON_H_
#define _RECON_H_

#include "def.h"
// #include "solver.h"
#include "solution.h"
#include "grid.h"
#include "boundary.h"
#include <forward_list>

// ===================================================================

class GRAD_PREDICT{
public:
   // 构造函数
   // INPUT
   // solution -- 数值解.
   GRAD_PREDICT(SOLUTION * solution) : _solution(solution) {};
   virtual void init_grad(ELEMENT * elem, VECT & grad) = 0;
protected:
   SOLUTION * _solution;
};

class GRAD_LIMITER{
public:
   // 构造函数
   // INPUT
   // solution -- 数值解.
   GRAD_LIMITER(SOLUTION * solution) : _solution(solution) {};
   // 将限制器应用于某单元. (纯虚函数)
   // INPUT
   // elem -- 当前单元
   // init_grad -- 初始梯度值
   // OUTPUT
   // limited_grad -- 限制之后的梯度值.
   virtual void apply_limiter(ELEMENT * elem, const VECT & init_grad,
                              VECT & limited_grad) = 0;
   // 预处理步, 为 MLP 预留. (纯虚函数)
   virtual void pretreatment() = 0;
protected:
   SOLUTION * _solution;
};

// ===================================================================

class GRAD_PREDICT_LS : public GRAD_PREDICT{
public:
   // 构造函数, 调用 GRAD_PREDICT 的构造函数.
   GRAD_PREDICT_LS(SOLUTION * solution) : GRAD_PREDICT(solution) {};
   // 梯度预估.
   virtual void init_grad(ELEMENT * elem, VECT & grad);
};

// ===================================================================

class GRAD_LIMITER_MLP : public GRAD_LIMITER{
public:
   // 构造函数
   GRAD_LIMITER_MLP(SOLUTION * solution) : GRAD_LIMITER(solution)
   {
      _max_min_vt.resize(solution->GetGrid()->n_vert);
   }

   // 限制器
   virtual void apply_limiter(ELEMENT * elem, const VECT & init_grad,
                              VECT & limited_grad);
   // 预处理步
   virtual void pretreatment();
private:
   // 节点处收集的最大最小值信息.
   std::vector<std::array<FLOAT, 2> > _max_min_vt;
};

class GRAD_LIMITER_MLP_LIMITE_ON_VALUE : public GRAD_LIMITER{
public:
   // 构造函数
   GRAD_LIMITER_MLP_LIMITE_ON_VALUE(SOLUTION * solution) : GRAD_LIMITER(solution)
   {
      _max_min_vt.resize(solution->GetGrid()->n_vert);
   }

   // 限制器
   virtual void apply_limiter(ELEMENT * elem, const VECT & init_grad,
                              VECT & limited_grad);
   // 预处理步
   virtual void pretreatment();
private:
   // 节点处收集的最大最小值信息.
   std::vector<std::array<FLOAT, 2> > _max_min_vt;
};


class GRAD_LIMITER_MLP_MIDPOINT : public GRAD_LIMITER{
public:
   // 构造函数
   GRAD_LIMITER_MLP_MIDPOINT(SOLUTION * solution) : GRAD_LIMITER(solution)
   {
      _max_min_vt.resize(solution->GetGrid()->n_vert);
   }

   // 限制器
   virtual void apply_limiter(ELEMENT * elem, const VECT & init_grad,
                              VECT & limited_grad);
   // 预处理步
   virtual void pretreatment();
private:
   // 节点处收集的最大最小值信息.
   std::vector<std::array<FLOAT, 2> > _max_min_vt;
};


class GRAD_LIMITER_BARTH_JESPERSEN : public GRAD_LIMITER{
public:
   GRAD_LIMITER_BARTH_JESPERSEN(SOLUTION * solution) : GRAD_LIMITER(solution)
   {};

   // 限制器
   virtual void apply_limiter(ELEMENT * elem, const VECT & init_grad,
                              VECT & limited_grad);
   // 预处理步
   virtual void pretreatment(){};
};

class GRAD_LIMITER_VERTEX_BARTH : public GRAD_LIMITER{
public:
   GRAD_LIMITER_VERTEX_BARTH(SOLUTION * solution) : GRAD_LIMITER(solution) {};

   // 限制器
   virtual void apply_limiter(ELEMENT * elem, const VECT & init_grad,
                              VECT & limited_grad);
   // 预处理步
   virtual void pretreatment(){};
};


class GRAD_LIMITER_LCD : public GRAD_LIMITER{
public:
   GRAD_LIMITER_LCD(SOLUTION * solution) : GRAD_LIMITER(solution) {};

   // 限制器
   virtual void apply_limiter(ELEMENT * elem, const VECT & init_grad,
                              VECT & limited_grad);
   // 预处理步
   virtual void pretreatment(){};
};

class GRAD_LIMITER_SHU_COCKBURN : public GRAD_LIMITER{
public:
   GRAD_LIMITER_SHU_COCKBURN(SOLUTION * solution) : GRAD_LIMITER(solution) {};

   // 限制器
   virtual void apply_limiter(ELEMENT * elem, const VECT & init_grad,
                              VECT & limited_grad);
   // 预处理步
   virtual void pretreatment(){};
};

class GRAD_LIMITER_NONE : public GRAD_LIMITER{
public:
   // 构造函数
   GRAD_LIMITER_NONE(SOLUTION * solution) : GRAD_LIMITER(solution) {};

   // 限制器
   virtual void apply_limiter(ELEMENT * elem, const VECT & init_grad,
                              VECT & limited_grad){limited_grad = init_grad;}
   // 预处理步
   virtual void pretreatment(){};
};


class GRAD_LIMITER_ORIGINAL_BARTH : public GRAD_LIMITER{
public:
   GRAD_LIMITER_ORIGINAL_BARTH(SOLUTION * solution) : GRAD_LIMITER(solution) {};

   // 限制器
   virtual void apply_limiter(ELEMENT * elem, const VECT & init_grad,
                              VECT & limited_grad);
   // 预处理步
   virtual void pretreatment(){};

private:

   void get_min_max(ELEMENT * elem, FLOAT & max, FLOAT & min);
};

// ===================================================================

class GRAD_LIMITER_MLP_ORIGINAL : public GRAD_LIMITER{
public:
   // 构造函数
   GRAD_LIMITER_MLP_ORIGINAL(SOLUTION * solution) : GRAD_LIMITER(solution) {};

   // 限制器
   virtual void apply_limiter(ELEMENT * elem, const VECT & init_grad,
                              VECT & limited_grad);
   // 预处理步
   virtual void pretreatment(){};
};

// ===================================================================

class GRAD_LIMITER_MLP_LIST : public GRAD_LIMITER{
public:
   // 构造函数
   GRAD_LIMITER_MLP_LIST(SOLUTION * solution);

   // 限制器
   virtual void apply_limiter(ELEMENT * elem, const VECT & init_grad,
                              VECT & limited_grad);
   // 预处理步
   virtual void pretreatment(){};
private:
   std::vector<std::forward_list<INT> > vertex_neighbour_elem;
};


// ===================================================================

class MUSCL_RECON{
public:
   // 构造函数
   // INPUT
   // predicter -- 梯度预测算法
   // limiter -- 限制器
   // solution -- 数值解对象 (在 solver.h 里面定义, 定义了数值解, 可以从这里获得网格信息)
   MUSCL_RECON(GRAD_PREDICT & predicter,
               GRAD_LIMITER & limiter,
               SOLUTION & solution) :
      _predicter(predicter),
      _limiter(limiter),
      _solution(solution)
   {
      _init_grad.x = 0;
      _init_grad.y = 0;
      _limit_grad.x = 0;
      _limit_grad.y = 0;
   }
   // 重构操作的接口.
   // INPUT
   // elem -- 当前循环到的单元. (指针)
   // OUTPUT
   // recon_value -- 在三条边中点的重构值.
   void Reconstruction(ELEMENT * elem, FLOAT * recon_value);
   // 依次进行单元重构之前的预先处理.
   void Pre_Treatment();

private:
   GRAD_PREDICT & _predicter;
   GRAD_LIMITER & _limiter;
   SOLUTION & _solution;
   VECT _init_grad, _limit_grad;
};

#endif /* _RECON_H_ */
