/*
  File:      solution.h
  Purpose:   
  
  Author: <kw.xu@hotmail.com>
  Last modified: "2016-04-25 22:54:45 kaiwen"
*/

#ifndef _SOLUTION_H_
#define _SOLUTION_H_

#include "def.h"
#include "grid.h"
#include "boundary.h"

class SOLUTION{
public:
   // 构造函数
   // INPUT
   // grid -- 内部网格
   // boundary -- 边界
   SOLUTION(GRID * grid, BOUNDARY * boundary);
   // 取值函数
   // INPUT
   // i -- 单元编号.
   FLOAT & operator[](int i){return _value[i];}
   // 拷贝构造函数
   SOLUTION & operator=(const SOLUTION & solution_c){
      _grid = solution_c._grid;
      _boundary = solution_c._boundary;
      _value = solution_c._value;
      _current_t = solution_c._current_t;
      return *this;
   }
   // 获取网格
   GRID * GetGrid(){return _grid;}
   // 获取边界
   BOUNDARY * GetBoundary(){return _boundary;}
   // 输出为 Tecplot 格式
   void ExportTecplot(const std::string & filename) const;
   // 当前时间
   FLOAT & GetCurrentTime(){return _current_t;}
   // 时间推进
   void IncrementDT(FLOAT dt){_current_t += dt;}
private:
   GRID * _grid;
   BOUNDARY * _boundary;
   std::vector<FLOAT> _value;
   // _current_t set to 0 at initializition.
   FLOAT _current_t;
};


#endif /* _SOLUTION_H_ */
