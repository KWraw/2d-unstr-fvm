/*
  File:      test_suite.h
  Purpose:   定义(标量函数)测试算例的, 初值, 通量, 场函数, 精确解.

  Author: <kw.xu@hotmail.com>
  Last modified: "2016-04-23 14:45:04 kaiwen"
*/

#ifndef _TEST_SUITE_H_
#define _TEST_SUITE_H_

#include "def.h"

// ===================================================================
// typedef FLOAT (* _INIT_FUNC_2D)(FLOAT xx, FLOAT yy);
typedef FLOAT (* _INIT_FUNC_2D)(const VECT & pos);
// define a functor.
class INIT_FUNC_2D{
public:
   INIT_FUNC_2D(_INIT_FUNC_2D init_function, const std::string & init_name)
      : _init_function(init_function), _init_name(init_name) {};
   FLOAT operator()(const VECT & pos) { return _init_function(pos); }
   const std::string & GetName(){return _init_name;}
private:
   _INIT_FUNC_2D _init_function;
   std::string _init_name;
};

// ===================================================================
// 接受空间变量 x,y, 返回 VECT.
// typedef VECT (* _FLUX_FUNC_2D)(FLOAT xx, FLOAT yy);
typedef VECT (* _FLUX_FUNC_2D)(const VECT & pos, FLOAT u);
class FLUX_FUNC_2D{
public:
   FLUX_FUNC_2D(_FLUX_FUNC_2D flux_function, const std::string & flux_name)
      : _flux_function(flux_function), _flux_name(flux_name) {};
   // VECT operator()(FLOAT x, FLOAT y) {return _flux_function(x, y); }
   VECT operator()(const VECT & pos, FLOAT u) {return _flux_function(pos, u); }
   const std::string & GetName(){return _flux_name;}
private:
   _FLUX_FUNC_2D _flux_function;
   std::string _flux_name;
};

// ===================================================================
// 速度场 (u,v),
// 常见的测试算例:
// 1. (u,v) = (1,2)
// 2. (u,v) = (-(y-0.5),x-0.5).
typedef VECT (* _FIELD_FUNC_2D)(const VECT & pos);
class FIELD_FUNC_2D{
public:
   FIELD_FUNC_2D(_FIELD_FUNC_2D field_function, const std::string & field_name)
      : _field_function(field_function), _field_name(field_name){};
   VECT operator()(const VECT & pos) {return _field_function(pos);}
   const std::string & GetName(){return _field_name;}
private:
   _FIELD_FUNC_2D _field_function;
   std::string _field_name;
};

// ===================================================================
typedef FLOAT (* _EXACT_FUNC_2D)(const VECT & pos, FLOAT t);
class EXACT_FUNC_2D{
public:
   EXACT_FUNC_2D(_EXACT_FUNC_2D exact_function, const std::string & exact_name)
      : _exact_function(exact_function), _exact_name(exact_name){};
   FLOAT operator()(const VECT & pos, FLOAT t) {return _exact_function(pos, t);}
   const std::string & GetName(){return _exact_name;}
private:
   _EXACT_FUNC_2D _exact_function;
   std::string _exact_name;
};


// class TEST_SUITE{
// };

// namespace TEST_SUITE{
//    extern INIT_FUNC_2D init_function;
//    extern FLUX_FUNC_2D flux_function;
//    extern FIELD_FUNC_2D field_function;
//    extern EXACT_FUNC_2D exact_function;
// }

// move to def.hpp
// enum TEST_SUITE_TYPE
// {TEST_SUITE_TYPE_ADV_DOUBLESINE=0, TEST_SUITE_TYPE_ADV_SQUARE=1,
//  TEST_SUITE_TYPE_ADV_SOLID_BODY_ROTATION=2};

class TEST_SUITE{
public:
   TEST_SUITE(TEST_SUITE_TYPE test_suite_type);
   INIT_FUNC_2D * ptr_init_function;
   FLUX_FUNC_2D * ptr_flux_function;
   FIELD_FUNC_2D * ptr_field_function;
   EXACT_FUNC_2D * ptr_exact_function;
   void ShowTestSuiteInfo();
private:
   TEST_SUITE_TYPE _test_suite_type;
};

// namespace TEST_SUITE_ADV_DOUBLESINE{
//    extern INIT_FUNC_2D init_function;
//    extern FLUX_FUNC_2D flux_function;
//    extern FIELD_FUNC_2D field_function;
//    extern EXACT_FUNC_2D exact_function;
// }
// namespace TEST_SUITE_ADV_SQUARE{
//    extern INIT_FUNC_2D init_function;
//    extern FLUX_FUNC_2D flux_function;
//    extern FIELD_FUNC_2D field_function;
//    extern EXACT_FUNC_2D exact_function;
// }
// namespace TEST_SUITE_ADV_SOLID_BODY_ROTATION{
//    extern INIT_FUNC_2D init_function;
//    extern FLUX_FUNC_2D flux_function;
//    extern FIELD_FUNC_2D field_function;
//    extern EXACT_FUNC_2D exact_function;
// }


#endif /* _TEST_SUITE_H_ */
