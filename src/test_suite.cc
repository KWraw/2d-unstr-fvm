/*
  File:      test_suite.cc
  Purpose:   测试所用函数

  Author: <kw.xu@hotmail.com>
  Last modified: "2016-04-26 16:38:24 kaiwen"
*/

#include "test_suite.h"
#include "def.h"
#include <cmath>

#define _USE_MATH_DEFINES       // for M_PI

// ===================================================================
// Cases
// ===================================================================
// ===================================================================
// Double Sine
namespace TEST_SUITE_ADV_DOUBLESINE{
   // 具体的函数.
   static FLOAT _init_function(const VECT & pos){
      // double-sine
      // u_0 = sin(2 \pi x)sin(2 \pi y)
      return std::sin(2 * M_PI * pos.x) * std::sin(2 * M_PI * pos.y);
   }
   static VECT _field_function(const VECT & pos){
      // (u, v) = (1, 2);
      VECT vc(1.0, 2.0);
      // VECT vc(1.0, 1.0);
      return vc;
   }
   static VECT _flux_function(const VECT & pos, FLOAT u){
      // const field
      // (f, g) = (1, 2) * u
      VECT vc(_field_function(pos));
      vc.x = vc.x * u;
      vc.y = vc.y * u;
      return vc;
   }
   static FLOAT _exact_function(const VECT & pos, FLOAT t){
      // 推进的时间设定为某个整数
      if(std::fabs(t - (int)t) > TOL){
         std::cerr<<"TEST_SUITE, ERROR! time must be integer."<<std::endl;
         exit(1);
      }
      return _init_function(pos);
   }

   // 封装的 functor.
   INIT_FUNC_2D init_function(_init_function, "u_0 = sin(2 \\pi x)sin(2 \\pi y)");
   FLUX_FUNC_2D flux_function(_flux_function, "(f, g) = (1, 2) * u");
   FIELD_FUNC_2D field_function(_field_function, "(u, v) = (1, 2)");
   EXACT_FUNC_2D exact_function(_exact_function, "u = sin(2 \\pi x)sin(2 \\pi y)");
}

// ===================================================================
// Square : with discontinuouty.
namespace TEST_SUITE_ADV_SQUARE{
   // 具体的函数.
   static FLOAT _init_function(const VECT & pos){

      // Remove the restriction of integer end time.
      FLOAT pos_x = pos[0];
      FLOAT pos_y = pos[1];

      // Some kind of "periodic"
      if(pos_x<0){
         do{
            pos_x += 1.0;
         } while (pos_x < 0);
      } else {
         while (pos_x - 1.0 > 0)
            pos_x -= 1.0;
      }
      if(pos_y<0){
         do{
            pos_y += 1.0;
         } while (pos_y < 0);
      } else {
         while (pos_y - 1.0 > 0)
            pos_y -= 1.0;
      }

      FLOAT val;
      if((pos_x<0.75&&pos_x>0.25) && (pos_y<0.75&&pos_y>0.25))
         val = 1;
      else
         val = 0;
      return val;
   }
   static VECT _field_function(const VECT & pos){
      // (u, v) = (1, 2);
      VECT vc(1.0, 2.0);
      // VECT vc(-1.0, -2.0);
      // VECT vc(1.0, 1.0);
      return vc;
   }
   static VECT _flux_function(const VECT & pos, FLOAT u){
      // const field
      // (f, g) = (1, 2) * u
      VECT vc(_field_function(pos));
      vc.x = vc.x * u;
      vc.y = vc.y * u;
      return vc;
   }
   static FLOAT _exact_function(const VECT & pos, FLOAT t){
      // 推进的时间设定为某个整数
      // 取消时间为整数的限制.
      // if(std::fabs(t - (int)t) > TOL){
      //    std::cerr<<"TEST_SUITE, ERROR! time must be integer."<<std::endl;
      //    exit(1);
      // }
      VECT new_pos(pos[0] - 1.0 * t, pos[1] - 2.0 * t);
      // VECT new_pos(pos[0] + 1.0 * t, pos[1] + 2.0 * t);
      // return _init_function(pos);
      return _init_function(new_pos);
   }

   // 封装的 functor.
   INIT_FUNC_2D init_function(_init_function,
                              "Square locate in the middle of the area");
   FLUX_FUNC_2D flux_function(_flux_function, "(f, g) = (1, 2) * u");
   FIELD_FUNC_2D field_function(_field_function, "(u, v) = (1, 2)");
   EXACT_FUNC_2D exact_function(_exact_function,
                                "Square locate in the middle of the area");
}

// ===================================================================
// Solid body rotation.
namespace TEST_SUITE_ADV_SOLID_BODY_ROTATION{
   // 具体的函数.
   static FLOAT _init_function(const VECT & pos){
      // double-sine
      // u_0 = sin(2 \pi x)sin(2 \pi y)
      // return std::sin(2 * M_PI * pos.x) * std::sin(2 * M_PI * pos.y);
      FLOAT val;
      VECT pos_0(0.5, 0.25);
      FLOAT norm_0 = 0.15;
      VECT vr = pos - pos_0;
      FLOAT norm = vr.GetNorm();
      if(norm > norm_0)
         val = 0.0;
      else
         val = 1.0 - norm/norm_0;
      return val;
   }
   static VECT _field_function(const VECT & pos){
      // (u, v) = (1, 2);
      // VECT vc(1.0, 2.0);
      // (u, v) = (-(y-0.5), (x-0.5)) * 2 * PI
      VECT vc(-pos[1] + 0.5, pos[0] - 0.5);
      // VECT vc(1.0, 1.0);
      return 2 * M_PI * vc;
   }
   static VECT _flux_function(const VECT & pos, FLOAT u){
      // const field
      // (f, g) = (1, 2) * u
      VECT vc(_field_function(pos));
      vc.x = vc.x * u;
      vc.y = vc.y * u;
      return vc;
   }
   static FLOAT _exact_function(const VECT & pos, FLOAT t){
      // 推进的时间设定为某个整数
      if(std::fabs(t - (int)t) > TOL){
         std::cerr<<"TEST_SUITE, ERROR! time must be integer."<<std::endl;
         exit(1);
      }
      return _init_function(pos);
   }
   // 封装的 functor.
   INIT_FUNC_2D init_function(_init_function, "Single Cone");
   FLUX_FUNC_2D flux_function(_flux_function, "(f, g) = (-y+0.5, x-0.5) * u");
   FIELD_FUNC_2D field_function(_field_function, "(u, v) = (-y+0.5, x-0.5)");
   EXACT_FUNC_2D exact_function(_exact_function,
                                "Single Cone");
}

// ===================================================================
// class TEST_SUITE
// ===================================================================
TEST_SUITE::TEST_SUITE(TEST_SUITE_TYPE test_suite_type){
   _test_suite_type = test_suite_type;
   switch(test_suite_type){
      case TEST_SUITE_TYPE_ADV_DOUBLESINE:{
         // init_function(TEST_SUITE_ADV_DOUBLESINE::_init_function,
         //               "u_0 = sin(2 \\pi x)sin(2 \\pi y)");
         // flux_function(TEST_SUITE_ADV_DOUBLESINE::_flux_function, "(f, g) = (1, 2) * u");
         // field_function(TEST_SUITE_ADV_DOUBLESINE::_field_function, "(u, v) = (1, 2)");
         // exact_function(TEST_SUITE_ADV_DOUBLESINE::_exact_function,
         //                "u = sin(2 \\pi x)sin(2 \\pi y)");
         ptr_init_function = &TEST_SUITE_ADV_DOUBLESINE::init_function;
         ptr_flux_function = &TEST_SUITE_ADV_DOUBLESINE::flux_function;
         ptr_field_function= &TEST_SUITE_ADV_DOUBLESINE::field_function;
         ptr_exact_function= &TEST_SUITE_ADV_DOUBLESINE::exact_function;
         break;
      }
      case TEST_SUITE_TYPE_ADV_SQUARE:{
         // init_function(TEST_SUITE_ADV_SQUARE::_init_function,
         //               "Square locate in the middle of the area");
         // flux_function(TEST_SUITE_ADV_SQUARE::_flux_function, "(f, g) = (1, 2) * u");
         // field_function(TEST_SUITE_ADV_SQUARE::_field_function, "(u, v) = (1, 2)");
         // exact_function(TEST_SUITE_ADV_SQUARE::_exact_function,
         //                "Square locate in the middle of the area");
         ptr_init_function = &TEST_SUITE_ADV_SQUARE::init_function;
         ptr_flux_function = &TEST_SUITE_ADV_SQUARE::flux_function;
         ptr_field_function= &TEST_SUITE_ADV_SQUARE::field_function;
         ptr_exact_function= &TEST_SUITE_ADV_SQUARE::exact_function;
         break;
      }
      case TEST_SUITE_TYPE_ADV_SOLID_BODY_ROTATION:{
         // init_function(TEST_SUITE_ADV_SOLID_BODY_ROTATION::_init_function,
         //               "Single Cone");
         // flux_function(TEST_SUITE_ADV_SOLID_BODY_ROTATION::_flux_function,
         //               "(f, g) = ((-y+0.5), (x-0.5)) * u");
         // field_function(TEST_SUITE_ADV_SOLID_BODY_ROTATION::_field_function,
         //                "(u, v) = (-y+0.5, x-0.5)");
         // exact_function(TEST_SUITE_ADV_SOLID_BODY_ROTATION::_exact_function,
         //                "Same with initial condition, Single Cone");
         ptr_init_function = &TEST_SUITE_ADV_SOLID_BODY_ROTATION::init_function;
         ptr_flux_function = &TEST_SUITE_ADV_SOLID_BODY_ROTATION::flux_function;
         ptr_field_function= &TEST_SUITE_ADV_SOLID_BODY_ROTATION::field_function;
         ptr_exact_function= &TEST_SUITE_ADV_SOLID_BODY_ROTATION::exact_function;
         break;
      }
      default:{
         std::cerr<<"Unkown Test Suite type: "<<test_suite_type<<std::endl;
         exit(1);
      }
   }
}

void TEST_SUITE::ShowTestSuiteInfo(){
   std::cout<<"TestSuite Info:"<<std::endl
            <<"   Initial Function: "<<ptr_init_function->GetName()<<std::endl
            <<"   Flux Function:    "<<ptr_flux_function->GetName()<<std::endl
            <<"   Field Function:   "<<ptr_field_function->GetName()<<std::endl
            <<"   Exact Solution:   "<<ptr_exact_function->GetName()<<std::endl;
   return;
}
