/*
  File:      def.h
  Purpose:   常量以及常用工具函数

  Author: <kw.xu@hotmail.com>
  Last modified: "2016-04-27 11:28:13 kaiwen"
*/

#ifndef _DEF_H_
#define _DEF_H_

#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <stdio.h>
#include <cstdlib>

/********************************************************
 * Prints the error message, the stack trace, and exits
 * ******************************************************/
#define FatalError(s) printf("Fatal error '%s' at %s:%d\n",s,__FILE__,__LINE__)

enum {Dim = 2, NVert = 3, NEdge = 3};
enum {BDMARK = -1};
enum {CHAR_WIDTH = 15};
typedef double      FLOAT;
typedef signed int  INT;
// typedef FLOAT       BARY[NVert];
// typedef std::vector<FLOAT, NVert> BARY_COOD;
typedef std::array<FLOAT, NVert> BARY_COOD;
const FLOAT TOL = 1.0e-10;

/* ---- collection of enum types ---- */

enum TEST_SUITE_TYPE
{
   TEST_SUITE_TYPE_ADV_DOUBLESINE=0,
   TEST_SUITE_TYPE_ADV_SQUARE=1,
   TEST_SUITE_TYPE_ADV_SOLID_BODY_ROTATION=2
};

enum LIMITER_TYPE
{
   LIMITER_TYPE_NONE=0,
   LIMITER_TYPE_BARTH=1,
   LIMITER_TYPE_MLP=2
};

enum RK_TYPE
{
   RK_TYPE_2=0,
   RK_TYPE_3=1
};

enum INIT_GRAD_TYPE
{
   GRAD_LS = 0,                 // 最小二乘
   GRAD_GS = 1,                 // Green-Gauss 重构
   GRAD_GD = 2                  // 梯度重构.
};

enum INPUT_TYPE
{
   INPUT_TYPE_EASYMESH=0
};

enum OUTPUT_TYPE
{
   OUTPUT_TYPE_TECPLOT=0
   // OUTPUT_TYPE_VTK=1
};

// ===================================================================
// VECT 定义了二维向量的运算.
class VECT
{
public:
   VECT() : x(0.0), y(0.0) {};
   VECT(FLOAT _x, FLOAT _y) : x(_x), y(_y) {};
   VECT(const VECT & vc) : x(vc.x), y(vc.y) {};
   // ~VECT();
public:
   FLOAT x,y;
   // VECT & operator=(const VECT & vc){x=vc.x; y=vc.y;}
   VECT & operator=(const VECT & vc){
      x = vc.x;
      y = vc.y;
      return *this;
   }
   FLOAT & operator[](int e) {
      if(e == 0)
         return x;
      else if (e == 1)
         return y;
      else {
         std::cout<<"VECT, ERROR! Visiting inexistence. e = " <<e<<std::endl;
         exit(1);
      }
   }
   FLOAT operator[](int e) const {
      if(e == 0)
         return x;
      else if (e == 1)
         return y;
      else {
         std::cout<<"VECT, ERROR! Visiting inexistence. e = " <<e<<std::endl;
         exit(1);
      }
   }
   friend std::ostream & operator<<(std::ostream & os, VECT & vc);
   friend VECT operator+(const VECT & vc1, const VECT & vc2);
   friend VECT operator-(const VECT & vc1, const VECT & vc2);
   // friend FLOAT operator*(const VECT & vc1, const VECT & vc2);
   friend VECT operator*(FLOAT a, const VECT & vc);
   friend FLOAT OuterProduct(const VECT & vc1, const VECT & vc2);
   friend FLOAT InerProduct(const VECT & vc1, const VECT & vc2);
   FLOAT GetNorm(){return std::sqrt(x*x+y*y);}
   VECT Normalize(){return (1.0/GetNorm()) * (*this);}
};


#endif /* _DEF_H_ */
