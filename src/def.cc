/*
  File:      def.cc
  Purpose:

  Author: <kw.xu@hotmail.com>
  Last modified: "2016-04-21 22:14:06 kaiwen"
*/

#include "def.h"
#include <cmath>
#include <iostream>
#include <iomanip>

std::ostream & operator<<(std::ostream & os, VECT & vc){
   enum {WIDTH = 30};
   os << std::scientific
      << std::setprecision(16)
      << std::setw(WIDTH) << vc.x
      << std::setw(WIDTH) << vc.y
      << std::endl;
   return os;
}

VECT operator+(const VECT & vc1, const VECT & vc2){
   VECT vc(vc1.x + vc2.x, vc1.y + vc2.y);
   return vc;
}
VECT operator-(const VECT & vc1, const VECT & vc2){
   VECT vc(vc1.x - vc2.x, vc1.y - vc2.y);
   return vc;
}
VECT operator*(FLOAT a, const VECT & vc){
   VECT temp_vc(a*vc.x, a*vc.y);
   return temp_vc;
}
FLOAT OuterProduct(const VECT & vc1, const VECT & vc2){
   return vc1.x * vc2.y - vc2.x * vc1.y;
}
FLOAT InerProduct(const VECT & vc1, const VECT & vc2){
   return vc1.x * vc2.x + vc1.y * vc2.y;
}

