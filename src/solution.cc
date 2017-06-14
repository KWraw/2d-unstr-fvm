/*
  File:      solution.cc
  Purpose:

  Author: <kw.xu@hotmail.com>
  Last modified: "2016-05-14 16:08:36 kaiwen"
*/

#include "solution.h"
#include <fstream>
#include <iomanip>

// ===================================================================
SOLUTION::SOLUTION(GRID * grid, BOUNDARY * boundary) :
   _grid(grid), _boundary(boundary), _current_t(0.0)
{
   _value.resize(grid->n_elem + boundary->n_boundary);
}

void SOLUTION::ExportTecplot(const std::string & filename) const
{
   FLOAT pt_value;
   std::string title_str("Advection equation");

   std::vector<FLOAT> pt_val(_grid->n_vert, 0.0);
   std::vector<INT> pt_num(_grid->n_vert, 0);

   for(int i = 0; i<_grid->n_elem; i++){
      for(int j = 0; j<NVert; j++){
         // pt_val[_grid->elem[i]->verts[j]] += _dofs[i]->GetAverage();
         pt_val[_grid->elem[i]->verts[j]] += _value[i];
         pt_num[_grid->elem[i]->verts[j]]++;
      }
   }

   std::ofstream f_out;
   f_out.open(filename.c_str());
   if( !f_out ){
      std::cerr << "SOLUTION::ExportTecplot, ERROR: unable to open \""
                << filename << "\"\n";
      // return false;
      exit(1);
   }
   std::cout<<"SOLUTION::Export solution in Tecplot format, to file : "<<filename<<"...";

   // header.
   f_out<<"TITLE =\""<<title_str<<"\"\n"
        <<"VARIABLES = \"X\", \"Y\", \"U\"\n";

   // zone record.
   // zone control line.
   f_out<<"ZONE T=\"U\", F=FEPOINT, "
        <<"N="<<_grid->n_vert<<", "
        <<"E="<<_grid->n_elem<<", "
        <<"ET=TRIANGLE"<<std::endl;
   // zone data lines.
   // f_out<<std::setw(20)<<std::setprecision(15);
   f_out<<std::setprecision(15);
   for(int i = 0; i<_grid->n_vert; i++)
   {
      pt_value = pt_val[i]/pt_num[i];

      // if(pt_value > 0.5) {pt_value = 1.0;}
      // else {pt_value = 0.0;}

      f_out<<_grid->vert_coord[i][0]<<" "
           <<_grid->vert_coord[i][1]<<" "
           << pt_value <<std::endl;
   }

   // more data zone...

   // Connectivity lists
   for(int i = 0; i<_grid->n_elem; i++)
      f_out<<_grid->elem[i]->verts[0] + 1<<" "
           <<_grid->elem[i]->verts[1] + 1<<" "
           <<_grid->elem[i]->verts[2] + 1<<"\n";

   f_out.close();
   std::cout<<"Finished!"<<std::endl;
   return;
}
