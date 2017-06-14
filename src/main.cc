/*
  File:      main.cc
  Purpose:

  Author: <kw.xu@hotmail.com>
  Last modified: "2016-05-17 17:00:12 kaiwen"
*/

#include "boundary.h"
#include "def.h"
#include "grid.h"
#include "recon.h"
#include "solver.h"
#include "test_suite.h"
#include <string>
#include <iomanip>

int
main(int argc, char *argv[]){

   FLOAT cfl = 0.5;
   RK_TYPE rk_type = RK_TYPE_2;
   FLOAT t_end = 1.0;
   // FLOAT t_end = 0.0;
   INT max_time_step = 10000;

   std::string mesh_file("./mesh/rect_trian_80");
   // std::string mesh_file("./mesh/rect_trian_80_typeA");
   // std::string mesh_file("./mesh/gmsh_square80");
   // std::string mesh_file("./mesh/gmsh_square80_front");

   // std::string outfile("./dat/rect_trian_20_typeC.dat");
   // std::string outfile("./dat/square_wave_Barth_Original_gmsh80.dat");
   // std::string outfile("./dat/square_wave_initial_rect80.dat");
   // std::string outfile("./dat/solid_body_rotation_initial_rect80.dat");
   // std::string outfile("./dat/solid_body_rotation_Barth_Original_gmsh80.dat");
   std::string outfile("./dat/test.dat");

   TEST_SUITE_TYPE test_suite_type = TEST_SUITE_TYPE_ADV_DOUBLESINE;
   // TEST_SUITE_TYPE test_suite_type = TEST_SUITE_TYPE_ADV_SQUARE;
   // TEST_SUITE_TYPE test_suite_type = TEST_SUITE_TYPE_ADV_SOLID_BODY_ROTATION;

   FLOAT L1_error;

   /* ---- Generate Gird. ---- */

   GRID grid;
   grid.ReadGridFile(mesh_file + ".e", mesh_file + ".n");
   grid.Precondition();
   PERIODIC_BOUNDARY boundary(&grid);
   boundary.ExpendBoundary();

   /* ---- Create solution ---- */

   SOLUTION solution(&grid, &boundary);

   /* ---- Construct the reconstructer ---- */

   GRAD_PREDICT_LS grad_predicter(&solution);

   // GRAD_LIMITER_NONE grad_limiter(&solution);
   // GRAD_LIMITER_VERTEX_BARTH grad_limiter(&solution);
   // GRAD_LIMITER_BARTH_JESPERSEN grad_limiter(&solution);
   GRAD_LIMITER_MLP grad_limiter(&solution);
   // GRAD_LIMITER_MLP_MIDPOINT grad_limiter(&svolution);
   // GRAD_LIMITER_LCD grad_limiter(&solution);
   // GRAD_LIMITER_SHU_COCKBURN grad_limiter(&solution);
   // GRAD_LIMITER_MLP_ORIGINAL grad_limiter(&solution);
   // GRAD_LIMITER_MLP_LIST grad_limiter(&solution);
   // GRAD_LIMITER_ORIGINAL_BARTH grad_limiter(&solution);
   // GRAD_LIMITER_MLP_LIMITE_ON_VALUE grad_limiter(&solution);

   MUSCL_RECON recon(grad_predicter, grad_limiter, solution);

   /* ---- Set testing case ---- */

   TEST_SUITE test_suite(test_suite_type);

   /* ---- Initialize Solver ---- */

   SOLVER solver(solution,
                 recon,
                 rk_type,
                 cfl,
                 // functors
                 test_suite.ptr_init_function,
                 test_suite.ptr_flux_function,
                 test_suite.ptr_field_function,
                 test_suite.ptr_exact_function
      );

   solver.ShowSolverInfo();

   /* ---- Advance in time ---- */

   solver.AdvanceSolution(t_end, max_time_step);

   /* ---- L1 error ---- */

   std::cout << "Calculating the L1 error ..." << std::endl;
   L1_error = solver.GetL1Error();
   std::cout << std::scientific;
   std::cout << std::setprecision(3);
   std::cout << "L1 error: " << L1_error << std::endl;
   std::cout << std::flush;

   /* ---- Peak ---- */

   std::cout << "Calculation the Peak ..." << std::endl;
   FLOAT max_peak, min_peak;
   std::cout << std::fixed;
   std::cout << std::setprecision(3);
   solver.GetPeak(max_peak, min_peak);
   std::cout << "Max: " << max_peak << std::endl;
   std::cout << "Min: " << min_peak << std::endl;
   std::cout << std::flush;

   /* ---- Output final profile ---- */

   solution.ExportTecplot(outfile);

   return 0;
} /* main */
