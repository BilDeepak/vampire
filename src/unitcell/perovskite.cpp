//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>

// Vampire headers
#include "unitcell.hpp"

// unitcell module headers
#include "internal.hpp"

namespace unitcell{
namespace internal{

//------------------------------------------------------------------------------
// Generates an idealised spinel structure
//
//------------------------------------------------------------------------------
void build_perovskite(unitcell::unit_cell_t& unit_cell){

   // Set basic unit cell properties
   unit_cell.dimensions[0] = 1.0;
	unit_cell.dimensions[1] = 1.003115;
	unit_cell.dimensions[2] = 1.412102;

	unit_cell.shape[0][0] = 1;
	unit_cell.shape[0][1] = 0;
	unit_cell.shape[0][2] = 0;

	unit_cell.shape[1][0] = 0;
	unit_cell.shape[1][1] = 1.003115;
	unit_cell.shape[1][2] = 0;

	unit_cell.shape[2][0] = 0;
	unit_cell.shape[2][1] = 0;
	unit_cell.shape[2][2] = 1.412102;

	unit_cell.lcsize = 20;
	unit_cell.hcsize = 8;
	unit_cell.interaction_range = 1;
	unit_cell.atom.resize(20);
	unit_cell.surface_threshold = 8;

	//-----------------------------
	unit_cell.atom[0].x   = 0;
	unit_cell.atom[0].y   = 0;
	unit_cell.atom[0].z   = 0;
	unit_cell.atom[0].mat = 0;
	unit_cell.atom[0].lc  = 0;
	unit_cell.atom[0].hc  = 0;
	unit_cell.atom[0].ni  = 6;
   unit_cell.atom[0].nm  = false;
	//-----------------------------
	unit_cell.atom[1].x   = 0;
	unit_cell.atom[1].y   = 0;
	unit_cell.atom[1].z   = 0.50;
	unit_cell.atom[1].mat = 0;
	unit_cell.atom[1].lc  = 1;
	unit_cell.atom[1].hc  = 4;
	unit_cell.atom[1].ni  = 6;
   unit_cell.atom[1].nm  = false;
	//-----------------------------
	unit_cell.atom[2].x   = 0.50;
	unit_cell.atom[2].y   = 0.50;
	unit_cell.atom[2].z   = 0.50;
	unit_cell.atom[2].mat = 0;
	unit_cell.atom[2].lc  = 4;
	unit_cell.atom[2].hc  = 4;
	unit_cell.atom[2].ni  = 6;
   unit_cell.atom[2].nm  = false;
	//-----------------------------
	unit_cell.atom[3].x   = 0.50;
	unit_cell.atom[3].y   = 0.50;
	unit_cell.atom[3].z   = 0;
	unit_cell.atom[3].mat = 0;
	unit_cell.atom[3].lc  = 3;
	unit_cell.atom[3].hc  = 0;
	unit_cell.atom[3].ni  = 6;
   unit_cell.atom[3].nm  = false;
	//-----------------------------
	unit_cell.atom[4].x   = 0.49247;
	unit_cell.atom[4].y   = 0.03403;
	unit_cell.atom[4].z   = 0.25;
	unit_cell.atom[4].mat = 1;
	unit_cell.atom[4].lc  = 4;
	unit_cell.atom[4].hc  = 2;
	unit_cell.atom[4].ni  = 16;
   unit_cell.atom[4].nm  = true;
	//-----------------------------
	unit_cell.atom[5].x   = 0.50753;
	unit_cell.atom[5].y   = 0.96597;
	unit_cell.atom[5].z   = 0.75;
	unit_cell.atom[5].mat = 1;
	unit_cell.atom[5].lc  = 5;
	unit_cell.atom[5].hc  = 6;
	unit_cell.atom[5].ni  = 16;
   unit_cell.atom[5].nm  = true;
	//-----------------------------
	unit_cell.atom[6].x   = 0.00753;
	unit_cell.atom[6].y   = 0.53403;
	unit_cell.atom[6].z   = 0.25;
	unit_cell.atom[6].mat = 1;
	unit_cell.atom[6].lc  = 6;
	unit_cell.atom[6].hc  = 2;
	unit_cell.atom[6].ni  = 16;
   unit_cell.atom[6].nm  = true;
   //-----------------------------
	unit_cell.atom[7].x   = 0.99247;
	unit_cell.atom[7].y   = 0.46597;
	unit_cell.atom[7].z   = 0.75;
	unit_cell.atom[7].mat = 1;
	unit_cell.atom[7].lc  = 7;
	unit_cell.atom[7].hc  = 6;
	unit_cell.atom[7].ni  = 16;
   unit_cell.atom[7].nm  = true;
	//-----------------------------
	unit_cell.atom[8].x   = 0.57590;
	unit_cell.atom[8].y   = 0.48582;
	unit_cell.atom[8].z   = 0.25;
	unit_cell.atom[8].mat = 2;
	unit_cell.atom[8].lc  = 8;
	unit_cell.atom[8].hc  = 2;
	unit_cell.atom[8].ni  = 12;
   unit_cell.atom[8].nm  = true;
	//-----------------------------
	unit_cell.atom[9].x   = 0.42410;
	unit_cell.atom[9].y   = 0.51418;
	unit_cell.atom[9].z   = 0.75;
	unit_cell.atom[9].mat = 2;
	unit_cell.atom[9].lc  = 9;
	unit_cell.atom[9].hc  = 6;
	unit_cell.atom[9].ni  = 12;
   unit_cell.atom[9].nm  = true;
	//-----------------------------
	unit_cell.atom[10].x   = 0.92410;
	unit_cell.atom[10].y   = 0.98582;
	unit_cell.atom[10].z   = 0.25;
	unit_cell.atom[10].mat = 2;
	unit_cell.atom[10].lc  = 10;
	unit_cell.atom[10].hc  = 2;
	unit_cell.atom[10].ni  = 12;
   unit_cell.atom[10].nm  = true;
	//-----------------------------
	unit_cell.atom[11].x   = 0.07590;
	unit_cell.atom[11].y   = 0.01418;
	unit_cell.atom[11].z   = 0.75;
	unit_cell.atom[11].mat = 2;
	unit_cell.atom[11].lc  = 11;
	unit_cell.atom[11].hc  = 6;
	unit_cell.atom[11].ni  = 12;
   unit_cell.atom[11].nm  = true;
	//-----------------------------
	unit_cell.atom[12].x   = 0.21972;
	unit_cell.atom[12].y   = 0.28406;
	unit_cell.atom[12].z   = 0.04114;
	unit_cell.atom[12].mat = 2;
	unit_cell.atom[12].lc  = 12;
	unit_cell.atom[12].hc  = 1;
	unit_cell.atom[12].ni  = 12;
   unit_cell.atom[12].nm  = true;
	//-----------------------------
	unit_cell.atom[13].x   = 0.78028;
	unit_cell.atom[13].y   = 0.71594;
	unit_cell.atom[13].z   = 0.54114;
	unit_cell.atom[13].mat = 2;
	unit_cell.atom[13].lc  = 13;
	unit_cell.atom[13].hc  = 5;
	unit_cell.atom[13].ni  = 12;
   unit_cell.atom[13].nm  = true;
	//-----------------------------
	unit_cell.atom[14].x   = 0.78028;
	unit_cell.atom[14].y   = 0.71594;
	unit_cell.atom[14].z   = 0.95886;
	unit_cell.atom[14].mat = 2;
	unit_cell.atom[14].lc  = 14;
	unit_cell.atom[14].hc  = 7;
	unit_cell.atom[14].ni  = 12;
   unit_cell.atom[14].nm  = true;
	//-----------------------------
	unit_cell.atom[15].x   = 0.21972;
	unit_cell.atom[15].y   = 0.28406;
	unit_cell.atom[15].z   = 0.45886;
	unit_cell.atom[15].mat = 2;
	unit_cell.atom[15].lc  = 15;
	unit_cell.atom[15].hc  = 3;
	unit_cell.atom[15].ni  = 12;
   unit_cell.atom[15].nm  = true;
	//-----------------------------
	unit_cell.atom[16].x   = 0.28028;
	unit_cell.atom[16].y   = 0.78406;
	unit_cell.atom[16].z   = 0.45886;
	unit_cell.atom[16].mat = 2;
	unit_cell.atom[16].lc  = 16;
	unit_cell.atom[16].hc  = 3;
	unit_cell.atom[16].ni  = 12;
   unit_cell.atom[16].nm  = true;
	//-----------------------------
	unit_cell.atom[17].x   = 0.71972;
	unit_cell.atom[17].y   = 0.21594;
	unit_cell.atom[17].z   = 0.95886;
	unit_cell.atom[17].mat = 2;
	unit_cell.atom[17].lc  = 17;
	unit_cell.atom[17].hc  = 7;
	unit_cell.atom[17].ni  = 12;
   unit_cell.atom[17].nm  = true;
	//-----------------------------
	unit_cell.atom[18].x   = 0.71972;
	unit_cell.atom[18].y   = 0.21594;
	unit_cell.atom[18].z   = 0.54114;
	unit_cell.atom[18].mat = 2;
	unit_cell.atom[18].lc  = 5;
	unit_cell.atom[18].hc  = 5;
	unit_cell.atom[18].ni  = 12;
   unit_cell.atom[18].nm  = true;
	//-----------------------------
	unit_cell.atom[19].x   = 0.28028;
	unit_cell.atom[19].y   = 0.78406;
	unit_cell.atom[19].z   = 0.04114;
	unit_cell.atom[19].mat = 2;
	unit_cell.atom[19].lc  = 19;
	unit_cell.atom[19].hc  = 1;
	unit_cell.atom[19].ni  = 12;
   unit_cell.atom[19].nm  = true;
	//-----------------------------
	

   unit_cell.cutoff_radius = 0.70856; // normalised to unit cell size

   uc::internal::calculate_interactions(unit_cell);

   // Set actual unit cell size after calculating interactions
   unit_cell.dimensions[0] *= unitcell::internal::unit_cell_size_x;
   unit_cell.dimensions[1] *= unitcell::internal::unit_cell_size_y;
   unit_cell.dimensions[2] *= unitcell::internal::unit_cell_size_z;

   return;

}

} // end of internal namespace
} // end of unitcell namespace
