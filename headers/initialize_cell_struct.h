/*
 * initialize_cell_structure.h
 *
 * Header file with the most initialize_cell_structure parameters declared
 *
 * Felipe Contreras
 * felipe.contrerass@postgrado.uv.cl
 */

/*
 * Copyright(c) 2022 Felipe Contreras
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __INITIALIZECELLSTRUCT__
#define __INITIALIZECELLSTRUCT__

#include "common.h"

struct cell_struct
{
    int *ptr_ptcl;  // Indexes of the Particles in the cell
    int ptcl_cap;   //Maximum capacity in the number of particles in the cell
    int ptcl_size;  //Number of particles in the cell
    vtype cell_mass;    // Cell mass
};

void initialize_cell_struct(struct cell_struct *ptr_cell_struct);

#endif
