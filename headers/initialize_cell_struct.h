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

/**
 * @file initialize_cell_struct.h ******************** Documented \e "initialize_cell_struct.h" header ******************** \n
 *
 * @brief This is the header file of the initialize_cell_struct.c script.
 *
 * **VERSION INFORMATION**: Felipe Contreras, 2022-10-01, version 1.0.
 *
 * **DESCRIPTION**: This is the header file of the initialize_cell_struct.c script.
 *
 * **PREREQUISITES**: Always used.
 */

#ifndef __INITIALIZECELLSTRUCT__
#define __INITIALIZECELLSTRUCT__

#include "common.h"

/**
 * @brief Used to store the mass of the cell and the particles inside of it.
 *
 * **SHORT DESCRIPTION**: Used to store the mass of the cell and the particles
 * inside of it.
 *
 * **LONG DESCRIPTION**:
 *
 * Used to store the mass of the cell and the particles inside of it.
 *
 * **ILLUSTRATIVE EXAMPLES**:
 * - [a]  Trivial.
 *
 * **RATIONALES**:
 * - [a]  a
 *
 * **NOTES**:
 * - [a]  This structure can be modified in the future to include other useful 
 * parameters, such as the density or the velocity of the cell.
 */

struct cell_struct
{
    int *ptr_ptcl;      /**< Global Particle Index of the particles in the cell */ 
    int ptcl_cap; /**< *Capacity* (\ref node "see key Concepts [a]") of the particles in the cell array ::ptr_ptcl */
    int ptcl_size; /**< *Size* (\ref node "see key Concepts [b]") of the particles in the cell array ::ptr_ptcl */ 
    vtype cell_mass;    /**< Total mass of the particles inside of the cell */
};

void initialize_cell_struct(struct cell_struct *ptr_cell_struct);

#endif
