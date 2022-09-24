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
 * @file main.h ******************** Documented \e "main.h" header ******************** \n
 *
 * @brief This is the header file of the main.c script.
 *
 * \b VERSION \b INFORMATION: Felipe Contreras, 2022-10-01, version 1.0.
 *
 * \b DESCRIPTION: This is the header file of the main.c script.
 *
 * \b PREREQUISITES: Always used.
 */

#ifndef __MAIN__
#define __MAIN__

//** >> Code Modules **/
#include "common.h"
#include "global_variables.h"
#include "garbage_collector.h"
#include "input.h"
#include "initialization.h"
#include "tree_construction.h"
#include "grid_density.h"
#include "potential.h"
#include "potential_head_node.h"
#include "poisson_error.h"
#include "grid_acceleration.h"
#include "particle_acceleration.h"
#include "timestep.h"
#include "particle_updating.h"
#include "observables.h"
#include "reset.h"
#include "tree_adaptation.h"


#include "output_main_parameters.h"
#include "output_snapshots.h"



#include "terminal_print.h"




#endif

