/*
 * common.h 
 *
 * Header file with the most commen parameters declared
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

#ifndef __COMMON__
#define __COMMON__

//** >> Standard Libraries **/
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "string.h"
#include <sys/stat.h>
#include <stdbool.h>
#include <stdarg.h>

#define _VTYPE_ 2
#if _VTYPE_ == 1
typedef float vtype;
#define myabs fabsf
#elif _VTYPE_ == 2
typedef double vtype;
#define myabs fabs
#elif _VTYPE_ == 3
typedef long float vtype;
#define myabs fabsfl
#else
typedef double vtype;
#define myabs fabs
#endif




//** Tree structure **/
#include "initialize_node.h"

#define _SUCCESS_ 0 /* integer returned after successful call of a function */
#define _FAILURE_ 1 /* integer returnd after failure in a function */

#define _TRUE_ 1 /* integer associated to true statement */
#define _FALSE_ 0 /* integer associated to false statement */

//** >> Terminal Colors **/
#define KNRM "\x1B[0m" // Normal
#define KRED "\x1B[31m" // Red
#define KGRN "\x1B[32m" // Green
#define KYEL "\x1B[33m" // Yellow
#define KBLU "\x1B[34m" // Blue
#define KMAG "\x1B[35m" // Magenta
#define KCYN "\x1B[36m" // Cyan 
#define KWHT "\x1B[37m" // White

//Abbreviation

//Constants
extern vtype _User_BoxSize_; //kpc
extern vtype _PI_;
extern vtype _Onesixth_;
extern vtype _kpc_to_m_;
extern vtype _Msolar_to_kg_;
extern vtype tt;
extern vtype _Mgyear_;
extern vtype _G_;
 

// Initial Parameters
extern vtype BoxSize;
extern int lmin;     //Coarset level of refinement
extern int lmax;  //Finest level of refinement
extern int no_lmin_cell; // Number of cells in the lmin level of refinement
extern int no_lmin_cell_pow2;
extern int no_lmin_cell_pow3;
extern int no_grid;
extern int GL_no_ptcl;
extern vtype Maxdt ;
extern vtype meanmass;
extern vtype total_mass;
extern int fr_output;
extern int MaxIterations;
extern int no_grid_pow2;
extern int no_grid_pow3;

//** >> Refinement criteria parameters **/
extern vtype ref_criterion_mass;
extern int n_exp;
extern vtype _CFL_; // CFL criteria 0.5
extern vtype _MAX_dt_;

//** >> Poisson parameters **/
// Relaxation solver at coarsest level
extern int _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_;
extern vtype _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_;
extern int check_poisson_error_method; // 0,1,2 = Method used
extern int multigrid; // 0 = V cycle, 1 = F cycle, 2 = W cycle.
extern int solver; // 0 = Gauss-Saidel, 1 = Jacobi
extern int _NiterPreS_; // _NiterPreS_:      Number of Solver iterations used for a Pre-Smoothing on multigrid method
extern int _NiterPostS_; //  _NiterPostS_:     Number of Solver iterations used for a Post-Smoothing on multigrid method
extern int _Niterfinddphic_; // Number of Solver iterations used in the coarsest level on the multigrid method
extern int _Iter_branches_SOR_; // Number of iterations in the potential computation in the branches of the tree using Successive over-relaxation
extern vtype _w_SOR_;          // The overrelaxation parameter

//** >> Particles **/
extern vtype *GL_ptcl_mass;
extern vtype *GL_ptcl_x;
extern vtype *GL_ptcl_y;
extern vtype *GL_ptcl_z;
extern vtype *GL_ptcl_vx;
extern vtype *GL_ptcl_vy;
extern vtype *GL_ptcl_vz;
extern vtype *GL_ptcl_ax;
extern vtype *GL_ptcl_ay;
extern vtype *GL_ptcl_az;
extern vtype **GL_ptcl;
extern bool *GL_ptcl_updating_flag;


//** >> Head or Main node **/
extern struct node *GL_ptr_tree;

//** >> Tentacles pointer **/
extern struct node ***GL_tentacles_old; // Array of arrays of pointers. Organized first by levels, then by pointers
extern struct node ***GL_tentacles_new;
extern int *GL_tentacles_cap;   // Capacity of pointers in each level
extern int *GL_tentacles_size; // Number of pointers in each level
extern int GL_tentacles_level_max; // Maximum level of refinement of the tentacles

//** >> Outputs **/
extern char file_data_name[1000];
extern char folder_name[1000];
extern bool folder_created;


//** >> Timer **/
extern clock_t GL_clock_begin;
extern double *GL_times;

//** >> Border of the simulation box **/
extern int bder_os_sim;        // Border outside the simulation box
extern int box_side_lmin; // Side of the coarsest box
extern int box_side_lmin_pow2;
extern int box_side_lmin_pow3;

//** >> OBSERVABLES **/
extern vtype _Min_Particle_Distance_For_Energy_Computation_;


// MEMORY
extern double TOTAL_MEMORY_NODES;
extern double TOTAL_MEMORY_CELDAS;
extern double TOTAL_MEMORY_PARTICULAS;
extern double TOTAL_MEMORY_CAJAS;
extern double TOTAL_MEMORY_AUX;
extern double TOTAL_MEMORY_TENTACLES;
extern double TOTAL_MEMORY_OTROS;

#endif