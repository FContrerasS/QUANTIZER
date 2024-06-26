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
 * @file common.h
 *
 * \f[{\color{magenta} \mathbf{ DOCUMENTED\ ``common.h"\ HEADER}}\f]
 *
 * @brief Header file with multiple commen paramter and global variables
 *
 * **VERSION INFORMATION**: Felipe Contreras, 2022-10-01, version 1.0.
 *
 * **DESCRIPTION**: This is the header file with multiple commen paramter and
 * global variables.
 *
 * **PREREQUISITES**: Always used.
 */

#ifndef __COMMON__
#define __COMMON__

//* >> Standard Libraries *//
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "string.h"
#include <sys/stat.h>
#include <stdbool.h>
#include <stdarg.h>
#include <omp.h>

#define _VTYPE_ 2 /**< Defining type of data: 1 = float; 2 = double; 3 = long double. other number = double */

#if _VTYPE_ == 1
typedef float vtype; /**< Type of date \f${\color{red} \mathbf{float}}\f$ assigned by the user in the macros \ref _VTYPE_ */
#define myabs fabsf 
#define myfloor floorf
#define mysqrt sqrtf
#define mycbrt cbrtf
#elif _VTYPE_ == 2
typedef double vtype; /**< Type of date \f${\color{red} \mathbf{double}}\f$ assigned by the user in the macros \ref _VTYPE_ */
#define myabs fabs
#define myfloor floor
#define mysqrt sqrt
#define mycbrt cbrt
#elif _VTYPE_ == 3
typedef long double vtype; /**< Type of date \f${\color{red} \mathbf{long double}}\f$ assigned by the user in the macros \ref _VTYPE_ */
#define myabs fabsl
#define myfloor floorl
#define mysqrt sqrtl
#define mycbrt cbrtl
#else
typedef double vtype; /**< Type of date \f${\color{red} \mathbf{double}}\f$ assigned by the user in the macros \ref _VTYPE_ */
#define myabs fabs
#define myfloor floor
#define mysqrt sqrt
#define mycbrt cbrt
#endif






#define _SUCCESS_ 0 /* integer returned after successful call of a function */
#define _FAILURE_ 1 /* integer returnd after failure in a function */

#define _TRUE_ 1  /* integer associated to true statement */
#define _FALSE_ 0 /* integer associated to false statement */

//* Tree structure *//
#include "initialize_cell_struct.h"
#include "initialize_node.h"
#include "ptcl_idx_to_box_idx.h"

//* >> Terminal Colors *//
#define KNRM "\x1B[0m" // Normal
#define KRED "\x1B[31m" // Red
#define KGRN "\x1B[32m" // Green
#define KYEL "\x1B[33m" // Yellow
#define KBLU "\x1B[34m" // Blue
#define KMAG "\x1B[35m" // Magenta
#define KCYN "\x1B[36m" // Cyan 
#define KWHT "\x1B[37m" // White

// Abbreviation

// Constants
extern vtype _User_BoxSize_; /**< Global Variable Constants: User boxsize in kpc unit.*/
extern vtype _PI_;           /**< Global Variable Constants: Pi value. */
extern vtype _Onesixth_;     /**< Global Variable Constants: The value of \f$1\over 6\f$. */
extern vtype _kpc_to_m_;     /**< Global Variable Constants: The dimensionless \f$x\f$ value in the relation \f$1\textrm{ kpc} = x \cdot 1\textrm{ m}\f$. */
extern vtype _conversion_dist_;
extern vtype _Msolar_to_kg_; /**< Global Variable Constants: The dimensionless \f$x\f$ value in the relation \f$1 \textup{ M}_\odot = x \cdot 1\textrm{ kg}\f$. */
extern vtype _year_to_s_;
extern vtype _conversion_time_;
extern vtype _conversion_velocity_;
//extern vtype tt;             /**< Global Variable Constants: The unit of time in code units, such that #_G_ = 1 in these code units. */
extern vtype _Myear_;       /**< Global Variable Constants: The value of one Mega year in seconds */
extern vtype _Gyear_;       /**< Global Variable Constants: The value of one Giga year in seconds */
extern vtype _G_;            /**< Global Variable Constants: The gravitational constant *G* in units of \f${\textrm{m}^3\over \textrm{kg s}^2} \f$ */

// Initial Parameters
extern int NUMBER_OF_THEADS;
extern vtype BoxSize;
extern int lmin;     //Coarset level of refinement
extern int lmax;  //Finest level of refinement
extern int no_lmin_cell; // Number of cells in the lmin level of refinement
extern int no_lmin_cell_pow2;
extern int no_lmin_cell_pow3;
extern int no_grid;
extern int GL_no_ptcl_initial;
extern int GL_no_ptcl_final;
extern vtype Maxdt ;
extern vtype meanmass;
extern vtype GL_total_mass_initial;
extern vtype GL_total_mass_final;
extern int fr_output;
extern int MaxIterations;
extern int no_grid_pow2;
extern int no_grid_pow3;
extern int min_box_extra_size_per_side_user;
extern int min_box_extra_size_per_side_real;

//* >> Refinement criteria parameters *//
extern vtype ref_criterion_mass;
extern int ref_criterion_ptcl;
extern int n_exp; /**< Global Variable Refinement Criteria: The n_exp parameter tells us how many neighboring cells should be marked for refinement regardless of whether or not they satisfy the refinement criterion (same as in the *RAMSES* code). */

extern vtype _CFL_; // CFL criteria 0.5
extern vtype _MAX_dt_;

//* >> Initial Center of Mass *//
//extern vtype GL_cm[3]; // Center of mass
extern vtype GL_cm_x; // Center of momentum pos x
extern vtype GL_cm_y; // Center of momentum pos y
extern vtype GL_cm_z; // Center of momentum pos z
extern vtype GL_cm_vx; // Center of momentum vel x
extern vtype GL_cm_vy; // Center of momentum vel y
extern vtype GL_cm_vz; // Center of momentum vel z

//* >> Poisson parameters *//
// Relaxation solver at coarsest level
extern int _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_;
extern vtype _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_;
extern vtype _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_2;
extern int check_poisson_error_method; // 0,1,2 = Method used
extern int multigrid_cycle;            // 0 = V cycle, 1 = F cycle, 2 = W cycle.
extern int solverPreS; // 0 = Gauss-Saidel, 1 = Jacobi // 2 Conjugate Gradient
extern int solverfinddphic;
extern int solverPostS;
extern int _NiterPreS_; // _NiterPreS_:      Number of Solver iterations used for a Pre-Smoothing on multigrid method
extern int _NiterPostS_; //  _NiterPostS_:     Number of Solver iterations used for a Post-Smoothing on multigrid method
extern int _Niterfinddphic_; // Number of Solver iterations used in the coarsest level on the multigrid method
extern int _Iter_branches_solver_; // Number of iterations in the potential computation in the branches of the tree using Successive over-relaxation
extern vtype _w_SOR_;          // The overrelaxation parameter
extern vtype _w_SOR_HEAD_;
extern int head_pot_method;
extern int branch_pot_method;
extern int iter_between_check_potential_solution;
extern int branches_maximal_node_number_to_activate_conjugate_gradient;

//* >> Time Step paramters *//
extern int time_step_method;

//* >> Force parameters *//
extern int force_stencil;

//* >> Initializing energy parameters *//
extern bool compute_observables_FLAG;
extern int potential_energy_type;
extern vtype *GL_energies;
extern vtype *GL_momentum;

//* >> Particles *//
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
extern int *GL_ptcl_ID;

//* >> Head or Main node *//
extern struct node *GL_ptr_tree;

//* >> Tentacles pointer *//
extern struct node ***GL_tentacles; // Array of arrays of pointers. Organized first by levels, then by pointers
extern int *GL_tentacles_cap;       /**< Global Variable Tentacles: Array of *Capacities* (see Key Concepts \ref Key_Concepts_Capacity "Capacity") at every level of the \ref tree_adaptation__REMINDER__Tentacles "Tentacles" */
extern int *GL_tentacles_size;  /**< Global Variable Tentacles: Array of *Sizes* (see Key Concepts \ref Key_Concepts_Size "Size") at every level of the \ref tree_adaptation__REMINDER__Tentacles "Tentacles" */
extern int GL_tentacles_level_max; /**< Maximum level of refinement of the \ref tree_adaptation__REMINDER__Tentacles "Tentacles" */

//* >> Pool of nodes *//
extern struct node *GL_pool_node_start;
extern struct node *GL_pool_node_end;

//* >> Outputs *//
extern char file_data_name[1000];
extern char folder_name[1000];
extern bool folder_created;


//* >> Timer *//
//extern clock_t GL_clock_begin;
extern struct timespec GL_start, GL_finish;
extern double *GL_times;

//* >> Border of the simulation box *//
extern int bder_os_sim;        // Border outside the simulation box
extern int box_side_lmin; // Side of the coarsest box
extern int box_side_lmin_pow2;
extern int box_side_lmin_pow3;

// MEMORY
extern double TOTAL_MEMORY_NODES;
extern double TOTAL_MEMORY_CELDAS;
extern double TOTAL_MEMORY_PARTICLES;
extern double TOTAL_MEMORY_CELL_STRUCT;
extern double TOTAL_MEMORY_CAJAS;
extern double TOTAL_MEMORY_GRID_POINTS;
extern double TOTAL_MEMORY_GRID_PROPERTIES;
extern double TOTAL_MEMORY_AUX;
extern double TOTAL_MEMORY_TENTACLES;
extern double TOTAL_MEMORY_OTROS;
extern double TOTAL_MEMORY_STACK;

// GARBAGE COLLECTOR
extern int Garbage_Collector_iter;

extern vtype **pp_phixx;
extern vtype **pp_rhoxx;
extern vtype **pp_restxx;
extern vtype **pp_dphixx;
extern vtype **zeros_xx;

//* >> Boundary parameters *//
extern int bdry_cond_type;
extern int GL_max_lv_ref_crosses_the_entire_simulation; // Maximum level of refinement such that it contains a zone of refinement which crosses the entire simulation

#endif
