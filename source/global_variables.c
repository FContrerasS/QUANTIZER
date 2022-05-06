/*
 * global_variables.c
 *
 * Define global variables
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

#include "global_variables.h"

    // Constants
    vtype _User_BoxSize_; //** >> User box size **/
vtype _PI_;
vtype _Onesixth_;
vtype _kpc_to_m_;
vtype _Msolar_to_kg_;
vtype tt;
vtype _Mgyear_;
vtype _G_;

// Initial Parameters
vtype BoxSize;
int lmin;     //Coarset level of refinement
int lmax;  //Finest level of refinement
int no_lmin_cell; // Number of cells in the lmin level of refinement
int no_lmin_cell_pow2;
int no_lmin_cell_pow3;
int no_grid;
int GL_no_ptcl;
vtype Maxdt ;
vtype meanmass;
vtype total_mass;
int fr_output;
int MaxIterations;
int no_grid_pow2;
int no_grid_pow3;

//** >> Refinement criteria parameters **/
vtype ref_criterion_mass;
int n_exp;
vtype _CFL_; 
vtype _MAX_dt_;

//** >> Poisson parameters **/
// Relaxation solver at coarsest level
/*
check_poisson_error_method: {0,1,2}
multigrid:      0 = V cycle, 1 = F cycle, 2 = W cycle.
solver:      0 = Gauss-Saidel, 1 = Jacobi
_NiterPreS_: Solver iterations in Pre-Smoothing of the multigrid method
_NiterPostS_: Solver iterations in Post-Smoothing of the multigrid method
_Niterfinddphic_: Solver iterations on coarsest level of the multigrid method
_Iter_branches_SOR_: Number of iterations in the potential computation in the
                     branches of the tree using Successive over-relaxation
vtype _w_SOR_: The overrelaxation parameter
*/
int _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_;
vtype _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_;
int check_poisson_error_method;
int multigrid;
int solver;
int _NiterPreS_;
int _NiterPostS_; 
int _Niterfinddphic_; 
int _Iter_branches_SOR_;
vtype _w_SOR_;         


//** >> Defining Particles Parameters **/
vtype *GL_ptcl_mass; // Mass 
vtype *GL_ptcl_x;    // X position
vtype *GL_ptcl_y;    // Y position
vtype *GL_ptcl_z;    // Z position
vtype *GL_ptcl_vx;   // X velocity
vtype *GL_ptcl_vy;   // Y velocity
vtype *GL_ptcl_vz;   // Z velocity
vtype *GL_ptcl_ax; // Accelerations
vtype *GL_ptcl_ay;
vtype *GL_ptcl_az;
vtype **GL_ptcl;     // All 10 ptcl parameters
bool *GL_ptcl_updating_flag;  // Particle updating state

//** >> Head or Main node **/
struct node *GL_ptr_tree; // Pointer to the tree head or  coarsest level

//** >> Tentacles struct pointer **/
struct node ***GL_tentacles;
int *GL_tentacles_cap;      // Capacity of pointers in each level
int *GL_tentacles_size;     // Number of pointers in each level
int GL_tentacles_level_max; // Maximum level of refinement of the tentacles

//** >> Initializing Pool node list **/
struct node *GL_pool_node_start = NULL;
struct node *GL_pool_node_end = NULL;

//** >> Border of the simulation box **/
int bder_os_sim;        // Border outside the simulation box
int box_side_lmin; // Side of the coarsest box
int box_side_lmin_pow2;
int box_side_lmin_pow3;


//** >> OBSERVABLES **/
vtype _Min_Particle_Distance_For_Energy_Computation_;

//** >> Creating folders **/
bool folder_created;

//** >> Timer **/
clock_t GL_clock_begin;
double *GL_times;

//** >> MEMORY **/
double TOTAL_MEMORY_NODES;
double TOTAL_MEMORY_CELDAS;
double TOTAL_MEMORY_PARTICULAS;
double TOTAL_MEMORY_CAJAS;
double TOTAL_MEMORY_GRID_POINTS;
double TOTAL_MEMORY_GRID_PROPERTIES;
double TOTAL_MEMORY_AUX;
double TOTAL_MEMORY_TENTACLES;
double TOTAL_MEMORY_STACK;

static void init_global_constants()
{
    //Constants
    _User_BoxSize_ = 0.1L; //kpc 
    _PI_ = 3.14159265358979323846L;
    _Onesixth_ = 1.0L / 6.0L;
    _kpc_to_m_ = 3.08568e19L;
    _Msolar_to_kg_ = 1.98847e30L;
    tt = sqrt(_kpc_to_m_ * _kpc_to_m_ * _kpc_to_m_ * _User_BoxSize_ * _User_BoxSize_ * _User_BoxSize_ / (6.67430e-11L * _Msolar_to_kg_));
    _Mgyear_ = 3.1556952e13L / tt;
    _G_ = 1.0L;
}

static void init_global_user_params()
{
    BoxSize = 1.0L;
    lmin = 5;     //Coarset level of refinement
    lmax = lmin + 4;  //Finest level of refinement
    no_lmin_cell = 1 << lmin; // Number of cells in the lmin level of refinement
    no_lmin_cell_pow2 = no_lmin_cell * no_lmin_cell;
    no_lmin_cell_pow3 = no_lmin_cell * no_lmin_cell * no_lmin_cell;
    no_grid = no_lmin_cell + 1;
    GL_no_ptcl = 10000;
    Maxdt = 3.0 * _Mgyear_;
    meanmass = 100;
    total_mass = GL_no_ptcl * meanmass;
    fr_output = 40;
    MaxIterations = 1000000;
    no_grid_pow2 = no_grid * no_grid;
    no_grid_pow3 = no_grid * no_grid * no_grid;

}

static void init_global_ref_crit()
{
    ref_criterion_mass = meanmass * 2;
    n_exp = 1;  //n_exp = 0 is corrupted because particles can move between more than 1 level of refinement
    _CFL_ = 0.5; // CFL criteria 0.5
    _MAX_dt_ = 6.7e-6;
}

static void init_global_poisson_params()
{
    //** >> Poisson parameters **/
    // Relaxation solver at coarsest level
    /*
    check_poisson_error_method: {0,1,2}
    multigrid:   0 = V cycle, 1 = F cycle, 2 = W cycle.
    solver:      0 = Gauss-Saidel, 1 = Jacobi
    _NiterPreS_: Solver iterations in Pre-Smoothing of the multigrid method
    _NiterPostS_: Solver iterations in Post-Smoothing of the multigrid method
    _Niterfinddphic_: Solver iterations on coarsest level of the multigrid method
    _Iter_branches_SOR_: Number of iterations in the potential computation in the
                         branches of the tree using Successive over-relaxation
    vtype _w_SOR_: The overrelaxation parameter
*/
    _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_ = 100;
    _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ = (1.5e-10);
    check_poisson_error_method = 0; 
    multigrid = 0; 
    solver = 0; 
    _NiterPreS_ = 2; 
    _NiterPostS_ = 2; 
    _Niterfinddphic_ = 2; 
    _Iter_branches_SOR_ = 20; 
    _w_SOR_ = 1.9;         

}

static void init_global_ptcl()
{
    GL_ptcl_mass = (vtype *)calloc(GL_no_ptcl , sizeof(vtype));
    GL_ptcl_x = (vtype *)calloc(GL_no_ptcl , sizeof(vtype));
    GL_ptcl_y = (vtype *)calloc(GL_no_ptcl , sizeof(vtype));
    GL_ptcl_z = (vtype *)calloc(GL_no_ptcl , sizeof(vtype));
    GL_ptcl_vx = (vtype *)calloc(GL_no_ptcl , sizeof(vtype));
    GL_ptcl_vy = (vtype *)calloc(GL_no_ptcl , sizeof(vtype));
    GL_ptcl_vz = (vtype *)calloc(GL_no_ptcl , sizeof(vtype));
    GL_ptcl_ax = (vtype *)calloc(GL_no_ptcl , sizeof(vtype));
    GL_ptcl_ay = (vtype *)calloc(GL_no_ptcl , sizeof(vtype));
    GL_ptcl_az = (vtype *)calloc(GL_no_ptcl , sizeof(vtype));

    GL_ptcl = (vtype **)malloc(10 * sizeof(vtype *));
    GL_ptcl[0] = GL_ptcl_mass;
    GL_ptcl[1] = GL_ptcl_x;
    GL_ptcl[2] = GL_ptcl_y;
    GL_ptcl[3] = GL_ptcl_z;
    GL_ptcl[4] = GL_ptcl_vx;
    GL_ptcl[5] = GL_ptcl_vy;
    GL_ptcl[6] = GL_ptcl_vz;
    GL_ptcl[7] = GL_ptcl_ax;
    GL_ptcl[8] = GL_ptcl_ay;
    GL_ptcl[9] = GL_ptcl_az;
}

static void init_tree_head()
{
    GL_ptr_tree = (struct node *)malloc(sizeof(struct node));
}

static void init_global_border_sim_box()
{
    bder_os_sim = 1 > (n_exp-1) ? 1 : (n_exp-1); // Border outside the simulation box
    box_side_lmin = no_lmin_cell + 2 * bder_os_sim; // Side of the coarsest box
    box_side_lmin_pow2 = box_side_lmin * box_side_lmin;
    box_side_lmin_pow3 = box_side_lmin * box_side_lmin * box_side_lmin;
}

static void init_global_observables()
{
    _Min_Particle_Distance_For_Energy_Computation_ = 1.0e-12L;

}

static void init_global_folder_params()
{
    folder_created = false;
}

static void init_global_timer()
{
    GL_times = (double *)calloc(100 , sizeof(double));
}

static void init_global_memory()
{
    TOTAL_MEMORY_NODES = 0;
    TOTAL_MEMORY_CELDAS = 0;
    TOTAL_MEMORY_PARTICULAS = 0;
    TOTAL_MEMORY_CAJAS = 0;
    TOTAL_MEMORY_GRID_POINTS = 0;
    TOTAL_MEMORY_GRID_PROPERTIES = 0;
    TOTAL_MEMORY_AUX = 0;
    TOTAL_MEMORY_TENTACLES = 0;
    TOTAL_MEMORY_STACK = 0;
}

void global_variables()
{
    //** >> Initializing constants **/ 
    init_global_constants();

    //** >> Initializing user parameters **/ 
    init_global_user_params();

    //** >> Initializing refinement criteria parameters **/
    init_global_ref_crit();

    //** >> Initializing Poisson parameters **/ 
    init_global_poisson_params();

    //** >> Initializing particles **/ 
    init_global_ptcl();

    //** >> Initializing tree head **/ 
    init_tree_head();

    //** >> Initializing border of simulation box parameters **/ 
    init_global_border_sim_box();
    
    //** >> Initializing observables parameters **/ 
    init_global_observables();

    //** >> Initializing output folder parameters **/ 
    init_global_folder_params();

    //** >> Initializing timer **/ 
    init_global_timer();

    //** >> Initializing memory parameters **/ 
    init_global_memory();
}
