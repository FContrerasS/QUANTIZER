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

//** >> Local Functions
static void init_global_constants(void);
static void init_global_user_params(void);
static void init_global_ref_crit(void);
static void init_global_poisson_params(void);
static void init_global_force_params(void);
static void init_global_energies_params(void);
static void init_global_ptcl(void);
static void init_tree_head(void);
static void init_global_border_sim_box(void);
static void init_global_folder_params(void);
static void init_global_timer(void);
static void init_global_memory(void);
static void init_global_garbage_collector_parameters(void);
static void init_multigrid2_parameters(void);

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
int Gl_no_ptcl_initial;
int GL_no_ptcl;
vtype Maxdt ;
vtype meanmass;
vtype total_mass;
int fr_output;
int MaxIterations;
int no_grid_pow2;
int no_grid_pow3;
int boundary_type;

//** >> Refinement criteria parameters **/
vtype ref_criterion_mass;
int ref_criterion_ptcl;
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
vtype _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_2;
int check_poisson_error_method;
int multigrid_cycle;
int solverPreS;
int solverfinddphic;
int solverPostS;
int _NiterPreS_;
int _NiterPostS_; 
int _Niterfinddphic_;
int _Iter_branches_solver_;
vtype _w_SOR_;
vtype _w_SOR_HEAD_;
int head_pot_method;
int branch_pot_method;
int iter_between_check_potential_solution;
int branches_maximal_node_number_to_activate_conjugate_gradient;

//** >> Force parameters **/
int force_stencil;

//** >> Initializing energy parameters **/
int potential_energy_type;

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

//** >> Creating folders **/
bool folder_created;

//** >> Timer **/
clock_t GL_clock_begin;
double *GL_times;

//** >> MEMORY **/
double TOTAL_MEMORY_NODES;
double TOTAL_MEMORY_CELDAS;
double TOTAL_MEMORY_PARTICLES;
double TOTAL_MEMORY_CELL_STRUCT;
double TOTAL_MEMORY_CAJAS;
double TOTAL_MEMORY_GRID_POINTS;
double TOTAL_MEMORY_GRID_PROPERTIES;
double TOTAL_MEMORY_AUX;
double TOTAL_MEMORY_TENTACLES;
double TOTAL_MEMORY_STACK;

//** >> GARBAGE COLLECTOR **/
int Garbage_Collector_iter;

vtype **pp_phixx;
vtype **pp_rhoxx;
vtype **pp_restxx;
vtype **pp_dphixx;
vtype **zeros_xx;

static void
init_global_constants(void)
{
    //Constants
    _User_BoxSize_ = 5.0L; //kpc
    //_User_BoxSize_ = 0.1L; //kpc
    _PI_ = 3.14159265358979323846L;
    _Onesixth_ = 1.0L / 6.0L;
    _kpc_to_m_ = 3.08568e19L;
    _Msolar_to_kg_ = 1.98847e30L;
    tt = sqrt(_kpc_to_m_ * _kpc_to_m_ * _kpc_to_m_ * _User_BoxSize_ * _User_BoxSize_ * _User_BoxSize_ / (6.67430e-11L * _Msolar_to_kg_));
    _Mgyear_ = 3.1556952e13L / tt;
    _G_ = 1.0L;
    _G_ = 6.67430e-11 / (_kpc_to_m_ * _kpc_to_m_ * _kpc_to_m_ * _User_BoxSize_ * _User_BoxSize_ * _User_BoxSize_) * _Msolar_to_kg_ * tt * tt;
    printf("_User_BoxSize_ = %f\n", (double) _User_BoxSize_);
    printf("G = %.16f\n", (double)_G_);
    printf("_Mgyear_ = %1.6e\n", (double)_Mgyear_);
    printf("tt = %1.6e\n", (double)tt);
}

static void init_global_user_params(void)
{
    BoxSize = 1.0L;
    lmin = 5;     //Coarset level of refinement
    lmax = lmin + 3;  //Finest level of refinement
    no_lmin_cell = 1 << lmin; // Number of cells in the lmin level of refinement
    no_lmin_cell_pow2 = no_lmin_cell * no_lmin_cell;
    no_lmin_cell_pow3 = no_lmin_cell * no_lmin_cell * no_lmin_cell;
    no_grid = no_lmin_cell + 1;
    Gl_no_ptcl_initial = 10;
    GL_no_ptcl = Gl_no_ptcl_initial;
    //GL_no_ptcl = 7550; // 2995865; // 299586; // 231299 // 298159
    // GL_no_ptcl = 10000;
    Maxdt = 1000.0 * _Mgyear_;
    //meanmass = 100; //Currently only used on input.c
    // total_mass = GL_no_ptcl * meanmass;
    // total_mass = 0;
    fr_output = 12;
    MaxIterations = 100000000;
    no_grid_pow2 = no_grid * no_grid;
    no_grid_pow3 = no_grid * no_grid * no_grid;
    boundary_type = 2; // 0 = Periodic, 1 = reflexive, 2 = outflow
}

static void init_global_ref_crit(void)
{
    ref_criterion_mass = 1.0e100; // meanmass * 7;
    ref_criterion_ptcl = 1;
    n_exp = 1;   // n_exp = 0 is corrupted because particles can move between more than 1 level of refinement
    _CFL_ = 0.5; // CFL criteria 0.5
    _MAX_dt_ = _Mgyear_ * 1.0;
}

static void init_global_poisson_params(void)
{
    //** >> Poisson parameters **/
    // Relaxation solver at coarsest level
    /*
    check_poisson_error_method: {0,1,2}
    multigrid:   0 = V cycle, 1 = F cycle, 2 = W cycle.
    solver PresS PostS finddphic:      // 0 = Gauss-Saidel, 1 = Jacobi // 2 Conjugate Gradient
    _NiterPreS_: Solver iterations in Pre-Smoothing of the multigrid method
    _NiterPostS_: Solver iterations in Post-Smoothing of the multigrid method
    _Niterfinddphic_: Solver iterations on coarsest level of the multigrid method
    _Iter_branches_SOR_: Number of iterations in the potential computation in the
                         branches of the tree using Successive over-relaxation
    vtype _w_SOR_: The overrelaxation parameter
*/
    _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_ = 1000;
    _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ = (1.0e-8);
    _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_2 = (1.0e-8);
    check_poisson_error_method = 1;  //Only used Gauss-Said or Jacobi in multigrid
    multigrid_cycle = 0; 
    solverPreS = 0;
    solverPostS = 0;
    solverfinddphic = 0;
    _NiterPreS_ = 2;
    _NiterPostS_ = 2; 
    _Niterfinddphic_ = 2;
    _Iter_branches_solver_ = 25; // 25;
    _w_SOR_ = 1.9;
    _w_SOR_HEAD_ = 1.0;

    iter_between_check_potential_solution = 5;

    branches_maximal_node_number_to_activate_conjugate_gradient = 513; // 513, 216 = node with minimum size of 1 (+1 n_exp) size, (1 + 2*n_exp)^3 * 8

    head_pot_method = 0; // 0 = Multygrid, 1 = Conjugate gradient
    branch_pot_method = 1; // 0 = SOR, 1 = Conjugate gradient
}

static void init_global_force_params(void)
{
    force_stencil = 1;  // 0 = 3-points, 1 = 5-points
}

static void init_global_energies_params(void)
{
    potential_energy_type = 1; // 0 = Exact, 1 = approximation using potential grid
}

static void init_global_ptcl(void)
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

static void init_tree_head(void)
{
    GL_ptr_tree = (struct node *)malloc(sizeof(struct node));
}

static void init_global_border_sim_box(void)
{
    bder_os_sim = 1 > (n_exp-1) ? 1 : (n_exp-1); // Border outside the simulation box
    box_side_lmin = no_lmin_cell + 2 * bder_os_sim; // Side of the coarsest box
    box_side_lmin_pow2 = box_side_lmin * box_side_lmin;
    box_side_lmin_pow3 = box_side_lmin * box_side_lmin * box_side_lmin;
}

static void init_global_folder_params(void)
{
    folder_created = false;
}

static void init_global_timer(void)
{
    GL_times = (double *)calloc(100 , sizeof(double));
}

static void init_global_memory(void)
{
    TOTAL_MEMORY_NODES = 0;
    TOTAL_MEMORY_CELDAS = 0;
    TOTAL_MEMORY_PARTICLES = 0;
    TOTAL_MEMORY_CELL_STRUCT = 0;
    TOTAL_MEMORY_CAJAS = 0;
    TOTAL_MEMORY_GRID_POINTS = 0;
    TOTAL_MEMORY_GRID_PROPERTIES = 0;
    TOTAL_MEMORY_AUX = 0;
    TOTAL_MEMORY_TENTACLES = 0;
    TOTAL_MEMORY_STACK = 0;
}

static void init_global_garbage_collector_parameters(void)
{
    Garbage_Collector_iter = 10000000; // Number of time-steps between each garbage collector
}

static void init_multigrid2_parameters(void)
{
    int size, size_pow3;
    pp_phixx = (vtype **)malloc((lmin - 1) * sizeof(vtype *));
    pp_rhoxx = (vtype **)malloc((lmin - 1) * sizeof(vtype *));
    pp_restxx = (vtype **)malloc((lmin - 2) * sizeof(vtype *));
    pp_dphixx = (vtype **)malloc((lmin - 2) * sizeof(vtype *));
    zeros_xx = (vtype **) malloc((lmin - 1) * sizeof(vtype *));
    size = ((1 << (lmin )) + 1);
    size_pow3 = size * size * size;
    for (int i = 0; i < lmin - 2;i++)
    {
        pp_restxx[i] = (vtype *)calloc(size_pow3, sizeof(vtype));
        pp_dphixx[i] = (vtype *)calloc(size_pow3, sizeof(vtype));
        zeros_xx[i] = (vtype *)calloc(size_pow3, sizeof(vtype));
        size = ((1 << (lmin - i - 1)) + 1);
        size_pow3 = size * size * size;
        pp_phixx[i+1] = (vtype *)malloc(size_pow3 * sizeof(vtype));
        pp_rhoxx[i+1] = (vtype *)calloc(size_pow3, sizeof(vtype)); 
    }

    zeros_xx[lmin-2] = (vtype *)calloc(size_pow3, sizeof(vtype));
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

    //** >> Initializing Force parameters **/
    init_global_force_params();

    //** >> Initializing energy parameters **/
    init_global_energies_params();

    //** >> Initializing particles **/
    init_global_ptcl();

    //** >> Initializing tree head **/ 
    init_tree_head();

    //** >> Initializing border of simulation box parameters **/ 
    init_global_border_sim_box();
    
    //** >> Initializing output folder parameters **/ 
    init_global_folder_params();

    //** >> Initializing timer **/ 
    init_global_timer();

    //** >> Initializing memory parameters **/ 
    init_global_memory();

    //** >> Garbage Collector **/
    init_global_garbage_collector_parameters();

    //** >> multigrid2 parameters **/
    init_multigrid2_parameters();
}
