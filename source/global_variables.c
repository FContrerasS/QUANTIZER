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

//* >> Local Functions
static void init_global_constants(void);
static void init_global_user_params(void);
static void init_global_ref_crit(void);
static void init_global_poisson_params(void);
static void init_global_time_step_params(void);
static void init_global_force_params(void);
static void init_global_energies_params(void);
static void init_global_ptcl(void);
static void init_tree_head(void);
static void init_global_border_sim_box(void);
static void init_global_folder_params(void);
static void init_global_timer(void);
static void init_global_memory(void);
static void init_global_garbage_collector_parameters(void);
static void init_global_boundary_parameters(void);
static void init_multigrid2_parameters(void);

// Constants
vtype _User_BoxSize_; //* >> User box size *//
vtype _PI_;
vtype _Onesixth_;
vtype _kpc_to_m_;
vtype _conversion_dist_;
vtype _Msolar_to_kg_;
vtype _year_to_s_;
vtype _conversion_time_;
vtype _conversion_velocity_;
//vtype tt;
vtype _Myear_;
vtype _Gyear_;
vtype _G_;

// Initial Parameters
int NUMBER_OF_THEADS;
vtype BoxSize;
int lmin;         // Coarset level of refinement
int lmax;         // Finest level of refinement
int no_lmin_cell; // Number of cells in the lmin level of refinement
int no_lmin_cell_pow2;
int no_lmin_cell_pow3;
int no_grid;
int GL_no_ptcl_initial;
int GL_no_ptcl_final;
vtype Maxdt;
vtype meanmass;
vtype GL_total_mass_initial;
vtype GL_total_mass_final;
int fr_output;
int MaxIterations;
int no_grid_pow2;
int no_grid_pow3;

int min_box_extra_size_per_side_user;
int min_box_extra_size_per_side_real;


//* >> Refinement criteria parameters *//
vtype ref_criterion_mass;
int ref_criterion_ptcl;
int n_exp;
vtype _CFL_;
vtype _MAX_dt_;

//* >> Initial Center of Mass *//
//vtype GL_cm[3]; // Center of mass
vtype GL_cm_x; // Center of momentum pos x
vtype GL_cm_y; // Center of momentum pos y
vtype GL_cm_z; // Center of momentum pos z
vtype GL_cm_vx; // Center of momentum vel x
vtype GL_cm_vy; // Center of momentum vel y
vtype GL_cm_vz; // Center of momentum vel z

//* >> Poisson parameters *//
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

//* >> Time Step paramters *//
int time_step_method;

//* >> Force parameters *//
int force_stencil;

//* >> Initializing energy parameters *//
bool compute_observables_FLAG;
int potential_energy_type;
vtype *GL_energies;
vtype *GL_momentum;

//* >> Defining Particles Parameters *//
vtype *GL_ptcl_mass; // Mass
vtype *GL_ptcl_x;    // X position
vtype *GL_ptcl_y;    // Y position
vtype *GL_ptcl_z;    // Z position
vtype *GL_ptcl_vx;   // X velocity
vtype *GL_ptcl_vy;   // Y velocity
vtype *GL_ptcl_vz;   // Z velocity
vtype *GL_ptcl_ax;   // Accelerations
vtype *GL_ptcl_ay;
vtype *GL_ptcl_az;
vtype **GL_ptcl;             // All 10 ptcl parameters
bool *GL_ptcl_updating_flag; // Particle updating state
int *GL_ptcl_ID;

//* >> Head or Main node *//
struct node *GL_ptr_tree; // Pointer to the tree head or  coarsest level

//* >> Tentacles struct pointer *//
struct node ***GL_tentacles;
int *GL_tentacles_cap;      // Capacity of pointers in each level
int *GL_tentacles_size;     // Number of pointers in each level
int GL_tentacles_level_max; // Maximum level of refinement of the tentacles

//* >> Initializing Pool node list *//
struct node *GL_pool_node_start = NULL;
struct node *GL_pool_node_end = NULL;

//* >> Border of the simulation box *//
int bder_os_sim;   // Border outside the simulation box
int box_side_lmin; // Side of the coarsest box
int box_side_lmin_pow2;
int box_side_lmin_pow3;

//* >> Creating folders *//
bool folder_created;

//* >> Timer *//
//clock_t GL_clock_begin;
struct timespec GL_start, GL_finish;
double *GL_times;

//* >> MEMORY *//
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

//* >> GARBAGE COLLECTOR *//
int Garbage_Collector_iter;

vtype **pp_phixx;
vtype **pp_rhoxx;
vtype **pp_restxx;
vtype **pp_dphixx;
vtype **zeros_xx;

//* >> Boundary parameters *//
int bdry_cond_type;
int GL_max_lv_ref_crosses_the_entire_simulation; // Maximum level of refinement such that it contains a zone of refinement which crosses the entire simulation

static void
init_global_constants(void)
{
  // Constants
  //_User_BoxSize_ = 10.0 * 4.8476302883992e-9; // kpc;;;;     1 au = 4.8476302883992e-9 kpc
  //_User_BoxSize_ = 110000.0; //40 Mpc for galaxy merger
  _User_BoxSize_ = 4000.0; //4 Mpc Plummer model
  _PI_ = 3.14159265358979323846;
  _Onesixth_ = 1.0 / 6.0;
  _kpc_to_m_ = 3.08567758128e19; 
  _conversion_dist_ = 1.0/_User_BoxSize_; // x [kpc] = x * _conversion_dist_ [\hat{kpc}], where, [\hat{kpc}] are the code distance units
  _Msolar_to_kg_ = 1.98847e30;
  _year_to_s_ = 31558149.756; // Year (how long it takes for the Earth to fully orbit the Sun) to second
  _conversion_time_ = mysqrt(6.67430e-11 / (_kpc_to_m_ * _kpc_to_m_ * _kpc_to_m_) * (_conversion_dist_*_conversion_dist_*_conversion_dist_)
       * _Msolar_to_kg_ * _year_to_s_ * _year_to_s_); //x [year] = x * _conversion_time_ [\hat{year}], where, [\hat{year}]  are the code time units
  //tt = mysqrt(_kpc_to_m_ * _kpc_to_m_ * _kpc_to_m_ * _User_BoxSize_ * _User_BoxSize_ * _User_BoxSize_ / (6.67430e-11L * _Msolar_to_kg_));
  _conversion_velocity_ = _conversion_dist_/_conversion_time_; // x [kpc/year] = x * _conversion_velocity_ [\hat{kpc}/[\hat{year}]] code units                                              
  //_Myear_ = 3.1556952e13L / tt;
  _Myear_ = 1.0e6 * _conversion_time_;  //1 mega year in code units
  _Gyear_ = 1.0e9 * _conversion_time_;  //1 Giga year in code units
  //_G_ = 1.0L;
  //_G_ = 6.67430e-11 / (_kpc_to_m_ * _kpc_to_m_ * _kpc_to_m_) * (_conversion_dist_*_conversion_dist_*_conversion_dist_)
  //     * _Msolar_to_kg_ * _year_to_s_ * _year_to_s_ * tt * tt;

  _G_ = 6.67430e-11 / (_kpc_to_m_ * _kpc_to_m_ * _kpc_to_m_) * (_conversion_dist_*_conversion_dist_*_conversion_dist_)
       * _Msolar_to_kg_ * _year_to_s_ * _year_to_s_ / (_conversion_time_ * _conversion_time_);      
  printf("_User_BoxSize_ = %f\n", (double)_User_BoxSize_);
  printf("G = %1.16e\n", (double)_G_);
  printf("kpc to code distance unit = %1.6e\n",(double)_conversion_dist_);
  printf("year to code time unit = %1.6e\n",(double)_conversion_time_); 
  printf("velocity kpc/year to code velocity unit = %1.6e\n",(double)_conversion_velocity_);
  printf("1 _Myear_ to code time units = %1.6e\n", (double)_Myear_);
  //printf("tt = %1.6e\n", (double)tt);
}

static void init_global_user_params(void)
{
  NUMBER_OF_THEADS = 1;
  BoxSize = 1.0;
  lmin = 5;                 // Coarset level of refinement
  lmax = lmin + 14;          // Finest level of refinement
  no_lmin_cell = 1 << lmin; // Number of cells in the lmin level of refinement
  no_lmin_cell_pow2 = no_lmin_cell * no_lmin_cell;
  no_lmin_cell_pow3 = no_lmin_cell * no_lmin_cell * no_lmin_cell;
  no_grid = no_lmin_cell + 1;
  //GL_no_ptcl_initial = 450000;
  GL_no_ptcl_initial = 100000;
  //GL_no_ptcl_initial = 10000;
  //GL_no_ptcl_initial = 1000;
  //GL_no_ptcl_initial = 2;
  GL_no_ptcl_final = GL_no_ptcl_initial;
  // GL_no_ptcl = 7550; // 2995865; // 299586; // 231299 // 298159
  //  GL_no_ptcl = 10000;
  //Maxdt = 100.0 *  _conversion_time_; // 1 year = 1.0e-6 Myear
  //Maxdt = 10.0 * _Gyear_; // 1 year = 1.0e-6 Myear
  Maxdt = 1.0 * _Gyear_; // 1 year = 1.0e-6 Myear
  printf("Code Max time to reach = %1.12e in code units of time\n",(double) Maxdt);
  //meanmass = 100; //Currently only used on input.c
  //  GL_total_mass_initial = GL_no_ptcl * meanmass;
  //  GL_total_mass_initial = 0;
  fr_output = 33; //114;
  MaxIterations = 100000000; // 100000000;
  no_grid_pow2 = no_grid * no_grid;
  no_grid_pow3 = no_grid * no_grid * no_grid; 
  bdry_cond_type = 1; // 0 = Periodic; 1 = Reflexive; 2 = Outflow   ##Note for Periodic conditions: Initial potential and acceleration are wrong, 
                                                                    // In fact, we shouldn need initial conditions because the periodicity 

  omp_set_num_threads(NUMBER_OF_THEADS);
}

static void init_global_ref_crit(void)
{
  ref_criterion_mass = 1.0e100; //1.0e100; // meanmass * 7;
  //ref_criterion_ptcl = 7;
  ref_criterion_ptcl = 4;
  //ref_criterion_ptcl = 1;
  n_exp = 1;   // n_exp = 0 is corrupted because particles can move between more than 1 level of refinement
  //_CFL_ = 0.9; // CFL criteria 0.5
  _CFL_ = 0.5; // CFL criteria 0.5
  _MAX_dt_ = 0.005 * Maxdt  ;//0.00333333333 * Maxdt ;// _Myear_ * 1;    Galaxy merger = 0.001 * Maxdt 
  //_MAX_dt_ = 1.0 * Maxdt  ;//0.00333333333 * Maxdt ;// _Myear_ * 1;    Galaxy merger = 0.001 * Maxdt 
  min_box_extra_size_per_side_user = 5;  // Minimum value in the extra cell per direcction (-x,+x,-y,+y,-z,+z) side of the boxes
  min_box_extra_size_per_side_real = min_box_extra_size_per_side_user > (n_exp - 1) ? 2 * min_box_extra_size_per_side_user : 2 * n_exp - 2;
}

static void init_global_poisson_params(void)
{
  //* >> Poisson parameters *//
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
  _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_ = 5000;
  // _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ = (1.0e-10);
  // _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_2 = (1.0e-10);


  _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ = (1.0e-10);
  _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_2 = (1.0e-10);
  check_poisson_error_method = 0; // Only used Gauss-Said or Jacobi in multigrid
  multigrid_cycle = 0;
  solverPreS = 0;
  solverPostS = 0;
  solverfinddphic = 0;
  _NiterPreS_ = 2;
  _NiterPostS_ = 2;
  _Niterfinddphic_ = 2;
  _Iter_branches_solver_ = 10; // 25;
  _w_SOR_ = 1.0;
  _w_SOR_HEAD_ = 1.0;

  iter_between_check_potential_solution = 5;

  branches_maximal_node_number_to_activate_conjugate_gradient = 0; //INT_MAX; // 513, 216 = node with minimum size of 1 (+1 n_exp) size, (1 + 2*n_exp)^3 * 8

  head_pot_method = 0;   // 0 = Multygrid, 1 = Conjugate gradient
  //branch_pot_method = 1; // 0 = SOR, 1 = Conjugate gradient
}

static void init_global_time_step_params(void)
{
  time_step_method = 1; // 0 = Only velocity, 1 = velocity + acceleration
}

static void init_global_force_params(void)
{
  force_stencil = 0; // 0 = 3-points, 1 = 5-points
}

static void init_global_energies_params(void)
{
  compute_observables_FLAG = true;
  potential_energy_type = 0; // 0 = Exact, 1 = approximation using potential grid
  // OBSERVABLES:
  //   Energies[0] = Kinetic
  //   Energies[1] = Potential
  //   Energies[2] = Total = K + P
  if (compute_observables_FLAG)
  {
    GL_energies = (vtype *)malloc(3 * sizeof(vtype));
    GL_momentum = (vtype *)malloc(3 * sizeof(vtype));
  }

  
}

static void init_global_ptcl(void)
{
  GL_ptcl_mass = (vtype *)calloc(GL_no_ptcl_initial, sizeof(vtype));
  GL_ptcl_x = (vtype *)calloc(GL_no_ptcl_initial, sizeof(vtype));
  GL_ptcl_y = (vtype *)calloc(GL_no_ptcl_initial, sizeof(vtype));
  GL_ptcl_z = (vtype *)calloc(GL_no_ptcl_initial, sizeof(vtype));
  GL_ptcl_vx = (vtype *)calloc(GL_no_ptcl_initial, sizeof(vtype));
  GL_ptcl_vy = (vtype *)calloc(GL_no_ptcl_initial, sizeof(vtype));
  GL_ptcl_vz = (vtype *)calloc(GL_no_ptcl_initial, sizeof(vtype));
  GL_ptcl_ax = (vtype *)calloc(GL_no_ptcl_initial, sizeof(vtype));
  GL_ptcl_ay = (vtype *)calloc(GL_no_ptcl_initial, sizeof(vtype));
  GL_ptcl_az = (vtype *)calloc(GL_no_ptcl_initial, sizeof(vtype));

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
  bder_os_sim = 1 > n_exp ? 1 : n_exp;            // Border outside the simulation box
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
  GL_times = (double *)calloc(100, sizeof(double));
}

static void init_global_memory(void)
{
  TOTAL_MEMORY_NODES = 0.0;
  TOTAL_MEMORY_CELDAS = 0.0;
  TOTAL_MEMORY_PARTICLES = 0.0;
  TOTAL_MEMORY_CELL_STRUCT = 0.0;
  TOTAL_MEMORY_CAJAS = 0.0;
  TOTAL_MEMORY_GRID_POINTS = 0.0;
  TOTAL_MEMORY_GRID_PROPERTIES = 0.0;
  TOTAL_MEMORY_AUX = 0.0;
  TOTAL_MEMORY_TENTACLES = 0.0;
  TOTAL_MEMORY_STACK = 0.0;
}

static void init_global_garbage_collector_parameters(void)
{
  Garbage_Collector_iter = 1000000; // Number of time-steps between each garbage collector
}

static void init_global_boundary_parameters(void)
{
  if (lmin < lmax && bdry_cond_type == 0)
  {
    GL_max_lv_ref_crosses_the_entire_simulation = lmin; // Maximum level of refinement such that it contains a zone of refinement which crosses the entire simulation
  }
}

static void init_multigrid2_parameters(void)
{
  int size, size_pow3;
  pp_phixx = (vtype **)malloc((lmin - 1) * sizeof(vtype *));
  pp_rhoxx = (vtype **)malloc((lmin - 1) * sizeof(vtype *));
  pp_restxx = (vtype **)malloc((lmin - 2) * sizeof(vtype *));
  pp_dphixx = (vtype **)malloc((lmin - 2) * sizeof(vtype *));
  zeros_xx = (vtype **)malloc((lmin - 1) * sizeof(vtype *));
  size = ((1 << (lmin)) + 1);
  size_pow3 = size * size * size;
  for (int i = 0; i < lmin - 2; i++)
  {
    pp_restxx[i] = (vtype *)calloc(size_pow3, sizeof(vtype));
    pp_dphixx[i] = (vtype *)calloc(size_pow3, sizeof(vtype));
    zeros_xx[i] = (vtype *)calloc(size_pow3, sizeof(vtype));
    size = ((1 << (lmin - i - 1)) + 1);
    size_pow3 = size * size * size;
    pp_phixx[i + 1] = (vtype *)malloc(size_pow3 * sizeof(vtype));
    pp_rhoxx[i + 1] = (vtype *)calloc(size_pow3, sizeof(vtype));
  }

  zeros_xx[lmin - 2] = (vtype *)calloc(size_pow3, sizeof(vtype));
}

void global_variables()
{
  //* >> Initializing constants *//
  init_global_constants();

  //* >> Initializing user parameters *//
  init_global_user_params();

  //* >> Initializing refinement criteria parameters *//
  init_global_ref_crit();

  //* >> Initializing Poisson parameters *//
  init_global_poisson_params();

  //* >> Initializing Time Step parameters *//
  init_global_time_step_params();

  //* >> Initializing Force parameters *//
  init_global_force_params();

  //* >> Initializing energy parameters *//
  init_global_energies_params();

  //* >> Initializing particles *//
  init_global_ptcl();

  //* >> Initializing tree head *//
  init_tree_head();

  //* >> Initializing border of simulation box parameters *//
  init_global_border_sim_box();

  //* >> Initializing output folder parameters *//
  init_global_folder_params();

  //* >> Initializing timer *//
  init_global_timer();

  //* >> Initializing memory parameters *//
  init_global_memory();

  //* >> Garbage Collector *//
  init_global_garbage_collector_parameters();

  //* >> Boundary parameters *//
  init_global_boundary_parameters();

  //* >> multigrid2 parameters *//
  init_multigrid2_parameters();
}
