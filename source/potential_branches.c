/*
 * potential_branches.c
 *
 * Compute the potential in the branch nodes of the level of refinement
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

#include "potential_branches.h"

//* >> Local Functions
static int fill_red_and_black_branch_node(int **pptr_red_black, int *ptr_red_black_cap, int *ptr_red_black_size, const struct node *ptr_node);
static void computing_simulation_boundary_grid_point_potential(const struct node *ptr_node);
static void initial_potential_branch_node(const struct node *ptr_node);
static int compute_potential_branch_node(int **pptr_red_black, const int *ptr_red_black_size, struct node *ptr_node);

static int fill_red_and_black_branch_node(int **pptr_red_black, int *ptr_red_black_cap, int *ptr_red_black_size, const struct node *ptr_node)
{
  int red_size;     // Number of red points
  int black_size;   // Number of black points
  int red_or_black; // Desition if the point is red or black

  int box_grid_idx; // box grid index

  red_size = 0;
  black_size = 0;

  //* >> Filling Red and Black arrays with the correspondig grid points *//
  for (int i = 0; i < ptr_node->grid_intr_size; i++)
  {
    box_grid_idx = ptr_node->ptr_intr_box_grid_idx[i];

    red_or_black = ptr_node->ptr_intr_grid_cell_idx_x[i] +
                   ptr_node->ptr_intr_grid_cell_idx_y[i] +
                   ptr_node->ptr_intr_grid_cell_idx_z[i];

    //* >> Space checking of red and black pointers *//
    if (red_or_black % 2 == 0) // If it is red
    {
      red_size++;

      //* >> Space checking of the particle capacity of the red array *//
      if (space_check(&(ptr_red_black_cap[0]), red_size, 2.0f, "p1i1", &(pptr_red_black[0])) == _FAILURE_)
      {
        printf("Error, in space_check function\n");
        return _FAILURE_;
      }

      //* >> Adding element to the red array *//
      pptr_red_black[0][red_size - 1] = box_grid_idx;
    }
    else // If it is black
    {
      black_size++;

      //* >> Space checking of the particle capacity of the black array *//
      if (space_check(&(ptr_red_black_cap[1]), black_size, 2.0f, "p1i1", &(pptr_red_black[1])) == _FAILURE_)
      {
        printf("Error, in space_check function\n");
        return _FAILURE_;
      }

      //* >> Adding element to the black array *//
      pptr_red_black[1][black_size - 1] = box_grid_idx;
    }
  }

  //* >> Size of the red and black arrays *//
  ptr_red_black_size[0] = red_size;
  ptr_red_black_size[1] = black_size;

  return _SUCCESS_;
}

static void computing_simulation_boundary_grid_point_potential(const struct node *ptr_node)
{
  int lv = ptr_node->lv;

  int box_grid_idx; // Box grid index

  vtype aux_i;
  vtype aux_j;
  vtype aux_k;

  vtype H = 1.0 / (1 << lv);

  vtype dist;

  vtype aux_coeff = -_G_ * GL_total_mass_initial;

  //* >> Simulation Boundary grid points *//
  for (int i = 0; i < ptr_node->grid_sim_bdry_size; i++)
  {
    box_grid_idx = ptr_node->ptr_sim_bdry_box_grid_idx[i];

    if (ptr_node->ptr_sim_bdry_grid_cell_idx_x[i] == 0 ||
        ptr_node->ptr_sim_bdry_grid_cell_idx_x[i] == (1 << lv))
    {
      aux_i = 0.0;
    }
    else
    {
      aux_i = ptr_node->ptr_sim_bdry_grid_cell_idx_x[i] * H;
    }

    if (ptr_node->ptr_sim_bdry_grid_cell_idx_y[i] == 0 ||
        ptr_node->ptr_sim_bdry_grid_cell_idx_y[i] == (1 << lv))
    {
      aux_j = 0.0 ;
    }
    else
    {
      aux_j = ptr_node->ptr_sim_bdry_grid_cell_idx_y[i] * H;
    }

    if (ptr_node->ptr_sim_bdry_grid_cell_idx_z[i] == 0 ||
        ptr_node->ptr_sim_bdry_grid_cell_idx_z[i] == (1 << lv))
    {
      //aux_k = GL_cm[2] < 0.5 ? 0 : (1 << lv) * H;
      aux_k = 0.0;
    }
    else
    {
      aux_k = ptr_node->ptr_sim_bdry_grid_cell_idx_z[i] * H;
    }

    //dist = (x - GL_cm[0]) * (aux_i - GL_cm[0]) + (aux_j - GL_cm[1]) * (aux_j - GL_cm[1]) + (aux_k - GL_cm[2]) * (aux_k - GL_cm[2]);
    dist = (aux_i - 0.5) * (aux_i - 0.5) + (aux_j - 0.5) * (aux_j - 0.5) + (aux_k - 0.5) * (aux_k - 0.5);
    dist = mysqrt(dist);

    ptr_node->ptr_pot[box_grid_idx] = aux_coeff / dist;
  }
}

static void initial_potential_branch_node(const struct node *ptr_node)
{

  //printf("\n\npotential branch node \n\n");

  struct node *ptr_pt = ptr_node->ptr_pt;

  int box_grid_idx_node; // Node box grid index

  int box_idx_1_x_pt; // Parent box index at X direcction
  int box_idx_1_y_pt;
  int box_idx_1_z_pt;
  int box_idx_2_x_pt;
  int box_idx_2_y_pt;
  int box_idx_2_z_pt;

  int box_grid_idx_1_pt; // Parent box grid index 1
  int box_grid_idx_2_pt; // Parent box grid index 2
  int box_grid_idx_3_pt; // Parent box grid index 3
  int box_grid_idx_4_pt; // Parent box grid index 4
  int box_grid_idx_5_pt; // Parent box grid index 5
  int box_grid_idx_6_pt; // Parent box grid index 6
  int box_grid_idx_7_pt; // Parent box grid index 7
  int box_grid_idx_8_pt; // Parent box grid index 8

  int aux_idx_1_x;
  int aux_idx_1_y;
  int aux_idx_1_z;
  int aux_idx_2_x;
  int aux_idx_2_y;
  int aux_idx_2_z;

  // vtype aux_pot; // Auxiliar Potential in the grid point

  int lv = ptr_node->lv;

  // vtype H = 1.0 / (1 << lv);
  // vtype aux_x, aux_y, aux_z, dist;

  int grid_box_real_dim_X_pt = (ptr_pt->box_real_dim_x + 1);
  int grid_box_real_dim_X_times_Y_pt = (ptr_pt->box_real_dim_x + 1) * (ptr_pt->box_real_dim_y + 1);

  // printf("lv = %d, child grid dim x, y , z = %d, %d, %d\n ",ptr_node->lv,(ptr_node->box_real_dim_x + 1),(ptr_node->box_real_dim_y + 1),(ptr_node->box_real_dim_z + 1));
  // printf("number of cells = %d\n",ptr_node->cell_size);

  //* >> Simulation Boundary grid points *//
  computing_simulation_boundary_grid_point_potential(ptr_node);

  //* >> Passing the potential from coarse parent to fine child *//
  // Border grid points
  //printf("booundary points \n");
  for (int i = 0; i < ptr_node->grid_bdry_size; i++)
  {
    //* >> Child box indices *//
    box_grid_idx_node = ptr_node->ptr_bdry_box_grid_idx[i];

    //* >> Parent box indices *//
    // box_idx_x0_node = (ptr_ch->ptr_bdry_grid_cell_idx_x[i] >> 1) - ptr_node->box_ts_x;
    // box_idx_x1_node = ((ptr_ch->ptr_bdry_grid_cell_idx_x[i] + 1) >> 1) - ptr_node->box_ts_x;
    // box_idx_y0_node = (ptr_ch->ptr_bdry_grid_cell_idx_y[i] >> 1) - ptr_node->box_ts_y;
    // box_idx_y1_node = ((ptr_ch->ptr_bdry_grid_cell_idx_y[i] + 1) >> 1) - ptr_node->box_ts_y;
    // box_idx_z0_node = (ptr_ch->ptr_bdry_grid_cell_idx_z[i] >> 1) - ptr_node->box_ts_z;
    // box_idx_z1_node = ((ptr_ch->ptr_bdry_grid_cell_idx_z[i] + 1) >> 1) - ptr_node->box_ts_z;

    //* >> Parent box indices *//
    aux_idx_1_x = ptr_node->ptr_bdry_grid_cell_idx_x[i] >> 1;
    aux_idx_1_y = ptr_node->ptr_bdry_grid_cell_idx_y[i] >> 1;
    aux_idx_1_z = ptr_node->ptr_bdry_grid_cell_idx_z[i] >> 1;
    aux_idx_2_x = (ptr_node->ptr_bdry_grid_cell_idx_x[i] + 1) >> 1;
    aux_idx_2_y = (ptr_node->ptr_bdry_grid_cell_idx_y[i] + 1) >> 1;
    aux_idx_2_z = (ptr_node->ptr_bdry_grid_cell_idx_z[i] + 1) >> 1;

    box_idx_1_x_pt = aux_idx_1_x - ptr_pt->box_ts_x;
    box_idx_1_y_pt = aux_idx_1_y - ptr_pt->box_ts_y;
    box_idx_1_z_pt = aux_idx_1_z - ptr_pt->box_ts_z;
    box_idx_2_x_pt = aux_idx_2_x - ptr_pt->box_ts_x;
    box_idx_2_y_pt = aux_idx_2_y - ptr_pt->box_ts_y;
    box_idx_2_z_pt = aux_idx_2_z - ptr_pt->box_ts_z;
    if (ptr_pt->pbc_crosses_sim_box_bdry == true)
    {
      
      if (aux_idx_1_x > ptr_pt->box_max_x)
      {
        box_idx_1_x_pt -= (1 << (lv - 1));
        box_idx_2_x_pt -= (1 << (lv - 1));
      }

      if (aux_idx_1_y > ptr_pt->box_max_y)
      {
        box_idx_1_y_pt -= (1 << (lv - 1));
        box_idx_2_y_pt -= (1 << (lv - 1));
      }

      if (aux_idx_1_z > ptr_pt->box_max_z)
      {
        box_idx_1_z_pt -= (1 << (lv - 1));
        box_idx_2_z_pt -= (1 << (lv - 1));
      }
    }

    //* >> Parent grid indices *//

    box_grid_idx_1_pt = box_idx_1_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
    box_grid_idx_2_pt = box_idx_2_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
    box_grid_idx_3_pt = box_idx_1_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
    box_grid_idx_4_pt = box_idx_2_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
    box_grid_idx_5_pt = box_idx_1_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;
    box_grid_idx_6_pt = box_idx_2_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;
    box_grid_idx_7_pt = box_idx_1_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;
    box_grid_idx_8_pt = box_idx_2_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;

    ptr_node->ptr_pot[box_grid_idx_node] = 0.125 * (ptr_pt->ptr_pot[box_grid_idx_1_pt] + ptr_pt->ptr_pot[box_grid_idx_2_pt] +
                                                    ptr_pt->ptr_pot[box_grid_idx_3_pt] + ptr_pt->ptr_pot[box_grid_idx_4_pt] +
                                                    ptr_pt->ptr_pot[box_grid_idx_5_pt] + ptr_pt->ptr_pot[box_grid_idx_6_pt] +
                                                    ptr_pt->ptr_pot[box_grid_idx_7_pt] + ptr_pt->ptr_pot[box_grid_idx_8_pt]);
    //printf("box grid idx = %d, pot = %.12f\n",box_grid_idx_node,(double)ptr_node->ptr_pot[box_grid_idx_node] );
  }

  // Interior grid points
  //printf("interior points\n");
  for (int i = 0; i < ptr_node->grid_intr_size; i++)
  {
    //* >> Child box indices *//
    box_grid_idx_node = ptr_node->ptr_intr_box_grid_idx[i];

    //* >> Parent box indices *//
    aux_idx_1_x = ptr_node->ptr_intr_grid_cell_idx_x[i] >> 1;
    aux_idx_1_y = ptr_node->ptr_intr_grid_cell_idx_y[i] >> 1;
    aux_idx_1_z = ptr_node->ptr_intr_grid_cell_idx_z[i] >> 1;
    aux_idx_2_x = (ptr_node->ptr_intr_grid_cell_idx_x[i] + 1) >> 1;
    aux_idx_2_y = (ptr_node->ptr_intr_grid_cell_idx_y[i] + 1) >> 1;
    aux_idx_2_z = (ptr_node->ptr_intr_grid_cell_idx_z[i] + 1) >> 1;

    box_idx_1_x_pt = aux_idx_1_x - ptr_pt->box_ts_x;
    box_idx_1_y_pt = aux_idx_1_y - ptr_pt->box_ts_y;
    box_idx_1_z_pt = aux_idx_1_z - ptr_pt->box_ts_z;
    box_idx_2_x_pt = aux_idx_2_x - ptr_pt->box_ts_x;
    box_idx_2_y_pt = aux_idx_2_y - ptr_pt->box_ts_y;
    box_idx_2_z_pt = aux_idx_2_z - ptr_pt->box_ts_z;

    if (ptr_pt->pbc_crosses_sim_box_bdry == true)
    {
      if (aux_idx_1_x > ptr_pt->box_max_x)
      {
        box_idx_1_x_pt -= (1 << (lv - 1));
        box_idx_2_x_pt -= (1 << (lv - 1));
      }

      if (aux_idx_1_y > ptr_pt->box_max_y)
      {
        box_idx_1_y_pt -= (1 << (lv - 1));
        box_idx_2_y_pt -= (1 << (lv - 1));
      }

      if (aux_idx_1_z > ptr_pt->box_max_z)
      {
        box_idx_1_z_pt -= (1 << (lv - 1));
        box_idx_2_z_pt -= (1 << (lv - 1));
      }
    }

    //* >> Parent grid indices *//

    box_grid_idx_1_pt = box_idx_1_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
    box_grid_idx_2_pt = box_idx_2_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
    box_grid_idx_3_pt = box_idx_1_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
    box_grid_idx_4_pt = box_idx_2_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_1_z_pt * grid_box_real_dim_X_times_Y_pt;
    box_grid_idx_5_pt = box_idx_1_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;
    box_grid_idx_6_pt = box_idx_2_x_pt + box_idx_1_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;
    box_grid_idx_7_pt = box_idx_1_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;
    box_grid_idx_8_pt = box_idx_2_x_pt + box_idx_2_y_pt * grid_box_real_dim_X_pt + box_idx_2_z_pt * grid_box_real_dim_X_times_Y_pt;

    ptr_node->ptr_pot[box_grid_idx_node] = 0.125 * (ptr_pt->ptr_pot[box_grid_idx_1_pt] + ptr_pt->ptr_pot[box_grid_idx_2_pt] +
                                                    ptr_pt->ptr_pot[box_grid_idx_3_pt] + ptr_pt->ptr_pot[box_grid_idx_4_pt] +
                                                    ptr_pt->ptr_pot[box_grid_idx_5_pt] + ptr_pt->ptr_pot[box_grid_idx_6_pt] +
                                                    ptr_pt->ptr_pot[box_grid_idx_7_pt] + ptr_pt->ptr_pot[box_grid_idx_8_pt]);
    
    //printf("box grid idx = %d, pot = %.12f\n",box_grid_idx_node,(double)ptr_node->ptr_pot[box_grid_idx_node] );
  }
}

static int compute_potential_branch_node(int **pptr_red_black, const int *ptr_red_black_size, struct node *ptr_node)
{

  //printf("calling\n");

  int red_size;
  int black_size;
  vtype aux_pot;

  red_size = ptr_red_black_size[0];
  black_size = ptr_red_black_size[1];

  int grid_box_real_dim_X_node = (ptr_node->box_real_dim_x + 1);
  int grid_box_real_dim_X_times_Y_node = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

  vtype H = 1.0 / (1 << ptr_node->lv);
  vtype H_pow2 = H * H;

  //* >> Cycle over the Successive over-relaxation *//
  for (int iter = 0; iter < _Iter_branches_solver_; iter++)
  {
    //* >> Cycle over red points *//
    int box_grid_idx;
    int box_grid_idxNbr_x_plus;  // Box grid index in the neigborhood on the right
    int box_grid_idxNbr_x_minus; // Box grid index in the neigborhood on the left
    int box_grid_idxNbr_y_plus;  // Box grid index in the neigborhood behind
    int box_grid_idxNbr_y_minus; // Box grid index in the neigborhood in front
    int box_grid_idxNbr_z_plus;  // Box grid index in the neigborhood up
    int box_grid_idxNbr_z_minus; // Box grid index in the neigborhood down

    #pragma omp parallel for private(aux_pot,box_grid_idx,box_grid_idxNbr_x_plus,box_grid_idxNbr_x_minus,box_grid_idxNbr_y_plus,box_grid_idxNbr_y_minus,box_grid_idxNbr_z_plus,box_grid_idxNbr_z_minus)
    for (int i = 0; i < red_size; i++)
    {
      box_grid_idx = pptr_red_black[0][i];
      box_grid_idxNbr_x_plus = box_grid_idx + 1;
      box_grid_idxNbr_x_minus = box_grid_idx - 1;
      box_grid_idxNbr_y_plus = box_grid_idx + grid_box_real_dim_X_node;
      box_grid_idxNbr_y_minus = box_grid_idx - grid_box_real_dim_X_node;
      box_grid_idxNbr_z_plus = box_grid_idx + grid_box_real_dim_X_times_Y_node;
      box_grid_idxNbr_z_minus = box_grid_idx - grid_box_real_dim_X_times_Y_node;
      
      aux_pot = 1.0L / 6.0L * (-H_pow2 * ptr_node->ptr_d[box_grid_idx] + ptr_node->ptr_pot[box_grid_idxNbr_x_plus] + ptr_node->ptr_pot[box_grid_idxNbr_x_minus] + ptr_node->ptr_pot[box_grid_idxNbr_y_plus] + ptr_node->ptr_pot[box_grid_idxNbr_y_minus] + ptr_node->ptr_pot[box_grid_idxNbr_z_plus] + ptr_node->ptr_pot[box_grid_idxNbr_z_minus]);
      aux_pot = _w_SOR_ * aux_pot + (1.0 - _w_SOR_) * ptr_node->ptr_pot[box_grid_idx];
      ptr_node->ptr_pot[box_grid_idx] = aux_pot;
    }

    //* >> Cycle over black points *//
    #pragma omp parallel for private(aux_pot,box_grid_idx,box_grid_idxNbr_x_plus,box_grid_idxNbr_x_minus,box_grid_idxNbr_y_plus,box_grid_idxNbr_y_minus,box_grid_idxNbr_z_plus,box_grid_idxNbr_z_minus)
    for (int i = 0; i < black_size; i++)
    {
      box_grid_idx = pptr_red_black[1][i];
      box_grid_idxNbr_x_plus = box_grid_idx + 1;
      box_grid_idxNbr_x_minus = box_grid_idx - 1;
      box_grid_idxNbr_y_plus = box_grid_idx + grid_box_real_dim_X_node;
      box_grid_idxNbr_y_minus = box_grid_idx - grid_box_real_dim_X_node;
      box_grid_idxNbr_z_plus = box_grid_idx + grid_box_real_dim_X_times_Y_node;
      box_grid_idxNbr_z_minus = box_grid_idx - grid_box_real_dim_X_times_Y_node;
      
      aux_pot = 1.0L / 6.0L * (-H_pow2 * ptr_node->ptr_d[box_grid_idx] + ptr_node->ptr_pot[box_grid_idxNbr_x_plus] + ptr_node->ptr_pot[box_grid_idxNbr_x_minus] + ptr_node->ptr_pot[box_grid_idxNbr_y_plus] + ptr_node->ptr_pot[box_grid_idxNbr_y_minus] + ptr_node->ptr_pot[box_grid_idxNbr_z_plus] + ptr_node->ptr_pot[box_grid_idxNbr_z_minus]);
      aux_pot = _w_SOR_ * aux_pot + (1.0 - _w_SOR_) * ptr_node->ptr_pot[box_grid_idx];
      ptr_node->ptr_pot[box_grid_idx] = aux_pot;
    }
  }

  return _SUCCESS_;
}

int potential_branches(void)
{

  //int number_of_threads = 8;

  //omp_set_num_threads(number_of_threads);

  //struct node *ptr_node;
  
  struct node *ptr_ch;

  //omp_set_num_threads(16);

  int no_pts; // Number of parents in the cycle

  int **pptr_red_black;
  int *ptr_red_black_cap;
  int *ptr_red_black_size;

  int iter;
  bool check;

  //clock_t aux_clock;
  //double aux_clock_omp;
  struct timespec start, finish;

  vtype *dummy_pvtype = NULL;
  vtype dummy_vtype = 0.0;

  // // Red and Black arrays contain the index element of the box grid of the potential
  pptr_red_black = (int **)malloc(2 * sizeof(int *));
  pptr_red_black[0] = NULL;
  pptr_red_black[1] = NULL;
  ptr_red_black_cap = (int *)malloc(2 * sizeof(int));
  ptr_red_black_cap[0] = 0;
  ptr_red_black_cap[1] = 0;
  ptr_red_black_size = (int *)malloc(2 * sizeof(int));
  ptr_red_black_size[0] = 0;
  ptr_red_black_size[1] = 0;

  for (int lv = 1; lv <= GL_tentacles_level_max; lv++)
  {
    no_pts = GL_tentacles_size[lv];
    
    for (int j = 0; j < no_pts; j++)
    {
      ptr_ch = GL_tentacles[lv][j];
      clock_gettime( CLOCK_REALTIME, &start);
      initial_potential_branch_node(ptr_ch);
      clock_gettime( CLOCK_REALTIME, &finish);
      GL_times[22] += ( finish.tv_sec - start.tv_sec ) + ( finish.tv_nsec - start.tv_nsec )/ 1000000000.;

      //memcpy(ptr_ch->ptr_pot_old, ptr_ch->ptr_pot, ptr_ch->grid_properties_cap * sizeof(vtype));
      #pragma omp parallel for
      for(int hh = 0; hh<ptr_ch->grid_properties_cap; hh++ )
      {
        ptr_ch->ptr_pot_old[hh] = ptr_ch->ptr_pot[hh];
      }

      if (ptr_ch->cell_size > branches_maximal_node_number_to_activate_conjugate_gradient) // 513, 216 = node with minimum size of 1 (+1 n_exp) size, (1 + 2*n_exp)^3 * 8
      {
        clock_gettime( CLOCK_REALTIME, &start);
        if (conjugate_gradient(ptr_ch) == _FAILURE_)
        {
          printf("Error at function conjugate_gradient()\n");
          return _FAILURE_;
        }
        clock_gettime( CLOCK_REALTIME, &finish);
        GL_times[24] += ( finish.tv_sec - start.tv_sec ) + ( finish.tv_nsec - start.tv_nsec )/ 1000000000.;
      }
      else
      {
        //* >> Computing the potential in the child node *//
        clock_gettime( CLOCK_REALTIME, &start);
        if (fill_red_and_black_branch_node(pptr_red_black, ptr_red_black_cap, ptr_red_black_size, ptr_ch) == _FAILURE_)
        {
          printf("Error at function fill_red_and_black()\n");
          return _FAILURE_;
        }
        clock_gettime( CLOCK_REALTIME, &finish);
        GL_times[23] += ( finish.tv_sec - start.tv_sec ) + ( finish.tv_nsec - start.tv_nsec )/ 1000000000.;

        check = false;

        iter = 0;
        while (iter < _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_ && check == false)
        {
          iter = iter + 1;
          //* >> Computing the potential in the child node *//
          clock_gettime( CLOCK_REALTIME, &start);
          if (compute_potential_branch_node(pptr_red_black, ptr_red_black_size, ptr_ch) == _FAILURE_)
          {
            printf("Error at function potential_branch_node()\n");
            return _FAILURE_;
          }
          clock_gettime( CLOCK_REALTIME, &finish);
          GL_times[24] += ( finish.tv_sec - start.tv_sec ) + ( finish.tv_nsec - start.tv_nsec )/ 1000000000.;
          if (iter % iter_between_check_potential_solution == 0)
          {
            //* >> CHEKING ERROR SOLUTION CONDITION *//
            clock_gettime( CLOCK_REALTIME, &start);
            check = poisson_error(ptr_ch, dummy_pvtype, dummy_vtype, 0);
            clock_gettime( CLOCK_REALTIME, &finish);
            GL_times[25] += ( finish.tv_sec - start.tv_sec ) + ( finish.tv_nsec - start.tv_nsec )/ 1000000000.;
          }
        }
        if (iter == _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_)
        {
          printf("\nERROR: The precision was not reached in the branch node. Too many SOR Iterations, plz choose a lower precision than %1.3e\n", (double)_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_);
          return _FAILURE_;
        }
      }

    } // End cycle over number of children
      //} //  End cycle over number of parents
  }

  // //* >> Free pointers *//
  if (pptr_red_black[0] != NULL)
  {
    free(pptr_red_black[0]); // Free red
  }
  if (pptr_red_black[1] != NULL)
  {
    free(pptr_red_black[1]); // Free black
  }
  free(pptr_red_black);     // Free red and black
  free(ptr_red_black_cap);  // Free capacity of red and black
  free(ptr_red_black_size); // Free size of red and black

  return _SUCCESS_;
}
