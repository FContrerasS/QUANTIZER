/*
 * potential_head_node.c
 *
 * Compute the potential in the head node or coarsest level of refinement
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

#include "potential_head_node.h"

//* >> Local Functions
static int compute_potential_head_node(struct node *ptr_node_pt);

static int compute_potential_head_node(struct node *ptr_node_pt)
{

  int grid_side; // Number of grid points per side
  int grid_side_pow2;
  int grid_limit; // Superior border in the internal grid box

  vtype *aux_pot;
  vtype *aux_d;
  int size;
  int size_pow2;
  int size_pow3;

  int aux_idx;
  int box_grid_idx;

  size = no_lmin_cell + 1;
  size_pow2 = size * size;
  size_pow3 = size_pow2 * size;
  aux_pot = (vtype *)malloc(size_pow3 * sizeof(vtype));
  aux_d = (vtype *)malloc(size_pow3 * sizeof(vtype));

  grid_side = box_side_lmin + 1;
  grid_side_pow2 = grid_side * grid_side;
  grid_limit = grid_side - bder_os_sim - 1;

  //* >> Pasing information from real parameters to auxiliary parameters *//
  for (int k = bder_os_sim; k < grid_limit + 1; k++)
  {
    for (int j = bder_os_sim; j < grid_limit + 1; j++)
    {
      for (int i = bder_os_sim; i < grid_limit + 1; i++)
      {
        box_grid_idx = i + j * grid_side + k * grid_side_pow2;
        aux_idx = (i - bder_os_sim) + (j - bder_os_sim) * size + (k - bder_os_sim) * size_pow2;
        aux_pot[aux_idx] = ptr_node_pt->ptr_pot[box_grid_idx];
        aux_d[aux_idx] = ptr_node_pt->ptr_d[box_grid_idx];
      }
    }
  }

  int type_multigrid = multigrid_cycle;
  //* >> Poisson Solver *//
  // Potential-Density-gridsize-type of multigrid
  multigrid(aux_pot, aux_d, size, type_multigrid);
  // multigrid2(aux_pot, aux_d, type_multigrid);

  //* >> Pasing information from auxiliary parameters to real parameters *//
  for (int k = bder_os_sim; k < grid_limit + 1; k++)
  {
    for (int j = bder_os_sim; j < grid_limit + 1; j++)
    {
      for (int i = bder_os_sim; i < grid_limit + 1; i++)
      {
        box_grid_idx = i + j * grid_side + k * grid_side_pow2;
        aux_idx = (i - bder_os_sim) + (j - bder_os_sim) * size + (k - bder_os_sim) * size_pow2;
        ptr_node_pt->ptr_pot[box_grid_idx] = aux_pot[aux_idx];
        ptr_node_pt->ptr_d[box_grid_idx] = aux_d[aux_idx];
      }
    }
  }
  free(aux_pot);
  free(aux_d);

  return _SUCCESS_;
}

int potential_head_node(void)
{

  //clock_t aux_clock;
  struct timespec start, finish;

  struct node *ptr_node = GL_ptr_tree;

  vtype *dummy_pvtype = NULL;
  vtype dummy_vtype = 0.0;

  memcpy(ptr_node->ptr_pot_old, ptr_node->ptr_pot, ptr_node->grid_properties_cap * sizeof(vtype));

  if (head_pot_method == 0)
  {
    int iter = 0;
    bool check;

    // //* >> CHEKING ERROR SOLUTION CONDITION *//
    // aux_clock = clock();
    // check = poisson_error(ptr_node, dummy_pvtype, dummy_vtype,0);
    // GL_times[21] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

    check = false;

    //* >> SOLVING POISSON EQUATION *//
    while (iter < _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_ && check == false)
    {
      iter = iter + 1;
      //* >> INITIAL POTENTIAL COMPUTATION *//
      //aux_clock = clock();
      clock_gettime( CLOCK_REALTIME, &start);
      if (compute_potential_head_node(ptr_node) == _FAILURE_)
      {
        printf("\n\n Error running potential_head_node() function \n\n ");
        return _FAILURE_;
      }
      //GL_times[20] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
      clock_gettime( CLOCK_REALTIME, &finish);
      GL_times[20] += ( finish.tv_sec - start.tv_sec ) + ( finish.tv_nsec - start.tv_nsec )/ 1000000000.;
      if (iter % iter_between_check_potential_solution == 0)
      {
        //* >> CHEKING ERROR SOLUTION CONDITION *//
        //aux_clock = clock();
        clock_gettime( CLOCK_REALTIME, &start);
        check = poisson_error(ptr_node, dummy_pvtype, dummy_vtype, 0);
        //GL_times[21] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
        clock_gettime( CLOCK_REALTIME, &finish);
        GL_times[21] += ( finish.tv_sec - start.tv_sec ) + ( finish.tv_nsec - start.tv_nsec )/ 1000000000.;
      }
    }

    if (iter == _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_)
    {
      printf("\nERROR: The precision was not reached in the parent node. Too many Multigrid Iterations, plz choose a lower precision than %1.3e\n", (double)_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_);
      return _FAILURE_;
    }
  }
  else if (head_pot_method == 1)
  {
    //printf("hola\n");
    //aux_clock = clock();
    clock_gettime( CLOCK_REALTIME, &start);
    //* >> SOLVING POISSON EQUATION *//
    if (conjugate_gradient(ptr_node) == _FAILURE_)
    {
      printf("\n\n Error running conjugate_gradient() function \n\n ");
      return _FAILURE_;
    }
    //GL_times[20] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
    clock_gettime( CLOCK_REALTIME, &finish);
    GL_times[20] += ( finish.tv_sec - start.tv_sec ) + ( finish.tv_nsec - start.tv_nsec )/ 1000000000.;
  }
  else
  {
    printf("Error, the Head potential method is equal to %d\n", head_pot_method);
    return _FAILURE_;
  }

  return _SUCCESS_;
}