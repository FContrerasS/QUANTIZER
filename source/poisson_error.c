/*
 * poisson_error.c
 *
 * Compute the error of the solution of the Poisson equation potential
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

#include "poisson_error.h"

//* >> Local Functions
static bool comparison_with_previous_solution(const struct node *ptr_node);
static bool poisson_error_mehod_0(const struct node *ptr_node);
static bool poisson_error_mehod_1(const struct node *ptr_node);
static bool poisson_error_mehod_2(const struct node *ptr_node);
static bool conj_grad_error(struct node *ptr_node, vtype *error, vtype error_total_pow2);

static bool comparison_with_previous_solution(const struct node *ptr_node)
{
  int box_grid_idx;
  vtype aux;

  vtype threshold_minus = 1 - _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_2;
  vtype threshold_plus = 1 + _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_2;

  bool check = true;

  for (int i = 0; i < ptr_node->grid_intr_size; i++)
  {
    box_grid_idx = ptr_node->ptr_intr_box_grid_idx[i];
    aux = ptr_node->ptr_pot[box_grid_idx] / ptr_node->ptr_pot_old[box_grid_idx];

    if ((threshold_minus > aux) || (threshold_plus < aux))
    {
      check = false;
      break;
    }
  }

  return check;
}

static bool poisson_error_mehod_0(const struct node *ptr_node)
{
  /* 	The condition used is option 1 (default one), and it
     uses the root mean square. In general, the error in the grid of the Poisson
     equation solution is computed as follows:

     error = Sqrt[ 1/Ng^3 * Sum[{( density_i - Laplacian(pot_i) )/rhomean }^2 ] ]

     The criterion is accepted when the error is less than the threshold. See the
     _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ parameter in
     Precision_Parameters.h. Note that error is an adimensional variable. */

  vtype error; // Error of the potential solution
  vtype diff;  // Diference between the Laplacian of the potential solution and the density

  int box_grid_idx;

  vtype H = 1.0L / (1 << ptr_node->lv);
  vtype one_over_H_pow_2 = 1.0L / (H * H);
  int size = ptr_node->cell_size; // Number of cells in the node

  vtype rhomean_times_4piG = 4 * _G_ * _PI_ * ptr_node->node_mass / (size * H * H * H);

  error = 0;

  int grid_box_real_dim_X = ptr_node->box_real_dim_x + 1;
  int grid_box_real_dim_X_times_Y = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

  // Computing the root mean square normalized to the mean density rhomean

  for (int i = 0; i < ptr_node->grid_intr_size; i++)
  {
    box_grid_idx = ptr_node->ptr_intr_box_grid_idx[i];
    diff = ptr_node->ptr_d[box_grid_idx] + one_over_H_pow_2 *
                                               (6.0 * ptr_node->ptr_pot[box_grid_idx] - ptr_node->ptr_pot[box_grid_idx + 1] - ptr_node->ptr_pot[box_grid_idx - 1] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X] - ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X_times_Y] - ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X_times_Y]);
    error += diff * diff;
  }

  error = error / ptr_node->grid_intr_size;
  error = sqrt(error);
  error = error / rhomean_times_4piG;
  //* >> If the precision condition satisfied *//
  if (error < _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_)
  {
    return true;
  }
  else
  {
    return false;
  }
}

static bool poisson_error_mehod_1(const struct node *ptr_node)
{

  /*    This module is responsible for checking the condition for the solution of the
     Poisson equation. The condition used is option 2, and it is the absolute
     error between the density and the laplacian for each grid point. In general,
     the error in each grid point of the Poisson equation solution is computed as
     follows:

     error_i =  |density_i - Laplacian(pot_i)| /rhomean

     The criterion is accepted only if all errors in the grid points are less than
     the threshold. See the _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ parameter in
     Precision_Parameters.h. Note that error is an adimensional variable. */

  vtype error; // Error of the potential solution

  int box_grid_idx;

  vtype H = 1.0L / (1 << ptr_node->lv);
  vtype one_over_H_pow_2 = 1.0L / (H * H);
  int size = ptr_node->cell_size; // Number of cells in the node

  vtype rhomean_times_4piG = 4 * _G_ * _PI_ * ptr_node->node_mass / (size * H * H * H);

  int grid_box_real_dim_X = ptr_node->box_real_dim_x + 1;
  int grid_box_real_dim_X_times_Y = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

  // Computing the root mean square normalized to the mean density rhomean
  for (int i = 0; i < ptr_node->grid_intr_size; i++)
  {
    box_grid_idx = ptr_node->ptr_intr_box_grid_idx[i];
    error = ptr_node->ptr_d[box_grid_idx] + one_over_H_pow_2 *
                                                (6.0 * ptr_node->ptr_pot[box_grid_idx] - ptr_node->ptr_pot[box_grid_idx + 1] - ptr_node->ptr_pot[box_grid_idx - 1] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X] - ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X_times_Y] - ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X_times_Y]);
    error = myabs(error) / rhomean_times_4piG;
    if (error >= _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_)
    {
      return false;
    }
  }

  return true;
}

static bool poisson_error_mehod_2(const struct node *ptr_node)
{
  /* 	   This module is responsible for checking the condition for the solution of the
   Poisson equation. The condition used is option 3, and it uses the root mean
   square. In general, the error in the grid of the Poisson equation solution
   is computed as follows:

   error =

   Sqrt[ 1/Ng^3 * Sum[{( density_i - Laplacian(pot_i) )/rhomean }^2 ] ] * 1/N^2

   The criterion is accepted when the error is less than the threshold. See the
   _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ parameter in
   Precision_Parameters.h. Note that error is an adimensional variable. */

  vtype error; // Error of the potential solution
  vtype diff;  // Diference between the Laplacian of the potential solution and the density

  int box_grid_idx;

  vtype H = 1.0L / (1 << ptr_node->lv);
  vtype one_over_H_pow_2 = 1.0L / (H * H);
  int size = ptr_node->cell_size; // Number of cells in the node

  vtype rhomean_times_4piG = 4 * _G_ * _PI_ * ptr_node->node_mass / (size * H * H * H);

  error = 0;

  int grid_box_real_dim_X = ptr_node->box_real_dim_x + 1;
  int grid_box_real_dim_X_times_Y = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

  // Computing the root mean square normalized to the mean density rhomean
  for (int i = 0; i < ptr_node->grid_intr_size; i++)
  {
    box_grid_idx = ptr_node->ptr_intr_box_grid_idx[i];
    diff = ptr_node->ptr_d[box_grid_idx] + one_over_H_pow_2 *
                                               (6.0 * ptr_node->ptr_pot[box_grid_idx] - ptr_node->ptr_pot[box_grid_idx + 1] - ptr_node->ptr_pot[box_grid_idx - 1] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X] - ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X_times_Y] - ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X_times_Y]);
    error += diff * diff;
  }

  error = error / ptr_node->grid_intr_size;
  error = sqrt(error);
  error = error / rhomean_times_4piG;
  // Following line is the only difference with poisson_error_mehod_0
  error = error / (cbrt(ptr_node->grid_intr_size) * cbrt(ptr_node->grid_intr_size));

  //* >> If the precision condition satisfied *//
  if (error < _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_)
  {
    return true;
  }
  else
  {
    return false;
  }
}

static bool conj_grad_error(struct node *ptr_node, vtype *error, vtype error_total_pow2)
{

  vtype H = 1.0L / (1 << ptr_node->lv);

  vtype rhomean_times_4piG = 4 * _G_ * _PI_ * ptr_node->node_mass / (ptr_node->cell_size * H * H * H);

  if (check_poisson_error_method == 0)
  {
    vtype aux_tol = ptr_node->grid_intr_size * (_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ * rhomean_times_4piG) * (_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ * rhomean_times_4piG); // sqrt(total_resid)/N/rhobar < error_tol
    // printf("lv = %d, ID = %d, cell_size = %d\n", ptr_node->lv, ptr_node->ID, ptr_node->cell_size);
    // printf("error_total_pow2 = %f, aux_tol = %f, grid_intr_size = %d, grid_bder_size = %d\n", (double)error_total_pow2, (double)aux_tol, ptr_node->grid_intr_size, ptr_node->grid_bder_size);

    if (error_total_pow2 < aux_tol)
    {
      return true;
    }
    return false;
  }
  else if (check_poisson_error_method == 1)
  {
    int box_grid_idx;
    vtype aux_error;
    // Computing the root mean square normalized to the mean density rhomean
    for (int i = 0; i < ptr_node->grid_intr_size; i++)
    {
      box_grid_idx = ptr_node->ptr_intr_box_grid_idx[i];
      aux_error = myabs(error[box_grid_idx]) / rhomean_times_4piG;
      if (aux_error >= _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_)
      {
        return false;
      }
    }
    return true;
  }
  else if (check_poisson_error_method == 2)
  {
    vtype aux_tol = ptr_node->grid_intr_size * (_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ * rhomean_times_4piG) * (_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_ * rhomean_times_4piG) * (cbrt(ptr_node->grid_intr_size) * cbrt(ptr_node->grid_intr_size));

    if (error_total_pow2 < aux_tol)
    {
      return true;
    }
    return false;
  }
  else
  {
    printf("Error: check_poisson_error_method method = %d is different to 0, 1 or 2\n", check_poisson_error_method);
    exit(EXIT_FAILURE);
  }

  return false;
}

bool poisson_error(struct node *ptr_node, vtype *error, vtype error_total_pow2, int type)
{

  bool check;

  check = comparison_with_previous_solution(ptr_node);

  // First checking with the previous solution
  if (check == true)
  {
    // Second checking with the error solution
    if (type == 0) // Gauss-Saidel and Jacobi
    {
      if (check_poisson_error_method == 0)
      {
        check = poisson_error_mehod_0(ptr_node);
      }
      else if (check_poisson_error_method == 1)
      {
        check = poisson_error_mehod_1(ptr_node);
      }
      else if (check_poisson_error_method == 2)
      {
        check = poisson_error_mehod_2(ptr_node);
      }
      else
      {
        printf("Error: check_poisson_error_method is different to 0, 1, 2 or 3\n");
        exit(EXIT_FAILURE);
      }
    }
    else if (type == 1) // Conjugate Gradient
    {
      check = conj_grad_error(ptr_node, error, error_total_pow2);
    }
    else
    {
      printf("Error: poisson_error type is different to 0, 1\n");
      exit(EXIT_FAILURE);
    }
  }

  if (check == false)
  {
    memcpy(ptr_node->ptr_pot_old, ptr_node->ptr_pot, ptr_node->grid_properties_cap * sizeof(vtype));
  }

  return check;
}