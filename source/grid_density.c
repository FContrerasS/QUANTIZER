/*
 * grid_density.c
 *
 * Compute the density of the grid in each node
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

#include "grid_density.h"

//* >> Local Functions
static void computing_grid_density_head(struct node *ptr_head);
static void computing_grid_density_branch(struct node *ptr_node);

static void computing_grid_density_head(struct node *ptr_head)
{
  int box_grid_idx_x; // grid_idx component in the grid box
  int box_grid_idx_y;
  int box_grid_idx_z;
  int box_grid_idx;    // Grid box grid_idx
  int box_grid_idxNbr; // Box index in the neigborhood

  vtype pos_x; // Particle position in the grid level
  vtype pos_y;
  vtype pos_z;

  int pos_x_floor; // floor of the particle position in the grid level
  int pos_y_floor;
  int pos_z_floor;

  vtype w_x_1; // Weight component of the CIC method. It corresponds to the distance between the particle and the grid point
  vtype w_y_1;
  vtype w_z_1;
  vtype w_x_2;
  vtype w_y_2;
  vtype w_z_2;

  vtype w[8]; // Weight of the CIC method

  int lv = ptr_head->lv;                                // Level of refinement
  vtype H = 1.0L / (1 << lv);                           // Size of the grid side
  vtype poisson_coeff = 4.0 * _G_ * _PI_ / (H * H * H); // Poisson coefficient

  int grid_box_real_dim_X = (ptr_head->box_real_dim_x + 1);
  int grid_box_real_dim_X_times_Y = (ptr_head->box_real_dim_x + 1) * (ptr_head->box_real_dim_y + 1);

  for (int i = 0; i < GL_no_ptcl_final; i++)
  {
    //* >> Position of the particles in the grid level *//
    pos_x = GL_ptcl_x[i] * (1 << lv);
    pos_y = GL_ptcl_y[i] * (1 << lv);
    pos_z = GL_ptcl_z[i] * (1 << lv);

    // //* >> floor of the particles positions in the grid level *//
    // pos_x_floor = myfloor(pos_x);
    // pos_y_floor = myfloor(pos_y);
    // pos_z_floor = myfloor(pos_z);
    pos_x_floor = (int)pos_x;
    pos_y_floor = (int)pos_y;
    pos_z_floor = (int)pos_z;

    //* >> Computing the weights of the nearest grid points of the particle *//
    // Each for cyle yields 2 options: X or 1-X, where X =  pos_x - pos_x_floor

    w_x_1 = pos_x - pos_x_floor;
    w_y_1 = pos_y - pos_y_floor;
    w_z_1 = pos_z - pos_z_floor;
    w_x_2 = 1 - w_x_1;
    w_y_2 = 1 - w_y_1;
    w_z_2 = 1 - w_z_1;
    w[0] = w_x_2 * w_y_2 * w_z_2;
    w[1] = w_x_1 * w_y_2 * w_z_2;
    w[2] = w_x_2 * w_y_1 * w_z_2;
    w[3] = w_x_1 * w_y_1 * w_z_2;
    w[4] = w_x_2 * w_y_2 * w_z_1;
    w[5] = w_x_1 * w_y_2 * w_z_1;
    w[6] = w_x_2 * w_y_1 * w_z_1;
    w[7] = w_x_1 * w_y_1 * w_z_1;

    box_grid_idx_x = pos_x_floor - ptr_head->box_ts_x;
    box_grid_idx_y = pos_y_floor - ptr_head->box_ts_y;
    box_grid_idx_z = pos_z_floor - ptr_head->box_ts_z;
    box_grid_idx = box_grid_idx_x + box_grid_idx_y * grid_box_real_dim_X + box_grid_idx_z * grid_box_real_dim_X_times_Y;
    //* >> Particle density contributes to 8 enclosure grid points *//
    for (int kk = 0; kk < 2; kk++)
    {
      for (int jj = 0; jj < 2; jj++)
      {
        for (int ii = 0; ii < 2; ii++)
        {
          box_grid_idxNbr = box_grid_idx + ii + jj * grid_box_real_dim_X + kk * grid_box_real_dim_X_times_Y;
          ptr_head->ptr_d[box_grid_idxNbr] += poisson_coeff * GL_ptcl_mass[i] * w[ii + 2 * jj + 4 * kk];
        }
      }
    }
  }
}

static void computing_grid_density_branch(struct node *ptr_node)
{
  int box_grid_idx_x; // grid_idx component in the grid box
  int box_grid_idx_y;
  int box_grid_idx_z;
  int box_grid_idx; // Grid box grid_idx

  int box_idx;  // Box index of the node cell
  int ptcl_idx; // Particle grid_idx in the node

  vtype pos_x; // Particle position in the grid level
  vtype pos_y;
  vtype pos_z;

  int pos_x_floor; // floor of the particle position in the grid level
  int pos_y_floor;
  int pos_z_floor;

  vtype pos_x_rel;
  vtype pos_y_rel;
  vtype pos_z_rel;

  vtype w_x_1; // Weight component of the CIC method. It corresponds to the distance between the particle and the grid point
  vtype w_y_1;
  vtype w_z_1;
  vtype w_x_2;
  vtype w_y_2;
  vtype w_z_2;

  vtype w[8]; // Weight of the CIC method

  int lv = ptr_node->lv;      // Level of refinement
  vtype H = 1.0L / (1 << lv); // Size of the grid side
  vtype one_over_H = (1 << lv);
  // vtype poisson_coeff = 4 * _G_ * _PI_ / (H * H * H); // Poisson coefficient
  vtype poisson_coeff_norm = 4.0 * _G_ * _PI_ / (H * H * H * H * H * H); // Poisson coefficient normalized

  int grid_box_real_dim_X = (ptr_node->box_real_dim_x + 1);
  int grid_box_real_dim_X_times_Y = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

  int box_grid_idxNbr[8];
  int aux_int;

  int counter_ptcl = 0;
  int total_ptcl = ptr_node->local_no_ptcl_full_node;
  int cell_ptcl;
  int cell_idx = -1;

  while (counter_ptcl < total_ptcl)
  {
    cell_idx++;
    box_idx = ptr_node->ptr_box_idx[cell_idx];
    cell_ptcl = ptr_node->ptr_cell_struct[box_idx].ptcl_size;

    if (cell_ptcl > 0)
    {
      ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[0];

      //* >> Position of the particles in the grid level *//
      pos_x = GL_ptcl_x[ptcl_idx] * one_over_H;
      pos_y = GL_ptcl_y[ptcl_idx] * one_over_H;
      pos_z = GL_ptcl_z[ptcl_idx] * one_over_H;

      // //* >> floor of the particles positions in the grid level *//
      // pos_x_floor = myfloor(pos_x);
      // pos_y_floor = myfloor(pos_y);
      // pos_z_floor = myfloor(pos_z);
      pos_x_floor = (int)pos_x;
      pos_y_floor = (int)pos_y;
      pos_z_floor = (int)pos_z;

      box_grid_idx_x = pos_x_floor - ptr_node->box_ts_x;
      box_grid_idx_y = pos_y_floor - ptr_node->box_ts_y;
      box_grid_idx_z = pos_z_floor - ptr_node->box_ts_z;

      if (ptr_node->pbc_crosses_the_boundary_simulation_box == true)
      {
        if (pos_x_floor > ptr_node->box_max_x)
        {
          box_grid_idx_x -= (1 << lv);
        }

        if (pos_y_floor > ptr_node->box_max_y)
        {
          box_grid_idx_y -= (1 << lv);
        }

        if (pos_z_floor > ptr_node->box_max_z)
        {
          box_grid_idx_z -= (1 << lv);
        }
      }

      box_grid_idx = box_grid_idx_x + box_grid_idx_y * grid_box_real_dim_X + box_grid_idx_z * grid_box_real_dim_X_times_Y;

      pos_x_rel = H * pos_x_floor;
      pos_y_rel = H * pos_y_floor;
      pos_z_rel = H * pos_z_floor;

      //* >> Computing the weights of the nearest grid points of the particle *//
      w_x_1 = GL_ptcl_x[ptcl_idx] - pos_x_rel;
      w_y_1 = GL_ptcl_y[ptcl_idx] - pos_y_rel;
      w_z_1 = GL_ptcl_z[ptcl_idx] - pos_z_rel;
      w_x_2 = H - w_x_1;
      w_y_2 = H - w_y_1;
      w_z_2 = H - w_z_1;
      w[0] = w_x_2 * w_y_2 * w_z_2;
      w[1] = w_x_1 * w_y_2 * w_z_2;
      w[2] = w_x_2 * w_y_1 * w_z_2;
      w[3] = w_x_1 * w_y_1 * w_z_2;
      w[4] = w_x_2 * w_y_2 * w_z_1;
      w[5] = w_x_1 * w_y_2 * w_z_1;
      w[6] = w_x_2 * w_y_1 * w_z_1;
      w[7] = w_x_1 * w_y_1 * w_z_1;

      //* >> Particle density contributes to 8 enclosure grid points *//
      for (int kk = 0; kk < 2; kk++)
      {
        for (int jj = 0; jj < 2; jj++)
        {
          for (int ii = 0; ii < 2; ii++)
          {
            aux_int = ii + 2 * jj + 4 * kk;
            box_grid_idxNbr[aux_int] = box_grid_idx + ii + jj * grid_box_real_dim_X + kk * grid_box_real_dim_X_times_Y;
            ptr_node->ptr_d[box_grid_idxNbr[aux_int]] += poisson_coeff_norm * GL_ptcl_mass[ptcl_idx] * w[aux_int];
          }
        }
      }

      for (int j = 1; j < cell_ptcl; j++)
      {
        ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

        //* >> Computing the weights of the nearest grid points of the particle *//
        w_x_1 = GL_ptcl_x[ptcl_idx] - pos_x_rel;
        w_y_1 = GL_ptcl_y[ptcl_idx] - pos_y_rel;
        w_z_1 = GL_ptcl_z[ptcl_idx] - pos_z_rel;
        w_x_2 = H - w_x_1;
        w_y_2 = H - w_y_1;
        w_z_2 = H - w_z_1;
        w[0] = w_x_2 * w_y_2 * w_z_2;
        w[1] = w_x_1 * w_y_2 * w_z_2;
        w[2] = w_x_2 * w_y_1 * w_z_2;
        w[3] = w_x_1 * w_y_1 * w_z_2;
        w[4] = w_x_2 * w_y_2 * w_z_1;
        w[5] = w_x_1 * w_y_2 * w_z_1;
        w[6] = w_x_2 * w_y_1 * w_z_1;
        w[7] = w_x_1 * w_y_1 * w_z_1;

        //* >> Particle density contributes to 8 enclosure grid points *//
        for (int k = 0; k < 8; k++)
        {
          ptr_node->ptr_d[box_grid_idxNbr[k]] += poisson_coeff_norm * GL_ptcl_mass[ptcl_idx] * w[k];
        }
      }
      counter_ptcl += cell_ptcl;
    }
  }

  // for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
  // {
  //     box_idx = ptr_node->ptr_box_idx[cell_idx];

  //     if (ptr_node->ptr_cell_struct[box_idx].ptcl_size > 0)
  //     {
  //         ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[0];

  //         //* >> Position of the particles in the grid level *//
  //         pos_x = GL_ptcl_x[ptcl_idx] * one_over_H;
  //         pos_y = GL_ptcl_y[ptcl_idx] * one_over_H;
  //         pos_z = GL_ptcl_z[ptcl_idx] * one_over_H;

  //         // //* >> floor of the particles positions in the grid level *//
  //         // pos_x_floor = myfloor(pos_x);
  //         // pos_y_floor = myfloor(pos_y);
  //         // pos_z_floor = myfloor(pos_z);
  //         pos_x_floor = (int)pos_x;
  //         pos_y_floor = (int)pos_y;
  //         pos_z_floor = (int)pos_z;

  //         box_grid_idx_x = pos_x_floor - ptr_node->box_ts_x;
  //         box_grid_idx_y = pos_y_floor - ptr_node->box_ts_y;
  //         box_grid_idx_z = pos_z_floor - ptr_node->box_ts_z;

  //         if (ptr_node->pbc_crosses_the_boundary_simulation_box == true)
  //         {
  //             if (pos_x_floor > ptr_node->box_max_x)
  //             {
  //                 box_grid_idx_x -= (1 << lv);
  //             }

  //             if (pos_y_floor > ptr_node->box_max_y)
  //             {
  //                 box_grid_idx_y -= (1 << lv);
  //             }

  //             if (pos_z_floor > ptr_node->box_max_z)
  //             {
  //                 box_grid_idx_z -= (1 << lv);
  //             }
  //         }

  //         box_grid_idx = box_grid_idx_x + box_grid_idx_y * grid_box_real_dim_X + box_grid_idx_z * grid_box_real_dim_X_times_Y;

  //         pos_x_rel = H * pos_x_floor;
  //         pos_y_rel = H * pos_y_floor;
  //         pos_z_rel = H * pos_z_floor;

  //         //* >> Computing the weights of the nearest grid points of the particle *//
  //         w_x_1 = GL_ptcl_x[ptcl_idx] - pos_x_rel;
  //         w_y_1 = GL_ptcl_y[ptcl_idx] - pos_y_rel;
  //         w_z_1 = GL_ptcl_z[ptcl_idx] - pos_z_rel;
  //         w_x_2 = H - w_x_1;
  //         w_y_2 = H - w_y_1;
  //         w_z_2 = H - w_z_1;
  //         w[0] = w_x_2 * w_y_2 * w_z_2;
  //         w[1] = w_x_1 * w_y_2 * w_z_2;
  //         w[2] = w_x_2 * w_y_1 * w_z_2;
  //         w[3] = w_x_1 * w_y_1 * w_z_2;
  //         w[4] = w_x_2 * w_y_2 * w_z_1;
  //         w[5] = w_x_1 * w_y_2 * w_z_1;
  //         w[6] = w_x_2 * w_y_1 * w_z_1;
  //         w[7] = w_x_1 * w_y_1 * w_z_1;

  //         //* >> Particle density contributes to 8 enclosure grid points *//
  //         for (int kk = 0; kk < 2; kk++)
  //         {
  //             for (int jj = 0; jj < 2; jj++)
  //             {
  //                 for (int ii = 0; ii < 2; ii++)
  //                 {
  //                     aux_int = ii + 2 * jj + 4 * kk;
  //                     box_grid_idxNbr[aux_int] = box_grid_idx + ii + jj * grid_box_real_dim_X + kk * grid_box_real_dim_X_times_Y;
  //                     ptr_node->ptr_d[box_grid_idxNbr[aux_int]] += poisson_coeff_norm * GL_ptcl_mass[ptcl_idx] * w[aux_int];
  //                 }
  //             }
  //         }
  //     }

  //     for (int j = 1; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
  //     {
  //         ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

  //         //* >> Computing the weights of the nearest grid points of the particle *//
  //         w_x_1 = GL_ptcl_x[ptcl_idx] - pos_x_rel;
  //         w_y_1 = GL_ptcl_y[ptcl_idx] - pos_y_rel;
  //         w_z_1 = GL_ptcl_z[ptcl_idx] - pos_z_rel;
  //         w_x_2 = H - w_x_1;
  //         w_y_2 = H - w_y_1;
  //         w_z_2 = H - w_z_1;
  //         w[0] = w_x_2 * w_y_2 * w_z_2;
  //         w[1] = w_x_1 * w_y_2 * w_z_2;
  //         w[2] = w_x_2 * w_y_1 * w_z_2;
  //         w[3] = w_x_1 * w_y_1 * w_z_2;
  //         w[4] = w_x_2 * w_y_2 * w_z_1;
  //         w[5] = w_x_1 * w_y_2 * w_z_1;
  //         w[6] = w_x_2 * w_y_1 * w_z_1;
  //         w[7] = w_x_1 * w_y_1 * w_z_1;

  //         //* >> Particle density contributes to 8 enclosure grid points *//
  //         for (int k = 0; k < 8; k++)
  //         {
  //             ptr_node->ptr_d[box_grid_idxNbr[k]] += poisson_coeff_norm * GL_ptcl_mass[ptcl_idx] * w[k];
  //         }
  //     }
  // }
}

void grid_density(void)
{

  //* >> Density in the grid *//

  // Head
  computing_grid_density_head(GL_tentacles[0][0]);

  // Branches
  for (int lv = 1; lv < GL_tentacles_level_max + 1; lv++)
  {
    // number of parent of the level = GL_tentacles_size[lv];
    //* >> For cycle over parent nodes *//
    for (int i = 0; i < GL_tentacles_size[lv]; i++)
    {
      // ptr_node = GL_tentacles[lv][i];
      computing_grid_density_branch(GL_tentacles[lv][i]);
    }
  }
}