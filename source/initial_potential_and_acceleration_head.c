/*
 * initial_potential_and_acceleration_head.c
 *
 * Compute the potential and acceleartion of the head node grid points in the
 * boundary of the simulation
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

#include "initial_potential_and_acceleration_head.h"

void initial_potential_and_acceleration_head(struct node *ptr_head)
{
  int grid_side;  // Number of grid points per side
  int grid_limit; // Superior border in the internal grid box

  vtype dist; // Distance between the center of mass and the border

  vtype H; // Size of the grid side

  int box_grid_idx; // Grid box index

  vtype aux_i, aux_j, aux_k;

  vtype aux_coeff; // Auxiliar coefficient

  grid_side = box_side_lmin + 1;
  

  grid_limit = grid_side - bder_os_sim - 1;

  // printf("grid_size = %d\n",grid_side);
  // printf(" bder_os_sim = %d, grid_limit + 1 = %d\n",bder_os_sim,grid_limit  + 1) ;

  H = 1.0 / no_lmin_cell;

  //* >> Computing the center of mass *//
  // GL_cm[0] = 0; // X position
  // GL_cm[1] = 0; // Y position
  // GL_cm[2] = 0; // Z position

  // for (int i = 0; i < GL_no_ptcl_final; i++)
  // {
  //   GL_cm[0] += GL_ptcl_mass[i] * GL_ptcl_x[i];
  //   GL_cm[1] += GL_ptcl_mass[i] * GL_ptcl_y[i];
  //   GL_cm[2] += GL_ptcl_mass[i] * GL_ptcl_z[i];
  // }

  // //* Normalizing per the total mass *//
  // GL_cm[0] = GL_cm[0] / GL_total_mass_initial;
  // GL_cm[1] = GL_cm[1] / GL_total_mass_initial;
  // GL_cm[2] = GL_cm[2] / GL_total_mass_initial;

  // printf(" cm_x = %f, cm_y = %f, cm_z =  %f\n",(double)GL_cm_x,(double)GL_cm_y,(double)GL_cm_z);
  // printf(" cm_vx = %f, cm_vy = %f, cm_vz =  %f\n",(double)GL_cm_vx,(double)GL_cm_vy,(double)GL_cm_vz);

  int border_grid_points = 0;
  int interior_grid_points = 0;

  if (bdry_cond_type == 0)
  {
    vtype dis_periodic[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    vtype aux_ijk = (grid_limit - bder_os_sim) * H - 0.5;
    int counter = 0;
    int box_grid_idx_periodic[8] = {0, 0, 0, 0, 0, 0, 0, 0};

    for (int k = bder_os_sim; k < grid_limit; k++)
    {
      aux_k = (k - bder_os_sim) * H - 0.5;
      for (int j = bder_os_sim; j < grid_limit; j++)
      {
        aux_j = (j - bder_os_sim) * H - 0.5;
        for (int i = bder_os_sim; i < grid_limit; i++)
        {
          aux_i = (i - bder_os_sim) * H - 0.5;

          if (i == bder_os_sim)
          {
            // Adding the box grid indices
            box_grid_idx_periodic[0] = i + j * grid_side + k * grid_side * grid_side;
            box_grid_idx_periodic[1] = grid_limit + j * grid_side + k * grid_side * grid_side;

            // Computing distance between the center of mass and the border point
            // i = bder_os_sim
            aux_coeff = aux_i * aux_i + aux_j * aux_j + aux_k* aux_k;
            dis_periodic[0] = mysqrt(aux_coeff);
            // i = grid_limit
            aux_coeff = aux_ijk * aux_ijk + aux_j * aux_j + aux_k * aux_k;
            dis_periodic[1] = mysqrt(aux_coeff);

            // Storing the minimum distance
            dist = dis_periodic[0] < dis_periodic[1] ? dis_periodic[0] : dis_periodic[1];

            counter++;
          }
          if (j == bder_os_sim)
          {
            // Adding the box grid indices
            if (counter == 0)
            {
              box_grid_idx_periodic[0] = i + j * grid_side + k * grid_side * grid_side;
              box_grid_idx_periodic[1] = 1 + grid_limit * grid_side + k * grid_side * grid_side;
              // Computing distance between the center of mass and the border point
              // j = bder_os_sim
              aux_coeff = aux_i * aux_i + aux_j * aux_j + aux_k * aux_k;
              dis_periodic[0] = mysqrt(aux_coeff);
              // j = grid_limit
              aux_coeff = aux_i* aux_i + aux_ijk* aux_ijk + aux_k* aux_k;
              dis_periodic[1] = mysqrt(aux_coeff);

              // Storing the minimum distance
              dist = dis_periodic[0] < dis_periodic[1] ? dis_periodic[0] : dis_periodic[1];
            }
            else
            {
              box_grid_idx_periodic[2] = 1 + grid_limit * grid_side + k * grid_side * grid_side;
              box_grid_idx_periodic[3] = grid_limit + grid_limit * grid_side + k * grid_side * grid_side;
              // Computing distance between the center of mass and the border point
              // i = bder_os_sim && j = grid_limit
              aux_coeff = aux_i* aux_i+ aux_ijk* aux_ijk+ aux_k * aux_k;
              dis_periodic[2] = mysqrt(aux_coeff);
              // i = grid_limit && j = grid_limit
              aux_coeff = aux_ijk* aux_ijk + aux_ijk* aux_ijk + aux_k* aux_k;
              dis_periodic[3] = mysqrt(aux_coeff);

              // Storing the minimum distance
              dist = dist < dis_periodic[2] ? dist : dis_periodic[2];
              dist = dist < dis_periodic[3] ? dist : dis_periodic[3];
            }
            counter++;
          }
          if (k == bder_os_sim)
          {
            // Adding the box grid indices
            if (counter == 0)
            {
              box_grid_idx_periodic[0] = i + j * grid_side + k * grid_side * grid_side;
              box_grid_idx_periodic[1] = 1 + j * grid_side + grid_limit * grid_side * grid_side;
              // Computing distance between the center of mass and the border point
              // k = bder_os_sim
              aux_coeff = aux_i * aux_i + aux_j * aux_j + aux_k * aux_k;
              dis_periodic[0] = mysqrt(aux_coeff);
              // k = grid_limit
              aux_coeff = aux_i * aux_i + aux_j* aux_j + aux_ijk * aux_ijk;
              dis_periodic[1] = mysqrt(aux_coeff);

              // Storing the minimum distance
              dist = dis_periodic[0] < dis_periodic[1] ? dis_periodic[0] : dis_periodic[1];
            }
            else if (counter == 1)
            {
              if (i == bder_os_sim)
              {
                box_grid_idx_periodic[2] = 1 + j * grid_side + grid_limit * grid_side * grid_side;
                box_grid_idx_periodic[3] = grid_limit + j * grid_side + grid_limit * grid_side * grid_side;
                // Computing distance between the center of mass and the border point
                // i = bder_os_sim && k = grid_limit
                aux_coeff = aux_i* aux_i+ aux_j * aux_j+ aux_ijk * aux_ijk;
                dis_periodic[2] = mysqrt(aux_coeff);
                // i = grid_limit && k = grid_limit
                aux_coeff = aux_ijk * aux_ijk + aux_j * aux_j + aux_ijk * aux_ijk;
                dis_periodic[3] = mysqrt(aux_coeff);
              }
              else
              {
                box_grid_idx_periodic[2] = 1 + j * grid_side + grid_limit * grid_side * grid_side;
                box_grid_idx_periodic[3] = i + grid_limit * grid_side + grid_limit * grid_side * grid_side;
                // Computing distance between the center of mass and the border point
                // j = bder_os_sim && k = grid_limit
                aux_coeff = aux_i * aux_i + aux_j * aux_j + aux_ijk * aux_ijk;
                dis_periodic[2] = mysqrt(aux_coeff);
                // j = grid_limit && k = grid_limit
                aux_coeff = aux_i * aux_i + aux_ijk * aux_ijk + aux_ijk * aux_ijk;
                dis_periodic[3] = mysqrt(aux_coeff);
              }

              // Storing the minimum distance
              dist = dist < dis_periodic[2] ? dist : dis_periodic[2];
              dist = dist < dis_periodic[3] ? dist : dis_periodic[3];
            }
            else if (counter == 2)
            {
              box_grid_idx_periodic[4] = 1 + j * grid_side + grid_limit * grid_side * grid_side;
              box_grid_idx_periodic[5] = grid_limit + j * grid_side + grid_limit * grid_side * grid_side;
              box_grid_idx_periodic[6] = 1 + grid_limit * grid_side + grid_limit * grid_side * grid_side;
              box_grid_idx_periodic[7] = grid_limit + grid_limit * grid_side + grid_limit * grid_side * grid_side;
              // Computing distance between the center of mass and the border point
              // i = bder_os_sim && j = bder_os_sim && k = grid_limit
              aux_coeff = aux_i * aux_i + aux_j * aux_j + aux_ijk * aux_ijk;
              dis_periodic[4] = mysqrt(aux_coeff);
              // i = grid_limit && j = bder_os_sim && k = grid_limit
              aux_coeff = aux_ijk * aux_ijk + aux_j * aux_j + aux_ijk * aux_ijk;
              dis_periodic[5] = mysqrt(aux_coeff);
              // i = bder_os_sim && j = grid_limit && k = grid_limit
              aux_coeff = aux_i * aux_i + aux_ijk * aux_ijk + aux_ijk * aux_ijk;
              dis_periodic[6] = mysqrt(aux_coeff);
              // i = grid_limit && j = grid_limit && k = grid_limit
              aux_coeff = aux_ijk * aux_ijk + aux_ijk * aux_ijk + aux_ijk * aux_ijk;
              dis_periodic[7] = mysqrt(aux_coeff);

              // Storing the minimum distance
              dist = dist < dis_periodic[4] ? dist : dis_periodic[4];
              dist = dist < dis_periodic[5] ? dist : dis_periodic[5];
              dist = dist < dis_periodic[6] ? dist : dis_periodic[6];
              dist = dist < dis_periodic[7] ? dist : dis_periodic[7];
            }

            counter++;
          }

          if (counter > 0)
          {
            aux_coeff = -_G_ * GL_total_mass_initial / dist;
            for (int l = 0; l < (1 << counter); l++)
            {
              ptr_head->ptr_pot[box_grid_idx_periodic[l]] = aux_coeff;
            }

            aux_coeff = -_G_ * GL_total_mass_initial / (dist * dist * dist);

            if (counter == 1)
            {
              //* Acceleration
              ptr_head->ptr_ax[box_grid_idx_periodic[0]] = aux_coeff * aux_i;
              ptr_head->ptr_ay[box_grid_idx_periodic[0]] = aux_coeff * aux_j;
              ptr_head->ptr_az[box_grid_idx_periodic[0]] = aux_coeff * aux_k;
              
            }
            else if (counter == 2)
            {
              //* Acceleration
              if (i == bder_os_sim)
              {
                ptr_head->ptr_ax[box_grid_idx_periodic[0]] = aux_coeff * aux_i;
                if (j == bder_os_sim)
                {
                  ptr_head->ptr_ay[box_grid_idx_periodic[0]] = aux_coeff * aux_j;
                  ptr_head->ptr_az[box_grid_idx_periodic[0]] = aux_coeff *  aux_k;
                }
                else
                {
                  ptr_head->ptr_ay[box_grid_idx_periodic[0]] = aux_coeff * aux_j;
                  ptr_head->ptr_az[box_grid_idx_periodic[0]] = aux_coeff * aux_k;
                }
              }
              else
              {
                ptr_head->ptr_ax[box_grid_idx_periodic[0]] = aux_coeff * aux_i;
                ptr_head->ptr_ay[box_grid_idx_periodic[0]] = aux_coeff * aux_j;
                ptr_head->ptr_az[box_grid_idx_periodic[0]] = aux_coeff * aux_k;
              }
            }
            else if (counter == 3)
            {
              //* Acceleration
              ptr_head->ptr_ax[box_grid_idx_periodic[0]] = aux_coeff * aux_i;
              ptr_head->ptr_ay[box_grid_idx_periodic[0]] = aux_coeff * aux_j;
              ptr_head->ptr_az[box_grid_idx_periodic[0]] = aux_coeff * aux_k;
            }
            for (int l = 1; l < (1 << counter); l++)
            {
              ptr_head->ptr_ax[box_grid_idx_periodic[l]] = ptr_head->ptr_ax[box_grid_idx_periodic[0]];
              ptr_head->ptr_ay[box_grid_idx_periodic[l]] = ptr_head->ptr_ay[box_grid_idx_periodic[0]];
              ptr_head->ptr_az[box_grid_idx_periodic[l]] = ptr_head->ptr_az[box_grid_idx_periodic[0]];
            }

            border_grid_points += (1 << counter);

            counter = 0;
          }
          else
          {
            interior_grid_points++;
          }
        }
      }
    }
  }
  else
  {

    for (int k = 0; k < box_side_lmin; k++)
    {
      aux_k = (k - bder_os_sim) * H - 0.5;
      //aux_k = (k - bder_os_sim) * H - GL_cm_z;
      for (int j = 0; j < box_side_lmin; j++)
      {
        aux_j = (j - bder_os_sim) * H - 0.5;
        //aux_j = (j - bder_os_sim) * H - GL_cm_y;
        for (int i = 0; i < box_side_lmin; i++)
        {
          aux_i = (i - bder_os_sim) * H - 0.5;
          //aux_i = (i - bder_os_sim) * H - GL_cm_x;
          dist = aux_i * aux_i + aux_j * aux_j + aux_k * aux_k;
          dist = mysqrt(dist);
          if(dist == 0.0)
          {
            dist = BoxSize/(1 << lmin) * 1.0e-3; //Putting an epsilon value
          }
          box_grid_idx = i + j * grid_side + k * grid_side * grid_side;
          aux_coeff = -_G_ * GL_total_mass_initial / (dist * dist * dist);
          aux_coeff = -_G_ * GL_total_mass_initial / (dist * dist * dist);
          ptr_head->ptr_pot[box_grid_idx] = -_G_ * GL_total_mass_initial / dist; 
          ptr_head->ptr_ax[box_grid_idx] = aux_coeff * aux_i;
          ptr_head->ptr_ay[box_grid_idx] = aux_coeff * aux_j;
          ptr_head->ptr_az[box_grid_idx] = aux_coeff * aux_k;

          // ptr_head->ptr_pot[box_grid_idx] = 0;
          // ptr_head->ptr_ax[box_grid_idx] =  0;
          // ptr_head->ptr_ay[box_grid_idx] =  0;
          // ptr_head->ptr_az[box_grid_idx] =  0;


          //printf("box_grid_idx = %d, i j k, ax ay az, pot = %d, %d, %d, %.12f, %.12f, %.12f, %.12f\n",box_grid_idx,i,j,k,(double)ptr_head->ptr_ax[box_grid_idx],(double)ptr_head->ptr_ay[box_grid_idx],(double)ptr_head->ptr_az[box_grid_idx],(double)ptr_head->ptr_pot[box_grid_idx] );
        }
      }
    }



    // for (int k = bder_os_sim; k < grid_limit + 1; k++)
    // {
    //   aux_k = (k - bder_os_sim) * H - 0.5L;
    //   for (int j = bder_os_sim; j < grid_limit + 1; j++)
    //   {
    //     aux_j = (j - bder_os_sim) * H - 0.5L;
    //     for (int i = bder_os_sim; i < grid_limit + 1; i++)
    //     {

    //       if (i == bder_os_sim || i == grid_limit || j == bder_os_sim || j == grid_limit || k == bder_os_sim || k == grid_limit)
    //       {
    //         aux_i = (i - bder_os_sim) * H - 0.5L;
    //         //printf("aux_ %d = %f\n",i,aux_i);
    //         // Computing distance between the center of mass and the border point
    //         dist = aux_i * aux_i + aux_j * aux_j + aux_k * aux_k;
    //         dist = mysqrt(dist);
    //         box_grid_idx = i + j * grid_side + k * grid_side * grid_side;
    //         aux_coeff = -_G_ * GL_total_mass_initial / (dist * dist * dist);
    //         ptr_head->ptr_pot[box_grid_idx] = -_G_ * GL_total_mass_initial / dist; 
    //         ptr_head->ptr_ax[box_grid_idx] =  aux_coeff * aux_i;
    //         ptr_head->ptr_ay[box_grid_idx] =  aux_coeff * aux_j;
    //         ptr_head->ptr_az[box_grid_idx] =  aux_coeff * aux_k;
    //         //printf(" i j k, ax ay az, pot = %d, %d, %d, %f, %f, %f, %f\n",i,j,k,(double)ptr_head->ptr_ax[box_grid_idx],(double)ptr_head->ptr_ay[box_grid_idx],(double)ptr_head->ptr_az[box_grid_idx],(double)ptr_head->ptr_pot[box_grid_idx] );
    //         // ptr_head->ptr_pot[box_grid_idx] = 0;
    //         // ptr_head->ptr_ax[box_grid_idx] = 0;
    //         // ptr_head->ptr_ay[box_grid_idx] = 0;
    //         // ptr_head->ptr_az[box_grid_idx] = 0;
    //       }
    //     }
    //   }
    // }




  }
}
