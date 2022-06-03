/*
 * observables.c
 *
 * Compute the energies of the simulation in a time-step given
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

#include "observables.h"
#include "space_check.h"

void kinetic_energy(vtype *energies)
{

    vtype _one_over_User_BoxSize = 1.0L / _User_BoxSize_;

    vtype particle_v_pow_2; // Particle velocity to the power of 2

    for (int i = 0; i < GL_no_ptcl; i++)
    {
        // Kinetic Energy
        particle_v_pow_2 = GL_ptcl_vx[i] * GL_ptcl_vx[i] + GL_ptcl_vy[i] * GL_ptcl_vy[i] + GL_ptcl_vz[i] * GL_ptcl_vz[i];
        // Adding Kinetic energy of the particle
        energies[0] += GL_ptcl_mass[i] * particle_v_pow_2;
    }

    energies[0] = energies[0] * 0.5 * _one_over_User_BoxSize;
}

void potential_energy_v1(vtype *energies)
{
    vtype distance_x, distance_y, distance_z; // Axial distance between 2 particles
    vtype distance;                           // Distance between 2 particles
    vtype _one_over_User_BoxSize = 1.0L / _User_BoxSize_;

    for (int i = 0; i < GL_no_ptcl; i++)
    {
        for (int j = 0; j < i; j++)
        {
            // Potential energy
            distance_x = GL_ptcl_x[i] - GL_ptcl_x[j];
            distance_y = GL_ptcl_y[i] - GL_ptcl_y[j];
            distance_z = GL_ptcl_z[i] - GL_ptcl_z[j];
            distance = distance_x * distance_x + distance_y * distance_y + distance_z * distance_z;
            distance = sqrt(distance);
            // Checking minimal distance accepted
            if (distance < _Min_Particle_Distance_For_Energy_Computation_)
            {
                distance = _Min_Particle_Distance_For_Energy_Computation_;
            }
            // Adding Potential energy of the particle
            energies[1] += GL_ptcl_mass[i] * GL_ptcl_mass[j] / (distance);
        }
    }

    energies[1] = -_G_ * energies[1] * _one_over_User_BoxSize;
}

void computing_partial_d(struct node *ptr_node, vtype *partial_d)
{
    int box_grid_idx_x; // grid_idx component in the grid box
    int box_grid_idx_y;
    int box_grid_idx_z;
    int box_grid_idx; // Grid box grid_idx

    int box_idx; // Box index of the node cell

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

    int grid_box_real_dim_X = (ptr_node->box_real_dim_x + 1);
    int grid_box_real_dim_X_times_Y = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

    int box_grid_idxNbr[8];
    int aux_int;

    if(ptr_node->chn_size > 0)
    {
        // Initial density:
        for (int i = 0; i < grid_box_real_dim_X_times_Y * (ptr_node->box_real_dim_z + 1); i++)
        {
            partial_d[i] = 0;
        }

        for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
        {
            box_idx = ptr_node->ptr_box_idx[cell_idx];

            if (ptr_node->ptr_box[box_idx] < 0)
            {
                if (ptr_node->ptr_cell_struct[box_idx].ptcl_size > 0)
                {
                    ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[0];

                    //** >> Position of the particles in the grid level **/
                    pos_x = GL_ptcl_x[ptcl_idx] * one_over_H;
                    pos_y = GL_ptcl_y[ptcl_idx] * one_over_H;
                    pos_z = GL_ptcl_z[ptcl_idx] * one_over_H;

                    // ** >> floor of the particles positions in the grid level **/
                    pos_x_floor = (int)pos_x;
                    pos_y_floor = (int)pos_y;
                    pos_z_floor = (int)pos_z;

                    box_grid_idx_x = pos_x_floor - ptr_node->box_ts_x;
                    box_grid_idx_y = pos_y_floor - ptr_node->box_ts_y;
                    box_grid_idx_z = pos_z_floor - ptr_node->box_ts_z;
                    box_grid_idx = box_grid_idx_x + box_grid_idx_y * grid_box_real_dim_X + box_grid_idx_z * grid_box_real_dim_X_times_Y;

                    pos_x_rel = H * pos_x_floor;
                    pos_y_rel = H * pos_y_floor;
                    pos_z_rel = H * pos_z_floor;

                    //** >> Computing the weights of the nearest grid points of the particle **/
                    w_x_1 = GL_ptcl_x[ptcl_idx] - pos_x_rel;
                    w_y_1 = GL_ptcl_y[ptcl_idx] - pos_y_rel;
                    w_z_1 = GL_ptcl_z[ptcl_idx] - pos_z_rel;
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

                    //** >> Particle density contributes to 8 enclosure grid points **/
                    for (int kk = 0; kk < 2; kk++)
                    {
                        for (int jj = 0; jj < 2; jj++)
                        {
                            for (int ii = 0; ii < 2; ii++)
                            {
                                aux_int = ii + 2 * jj + 4 * kk;
                                box_grid_idxNbr[aux_int] = box_grid_idx + ii + jj * grid_box_real_dim_X + kk * grid_box_real_dim_X_times_Y;
                                partial_d[box_grid_idxNbr[aux_int]] +=  GL_ptcl_mass[ptcl_idx] * w[aux_int];
                            }
                        }
                    }
                }

                for (int j = 1; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
                {
                    ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

                    //** >> Computing the weights of the nearest grid points of the particle **/
                    w_x_1 = GL_ptcl_x[ptcl_idx] - pos_x_rel;
                    w_y_1 = GL_ptcl_y[ptcl_idx] - pos_y_rel;
                    w_z_1 = GL_ptcl_z[ptcl_idx] - pos_z_rel;
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

                    //** >> Particle density contributes to 8 enclosure grid points **/
                    for (int k = 0; k < 8; k++)
                    {
                        partial_d[box_grid_idxNbr[k]] += GL_ptcl_mass[ptcl_idx] * w[k];
                    }
                }
            }
        }
    }
    else
    {
        vtype aux_C = H * H * H / (4 * _PI_ * _G_);
        for (int i = 0; i < grid_box_real_dim_X_times_Y * (ptr_node->box_real_dim_z + 1); i++)
        {
            partial_d[i] = ptr_node->ptr_d[i] * aux_C;
        }
    }

}

int potential_energy_v2(vtype *energies)
{
    vtype _one_over_User_BoxSize = 1.0L / _User_BoxSize_;
    struct node *ptr_node;

    int no_pts;

    int box_grid_idx;

    vtype *partial_d = NULL;
    int size;
    int partial_d_cap = 0;

    for (int lv = GL_tentacles_level_max; lv > -1; lv-- )
    {
        no_pts = GL_tentacles_size[lv];
        for(int i = 0; i< no_pts; i++)
        {
            ptr_node = GL_tentacles[lv][i];

            //size = ptr_node->grid_properties_cap;
            size = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1) * (ptr_node->box_real_dim_z + 1);
            if (space_check(&(partial_d_cap), size, 2.0f, "p1v1", &(partial_d)) == _FAILURE_)
            {
                printf("Error, in space_check function\n");
                return _FAILURE_;
            }

            computing_partial_d(ptr_node, partial_d);
            for (int j = 0; j < ptr_node->grid_intr_size; j++)
            {
                box_grid_idx = ptr_node->ptr_intr_grid_idx[j];
                energies[1] += partial_d[box_grid_idx] * ptr_node->ptr_pot[box_grid_idx];
            }

            for (int j = 0; j < ptr_node->grid_bder_size; j++)
            {
                box_grid_idx = ptr_node->ptr_bder_grid_idx[j];
                energies[1] += partial_d[box_grid_idx] * ptr_node->ptr_pot[box_grid_idx];
            }
        }
    }

    energies[1] = energies[1] * 0.5 * _one_over_User_BoxSize;

    free(partial_d);

    return _SUCCESS_;
}

static void computing_particle_potential_head_only(const struct node *ptr_head, vtype *energies)
{
    int aux_idx;

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

    int box_grid_idx_x; // grid_idx component in the grid box
    int box_grid_idx_y;
    int box_grid_idx_z;
    int box_grid_idx;    // Grid box grid_idx
    int box_grid_idxNbr; // Box index in the neigborhood

    vtype w[8]; // Weight of the CIC method

    int lv = ptr_head->lv; // Level of refinement

    vtype aux_pot;

    int grid_box_real_dim_X = (ptr_head->box_real_dim_x + 1);
    int grid_box_real_dim_X_times_Y = (ptr_head->box_real_dim_x + 1) * (ptr_head->box_real_dim_y + 1);

    for (int i = 0; i < GL_no_ptcl; i++)
    {
        //** >> Position of the particles in the grid level **/
        pos_x = GL_ptcl_x[i] * (1 << lv);
        pos_y = GL_ptcl_y[i] * (1 << lv);
        pos_z = GL_ptcl_z[i] * (1 << lv);

        //** >> Floor of the particles positions in the grid level **/
        pos_x_floor = (int)pos_x;
        pos_y_floor = (int)pos_y;
        pos_z_floor = (int)pos_z;

        //** >> Computing the weights of the nearest grid points of the particle **/
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

        // //** >> Computing the weights of the nearest grid points of the particle **/
        // // Each for cyle yields 2 options: X or 1-X, where X =  pos_x - pos_x_floor
        // for (int kk = 0; kk < 2; kk++)
        // {
        //     w_z = kk + (1 - 2 * kk) * (pos_z - pos_z_floor);
        //     for (int jj = 0; jj < 2; jj++)
        //     {
        //         w_y = jj + (1 - 2 * jj) * (pos_y - pos_y_floor);
        //         for (int ii = 0; ii < 2; ii++)
        //         {
        //             w_x = ii + (1 - 2 * ii) * (pos_x - pos_x_floor);
        //             w[kk * 4 + jj * 2 + ii] = (1 - w_x) * (1 - w_y) * (1 - w_z);
        //         }
        //     }
        // }

        box_grid_idx_x = pos_x_floor - ptr_head->box_ts_x;
        box_grid_idx_y = pos_y_floor - ptr_head->box_ts_y;
        box_grid_idx_z = pos_z_floor - ptr_head->box_ts_z;
        box_grid_idx = box_grid_idx_x + box_grid_idx_y * grid_box_real_dim_X + box_grid_idx_z * grid_box_real_dim_X_times_Y;
        //** >> Particle density contributes to 8 enclosure grid points **/
        aux_pot = 0;
        for (int kk = 0; kk < 2; kk++)
        {
            for (int jj = 0; jj < 2; jj++)
            {
                for (int ii = 0; ii < 2; ii++)
                {
                    // box_grid_idx_x = (pos_x_floor + ii) - ptr_head->box_ts_x;
                    // box_grid_idx_y = (pos_y_floor + jj) - ptr_head->box_ts_y;
                    // box_grid_idx_z = (pos_z_floor + kk) - ptr_head->box_ts_z;
                    // box_grid_idx = box_grid_idx_x + box_grid_idx_y * (ptr_head->box_real_dim_x + 1) + box_grid_idx_z * (ptr_head->box_real_dim_x + 1) * (ptr_head->box_real_dim_y + 1);
                    aux_idx = ii + 2 * jj + 4 * kk;
                    box_grid_idxNbr = box_grid_idx + ii + jj * grid_box_real_dim_X + kk * grid_box_real_dim_X_times_Y;
                    aux_pot += ptr_head->ptr_pot[box_grid_idxNbr] * w[aux_idx];
                }
            }
        }
        //Energy computation:
        energies[1] += GL_ptcl_mass[i] * aux_pot;
    }
}

static void computing_particle_potential_head_plus_branches(const struct node *ptr_node, vtype *energies)
{
    int box_idx; // Box index of the node cell

    int ptcl_idx; // Particle grid_idx in the node

    int aux_idx;

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

    int box_grid_idx_x; // grid_idx component in the grid box
    int box_grid_idx_y;
    int box_grid_idx_z;
    int box_grid_idx;    // Grid box grid_idx
    int box_grid_idxNbr; // Box index in the neigborhood

    vtype aux_pot;

    vtype w[8]; // Weight of the CIC method

    int lv = ptr_node->lv; // Level of refinement

    int grid_box_real_dim_X = (ptr_node->box_real_dim_x + 1);
    int grid_box_real_dim_X_times_Y = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

    //** >> Case no more child, the node is a leaf **/
    if (ptr_node->chn_size == 0)
    {
        for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
        {
            box_idx = ptr_node->ptr_box_idx[cell_idx];
            for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
            {
                ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

                //** >> Position of the particles in the grid level **/
                pos_x = GL_ptcl_x[ptcl_idx] * (1 << lv);
                pos_y = GL_ptcl_y[ptcl_idx] * (1 << lv);
                pos_z = GL_ptcl_z[ptcl_idx] * (1 << lv);

                //** >> Floor of the particles positions in the grid level **/
                pos_x_floor = (int)pos_x;
                pos_y_floor = (int)pos_y;
                pos_z_floor = (int)pos_z;

                //** >> Computing the weights of the nearest grid points of the particle **/
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
                // Each for cyle yields 2 options: X or 1-X, where X =  pos_x - pos_x_floor
                // for (int kk = 0; kk < 2; kk++)
                // {
                //     w_z = kk + (1 - 2 * kk) * (pos_z - pos_z_floor);
                //     for (int jj = 0; jj < 2; jj++)
                //     {
                //         w_y = jj + (1 - 2 * jj) * (pos_y - pos_y_floor);
                //         for (int ii = 0; ii < 2; ii++)
                //         {
                //             w_x = ii + (1 - 2 * ii) * (pos_x - pos_x_floor);
                //             w[kk * 4 + jj * 2 + ii] = (1 - w_x) * (1 - w_y) * (1 - w_z);
                //         }
                //     }
                // }

                box_grid_idx_x = pos_x_floor - ptr_node->box_ts_x;
                box_grid_idx_y = pos_y_floor - ptr_node->box_ts_y;
                box_grid_idx_z = pos_z_floor - ptr_node->box_ts_z;
                box_grid_idx = box_grid_idx_x + box_grid_idx_y * grid_box_real_dim_X + box_grid_idx_z * grid_box_real_dim_X_times_Y;
                //** >> Particle density contributes to 8 enclosure grid points **/
                aux_pot = 0;
                for (int kk = 0; kk < 2; kk++)
                {
                    for (int jj = 0; jj < 2; jj++)
                    {
                        for (int ii = 0; ii < 2; ii++)
                        {
                            // box_grid_idx_x = (pos_x_floor + ii) - ptr_node->box_ts_x;
                            // box_grid_idx_y = (pos_y_floor + jj) - ptr_node->box_ts_y;
                            // box_grid_idx_z = (pos_z_floor + kk) - ptr_node->box_ts_z;
                            // box_grid_idx = box_grid_idx_x + box_grid_idx_y * (ptr_node->box_real_dim_x + 1) + box_grid_idx_z * (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);
                            aux_idx = ii + 2 * jj + 4 * kk;
                            box_grid_idxNbr = box_grid_idx + ii + jj * grid_box_real_dim_X + kk * grid_box_real_dim_X_times_Y;
                            aux_pot += ptr_node->ptr_pot[box_grid_idxNbr] * w[aux_idx];
                        }
                    }
                }

                // Energy computation:
                energies[1] += GL_ptcl_mass[ptcl_idx] * aux_pot;
            }
        }
    }
    //** >> Case there are more children, the node is a branch **/
    else
    {
        for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
        {
            box_idx = ptr_node->ptr_box_idx[cell_idx];

            if (ptr_node->ptr_box[box_idx] < 0)
            {
                for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
                {
                    ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

                    //** >> Position of the particles in the grid level **/
                    pos_x = GL_ptcl_x[ptcl_idx] * (1 << lv);
                    pos_y = GL_ptcl_y[ptcl_idx] * (1 << lv);
                    pos_z = GL_ptcl_z[ptcl_idx] * (1 << lv);

                    //** >> Floor of the particles positions in the grid level **/
                    pos_x_floor = (int)pos_x;
                    pos_y_floor = (int)pos_y;
                    pos_z_floor = (int)pos_z;

                    //** >> Computing the weights of the nearest grid points of the particle **/
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
                    // // Each for cyle yields 2 options: X or 1-X, where X =  pos_x - pos_x_floor
                    // for (int kk = 0; kk < 2; kk++)
                    // {
                    //     //w_z = kk + (1 - 2 * kk) * (pos_z - pos_z_floor);
                    //     for (int jj = 0; jj < 2; jj++)
                    //     {
                    //         //w_y = jj + (1 - 2 * jj) * (pos_y - pos_y_floor);
                    //         for (int ii = 0; ii < 2; ii++)
                    //         {
                    //             //w_x = ii + (1 - 2 * ii) * (pos_x - pos_x_floor);
                    //             w[kk * 4 + jj * 2 + ii] = (1 - w_x) * (1 - w_y) * (1 - w_z);
                    //         }
                    //     }
                    // }

                    box_grid_idx_x = pos_x_floor - ptr_node->box_ts_x;
                    box_grid_idx_y = pos_y_floor - ptr_node->box_ts_y;
                    box_grid_idx_z = pos_z_floor - ptr_node->box_ts_z;
                    box_grid_idx = box_grid_idx_x + box_grid_idx_y * grid_box_real_dim_X + box_grid_idx_z * grid_box_real_dim_X_times_Y;
                    //** >> Particle density contributes to 8 enclosure grid points **/
                    aux_pot = 0;
                    for (int kk = 0; kk < 2; kk++)
                    {
                        for (int jj = 0; jj < 2; jj++)
                        {
                            for (int ii = 0; ii < 2; ii++)
                            {
                                // box_grid_idx_x = (pos_x_floor + ii) - ptr_node->box_ts_x;
                                // box_grid_idx_y = (pos_y_floor + jj) - ptr_node->box_ts_y;
                                // box_grid_idx_z = (pos_z_floor + kk) - ptr_node->box_ts_z;
                                // box_grid_idx = box_grid_idx_x + box_grid_idx_y * (ptr_node->box_real_dim_x + 1) + box_grid_idx_z * (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);
                                aux_idx = ii + 2 * jj + 4 * kk;
                                box_grid_idxNbr = box_grid_idx + ii + jj * grid_box_real_dim_X + kk * grid_box_real_dim_X_times_Y;
                                aux_pot += ptr_node->ptr_pot[box_grid_idxNbr] * w[aux_idx];
                            }
                        }
                    }

                    // Energy computation:
                    energies[1] += GL_ptcl_mass[ptcl_idx] * aux_pot;
                }
            }
        }
    }
}

void potential_energy_v3(vtype *energies)
{
    vtype _one_over_User_BoxSize = 1.0L / _User_BoxSize_;

    if (lmin < lmax)
    {
        for (int lv = GL_tentacles_level_max; lv > -1; lv--)
        {
            //** >> For cycle over parent nodes **/
            for (int i = 0; i < GL_tentacles_size[lv]; i++)
            {
                // ptr_node = GL_tentacles[lv][i];

                computing_particle_potential_head_plus_branches(GL_tentacles[lv][i], energies);
            }
        }
    }
    else
    {
        computing_particle_potential_head_only(GL_tentacles[0][0], energies);
    }

    energies[1] = energies[1] * 0.5 * _one_over_User_BoxSize;
}

void observables(vtype *energies)
{

    energies[0] = 0; // Kinetic
    energies[1] = 0; // Potential
    energies[2] = 0; // Total

    kinetic_energy(energies);
    //potential_energy_v1(energies);
    // if(potential_energy_v2(energies) == _FAILURE_)
    // {
    //     printf("Error\n");
    // }
    potential_energy_v3(energies);

    energies[2] = energies[0] + energies[1];

}
