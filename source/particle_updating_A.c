/*
 * particle_updating_A.c
 *
 * Upadte global velocity and position particles by a "Predictor step"
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

#include "particle_updating_A.h"

static int computing_particles_updating_A(struct node *ptr_node, vtype dt, bool status)
{

    int no_ptcl; // Total number of particles in the node

    struct node *ptr_node_ch; // Child node of the node ptr_node
    struct node *ptr_node_pt; // parent node of the node ptr_node
    struct node *ptr_node_sib; // sibling node of the node ptr_node

    int box_idx_x_node; // Box index in X direcction of the node cell
    int box_idx_y_node; // Box index in Y direcction of the node cell
    int box_idx_z_node; // Box index in Z direcction of the node cell

    int ptcl_idx; // Particle grid_idx in the node

    int zone_idx; // Index of the refinement zone

    int box_idx_node_old;   // Box index of the current node
    int box_idx_node_new;   // Box index of the current node
    int box_idx_pt;   // Box index of the parent node
    int box_idx_sib;   // Box index of the sibling node
    int box_idx_ch;   // Box index of the child node


    if(lmin < lmax)
    {
        if (ptr_node->chn_size == 0)
        {
            for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
            {
                box_idx_x_node = ptr_node->ptr_cell_idx_x[cell_idx] - ptr_node->box_ts_x;
                box_idx_y_node = ptr_node->ptr_cell_idx_y[cell_idx] - ptr_node->box_ts_y;
                box_idx_z_node = ptr_node->ptr_cell_idx_z[cell_idx] - ptr_node->box_ts_z;
                box_idx_node_old = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

                no_ptcl = ptr_node->ptr_cell_struct[box_idx_node_old].ptcl_size;
                for (int j = 0; j < no_ptcl; j++)
                {
                    ptcl_idx = ptr_node->ptr_cell_struct[box_idx_node_old].ptr_ptcl[j];
                    if (GL_ptcl_updating_flag[ptcl_idx] != status)
                    {

                        //** >> Updating the new position of the particle **/
                        //** >> Velocities **/
                        GL_ptcl_vx[ptcl_idx] += GL_ptcl_ax[ptcl_idx] * dt / 2;
                        GL_ptcl_vy[ptcl_idx] += GL_ptcl_ay[ptcl_idx] * dt / 2;
                        GL_ptcl_vz[ptcl_idx] += GL_ptcl_az[ptcl_idx] * dt / 2;

                        //** >> Positions **/
                        GL_ptcl_x[ptcl_idx] += GL_ptcl_vx[ptcl_idx] * dt;
                        GL_ptcl_y[ptcl_idx] += GL_ptcl_vy[ptcl_idx] * dt;
                        GL_ptcl_z[ptcl_idx] += GL_ptcl_vz[ptcl_idx] * dt;

                        //** >> Checking if the particle exits the simulation **/
                        if (GL_ptcl_x[ptcl_idx] < 0 || GL_ptcl_x[ptcl_idx] > 1 || GL_ptcl_y[ptcl_idx] < 0 || GL_ptcl_y[ptcl_idx] > 1 || GL_ptcl_z[ptcl_idx] < 0 || GL_ptcl_z[ptcl_idx] > 1)
                        {
                            printf("Error, Partícula %d, sale de la simulación at positions:\n", ptcl_idx);
                            printf("x = %f\n", (double)GL_ptcl_x[ptcl_idx]);
                            printf("y = %f\n", (double)GL_ptcl_y[ptcl_idx]);
                            printf("z = %f\n", (double)GL_ptcl_z[ptcl_idx]);
                            printf("vx = %f\n", (double)GL_ptcl_vx[ptcl_idx]);
                            printf("vy = %f\n", (double)GL_ptcl_vy[ptcl_idx]);
                            printf("vz = %f\n", (double)GL_ptcl_vz[ptcl_idx]);
                            printf("ax = %f\n", (double)GL_ptcl_ax[ptcl_idx]);
                            printf("ay = %f\n", (double)GL_ptcl_ay[ptcl_idx]);
                            printf("az = %f\n", (double)GL_ptcl_az[ptcl_idx]);
                            printf("index = %d\n", ptcl_idx);

                            return _FAILURE_;
                        }

                        //** >> Moving the particle to the new node if it is necessary **/
                        box_idx_node_new = ptcl_idx_to_box_idx(ptr_node, ptcl_idx);

                        //** We ask if the particle leaves the node
                        //** >> The particle moves towards its parent node or towards some sibling node  **/

                        if (box_idx_node_old != box_idx_node_new)
                        {
                            if (ptr_node->ptr_box[box_idx_node_new] < -3)
                            {
                                ptr_node_pt = ptr_node->ptr_pt;
                                //** >> Box index in the parent node **/
                                box_idx_pt = ptcl_idx_to_box_idx(ptr_node_pt, ptcl_idx);

                                //** >> If the particle moves towards a sibling node **/
                                if (ptr_node_pt->ptr_box[box_idx_pt] >= 0)
                                {
                                    zone_idx = ptr_node_pt->ptr_box[box_idx_pt];
                                    ptr_node_sib = ptr_node_pt->pptr_chn[zone_idx];

                                    box_idx_sib = ptcl_idx_to_box_idx(ptr_node_sib, ptcl_idx);

                                    //** >> Space checking of the particle capacity in the sibling node **/
                                    if (space_check(&(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_cap), ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_size + 1, 2.0f, "p1i1", &(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptr_ptcl)) == _FAILURE_)
                                    {
                                        printf("Error, in space_check function\n");
                                        return _FAILURE_;
                                    }

                                    //** >> Adding the particle in the sibling node **/
                                    ptr_node_sib->ptr_cell_struct[box_idx_sib].ptr_ptcl[ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_size] = ptcl_idx;
                                    ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_size += 1;
                                    ptr_node_sib->ptr_cell_struct[box_idx_sib].cell_mass += GL_ptcl_mass[ptcl_idx]; // Cell mass
                                    ptr_node_sib->local_mass += GL_ptcl_mass[ptcl_idx];                             // Local mass
                                }
                                //** If the particle is only in the parent node **/
                                else
                                {
                                    //** >> Space checking of the particle capacity in the sibling cell node **/
                                    if (space_check(&(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_cap), ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_size + 1, 2.0f, "p1i1", &(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptr_ptcl)) == _FAILURE_)
                                    {
                                        printf("Error, in space_check function\n");
                                        return _FAILURE_;
                                    }
                                    //** >> Adding the particle in the parent cell node **/
                                    ptr_node_pt->ptr_cell_struct[box_idx_pt].ptr_ptcl[ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_size] = ptcl_idx;
                                    ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_size += 1;
                                    ptr_node_pt->ptr_cell_struct[box_idx_pt].cell_mass += GL_ptcl_mass[ptcl_idx]; // Cell mass
                                }

                                //** The local mass is reduced **/
                                ptr_node->local_mass -= GL_ptcl_mass[ptcl_idx];
                            }
                            //** >> The particle stay in the node **/
                            else
                            {
                                //** >> Space checking of the particle capacity in the sibling cell **/
                                if (space_check(&(ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_cap), ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_size + 1, 2.0f, "p1i1", &(ptr_node->ptr_cell_struct[box_idx_node_new].ptr_ptcl)) == _FAILURE_)
                                {
                                    printf("Error, in space_check function\n");
                                    return _FAILURE_;
                                }

                                //** >> Adding the particle in the sibling cell **/
                                ptr_node->ptr_cell_struct[box_idx_node_new].ptr_ptcl[ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_size] = ptcl_idx;
                                ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_size += 1;
                                ptr_node->ptr_cell_struct[box_idx_node_new].cell_mass += GL_ptcl_mass[ptcl_idx]; // Cell mass
                            }

                            // Whether the particle stays at the parent node or moves to a sibling node, it must be removed from the current cell of the old node
                            //** >> Removing the particle index of the current node ptr_node **/
                            // We move the last element of the array to the current position
                            ptr_node->ptr_cell_struct[box_idx_node_old].ptr_ptcl[j] = ptr_node->ptr_cell_struct[box_idx_node_old].ptr_ptcl[no_ptcl - 1];
                            no_ptcl--; // The total number of particle decrease
                            j--;       // The last element that was just moved to the current position should also must be analized
                            ptr_node->ptr_cell_struct[box_idx_node_old].cell_mass -= GL_ptcl_mass[ptcl_idx];
                        }
                        //** >> The status of the particle is changed from not updated to updated **/
                        GL_ptcl_updating_flag[ptcl_idx] = status;
                    }
                }
                ptr_node->ptr_cell_struct[box_idx_node_old].ptcl_size = no_ptcl;
            }
        } // End cycle over particles in the node
        else
        {
            for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
            {
                box_idx_x_node = ptr_node->ptr_cell_idx_x[cell_idx] - ptr_node->box_ts_x;
                box_idx_y_node = ptr_node->ptr_cell_idx_y[cell_idx] - ptr_node->box_ts_y;
                box_idx_z_node = ptr_node->ptr_cell_idx_z[cell_idx] - ptr_node->box_ts_z;
                box_idx_node_old = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

                no_ptcl = ptr_node->ptr_cell_struct[box_idx_node_old].ptcl_size;
                for (int j = 0; j < no_ptcl; j++)
                {
                    ptcl_idx = ptr_node->ptr_cell_struct[box_idx_node_old].ptr_ptcl[j];
                    if (GL_ptcl_updating_flag[ptcl_idx] != status)
                    {

                        //** >> Updating the new position of the particle **/
                        //** >> Velocities **/
                        GL_ptcl_vx[ptcl_idx] += GL_ptcl_ax[ptcl_idx] * dt / 2;
                        GL_ptcl_vy[ptcl_idx] += GL_ptcl_ay[ptcl_idx] * dt / 2;
                        GL_ptcl_vz[ptcl_idx] += GL_ptcl_az[ptcl_idx] * dt / 2;

                        //** >> Positions **/
                        GL_ptcl_x[ptcl_idx] += GL_ptcl_vx[ptcl_idx] * dt;
                        GL_ptcl_y[ptcl_idx] += GL_ptcl_vy[ptcl_idx] * dt;
                        GL_ptcl_z[ptcl_idx] += GL_ptcl_vz[ptcl_idx] * dt;

                        //** >> Checking if the particle exits the simulation **/
                        if (GL_ptcl_x[ptcl_idx] < 0 || GL_ptcl_x[ptcl_idx] > 1 || GL_ptcl_y[ptcl_idx] < 0 || GL_ptcl_y[ptcl_idx] > 1 || GL_ptcl_z[ptcl_idx] < 0 || GL_ptcl_z[ptcl_idx] > 1)
                        {
                            printf("Error, Partícula %d, sale de la simulación at positions:\n", ptcl_idx);
                            printf("x = %f\n", (double)GL_ptcl_x[ptcl_idx]);
                            printf("y = %f\n", (double)GL_ptcl_y[ptcl_idx]);
                            printf("z = %f\n", (double)GL_ptcl_z[ptcl_idx]);
                            printf("vx = %f\n", (double)GL_ptcl_vx[ptcl_idx]);
                            printf("vy = %f\n", (double)GL_ptcl_vy[ptcl_idx]);
                            printf("vz = %f\n", (double)GL_ptcl_vz[ptcl_idx]);
                            printf("ax = %f\n", (double)GL_ptcl_ax[ptcl_idx]);
                            printf("ay = %f\n", (double)GL_ptcl_ay[ptcl_idx]);
                            printf("az = %f\n", (double)GL_ptcl_az[ptcl_idx]);
                            printf("index = %d\n", ptcl_idx);

                            return _FAILURE_;
                        }

                        //** >> Moving the particle to the new node if it is necessary **/
                        box_idx_node_new = ptcl_idx_to_box_idx(ptr_node, ptcl_idx);

                        //** We ask if the particle leaves the node
                        //** >> The particle moves towards its parent node or towards some sibling node  **/

                        if (box_idx_node_old != box_idx_node_new)
                        {
                            //** >> The particle moves towards one of its child nodes **/
                            if (ptr_node->ptr_box[box_idx_node_new] >= 0)
                            {
                                zone_idx = ptr_node->ptr_box[box_idx_node_new];
                                ptr_node_ch = ptr_node->pptr_chn[zone_idx];

                                box_idx_ch = ptcl_idx_to_box_idx(ptr_node_ch, ptcl_idx);

                                //** >> Space checking of the particle capacity in the child node cell **/
                                if (space_check(&(ptr_node_ch->ptr_cell_struct[box_idx_ch].ptcl_cap), ptr_node_ch->ptr_cell_struct[box_idx_ch].ptcl_size + 1, 2.0f, "p1i1", &(ptr_node_ch->ptr_cell_struct[box_idx_ch].ptr_ptcl)) == _FAILURE_)
                                {
                                    printf("Error, in space_check function\n");
                                    return _FAILURE_;
                                }

                                //** >> Adding the particle in the child cell node **/
                                ptr_node_ch->ptr_cell_struct[box_idx_ch].ptr_ptcl[ptr_node_ch->ptr_cell_struct[box_idx_ch].ptcl_size] = ptcl_idx;
                                ptr_node_ch->ptr_cell_struct[box_idx_ch].ptcl_size += 1;
                                ptr_node_ch->ptr_cell_struct[box_idx_ch].cell_mass += GL_ptcl_mass[ptcl_idx]; // Cell mass
                                ptr_node_ch->local_mass += GL_ptcl_mass[ptcl_idx];
                            }
                            //** >> The particle moves towards its parent node or towards some sibling node  **/
                            else if (ptr_node->ptr_box[box_idx_node_new] < -3)
                            {
                                ptr_node_pt = ptr_node->ptr_pt;
                                //** >> Box index in the parent node **/
                                box_idx_pt = ptcl_idx_to_box_idx(ptr_node_pt, ptcl_idx);

                                //** >> If the particle moves towards a sibling node **/
                                if (ptr_node_pt->ptr_box[box_idx_pt] >= 0)
                                {
                                    zone_idx = ptr_node_pt->ptr_box[box_idx_pt];
                                    ptr_node_sib = ptr_node_pt->pptr_chn[zone_idx];

                                    box_idx_sib = ptcl_idx_to_box_idx(ptr_node_sib, ptcl_idx);

                                    //** >> Space checking of the particle capacity in the sibling node **/
                                    if (space_check(&(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_cap), ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_size + 1, 2.0f, "p1i1", &(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptr_ptcl)) == _FAILURE_)
                                    {
                                        printf("Error, in space_check function\n");
                                        return _FAILURE_;
                                    }

                                    //** >> Adding the particle in the sibling node **/
                                    ptr_node_sib->ptr_cell_struct[box_idx_sib].ptr_ptcl[ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_size] = ptcl_idx;
                                    ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_size += 1;
                                    ptr_node_sib->ptr_cell_struct[box_idx_sib].cell_mass += GL_ptcl_mass[ptcl_idx]; // Cell mass
                                    ptr_node_sib->local_mass += GL_ptcl_mass[ptcl_idx];                             // Local mass
                                }
                                //** If the particle is only in the parent node **/
                                else
                                {
                                    //** >> Space checking of the particle capacity in the sibling cell node **/
                                    if (space_check(&(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_cap), ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_size + 1, 2.0f, "p1i1", &(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptr_ptcl)) == _FAILURE_)
                                    {
                                        printf("Error, in space_check function\n");
                                        return _FAILURE_;
                                    }
                                    //** >> Adding the particle in the parent cell node **/
                                    ptr_node_pt->ptr_cell_struct[box_idx_pt].ptr_ptcl[ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_size] = ptcl_idx;
                                    ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_size += 1;
                                    ptr_node_pt->ptr_cell_struct[box_idx_pt].cell_mass += GL_ptcl_mass[ptcl_idx]; // Cell mass
                                }

                                //** The local mass of the node is reduced **/
                                ptr_node->local_mass -= GL_ptcl_mass[ptcl_idx];
                            }
                            //** >> The particle stay in the node **/
                            else
                            {
                                //** >> Space checking of the particle capacity in the sibling cell **/
                                if (space_check(&(ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_cap), ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_size + 1, 2.0f, "p1i1", &(ptr_node->ptr_cell_struct[box_idx_node_new].ptr_ptcl)) == _FAILURE_)
                                {
                                    printf("Error, in space_check function\n");
                                    return _FAILURE_;
                                }

                                //** >> Adding the particle in the sibling cell **/
                                ptr_node->ptr_cell_struct[box_idx_node_new].ptr_ptcl[ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_size] = ptcl_idx;
                                ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_size += 1;
                                ptr_node->ptr_cell_struct[box_idx_node_new].cell_mass += GL_ptcl_mass[ptcl_idx]; // Cell mass
                            }

                            // Whether the particle stays at the parent node or moves to a sibling node, it must be removed from the current cell of the old node
                            //** >> Removing the particle index of the current node ptr_node **/
                            // We move the last element of the array to the current position
                            ptr_node->ptr_cell_struct[box_idx_node_old].ptr_ptcl[j] = ptr_node->ptr_cell_struct[box_idx_node_old].ptr_ptcl[no_ptcl - 1];
                            no_ptcl--; // The total number of particle decrease
                            j--;       // The last element that was just moved to the current position should also must be analized
                            ptr_node->ptr_cell_struct[box_idx_node_old].cell_mass -= GL_ptcl_mass[ptcl_idx];
                        }
                        //** >> The status of the particle is changed from not updated to updated **/
                        GL_ptcl_updating_flag[ptcl_idx] = status;
                    }
                }
                ptr_node->ptr_cell_struct[box_idx_node_old].ptcl_size = no_ptcl;
            }
        }
    }   //End if(lmin < lmax)
    else // Case only 1 level, i.e. lmin = lmax
    {
        for (int i = 0; i < GL_no_ptcl; i++)
        {
            //** >> Updating the new position of the particle **/
            //** >> Velocities **/
            GL_ptcl_vx[i] += GL_ptcl_ax[i] * dt / 2;
            GL_ptcl_vy[i] += GL_ptcl_ay[i] * dt / 2;
            GL_ptcl_vz[i] += GL_ptcl_az[i] * dt / 2;

            //** >> Positions **/
            GL_ptcl_x[i] += GL_ptcl_vx[i] * dt;
            GL_ptcl_y[i] += GL_ptcl_vy[i] * dt;
            GL_ptcl_z[i] += GL_ptcl_vz[i] * dt;

            //** >> Checking if the particle exits the simulation **/
            if (GL_ptcl_x[i] < 0 || GL_ptcl_x[i] > 1 || GL_ptcl_y[i] < 0 || GL_ptcl_y[i] > 1 || GL_ptcl_z[i] < 0 || GL_ptcl_z[i] > 1)
            {
                printf("Error, Partícula %d, sale de la simulación at positions:\n", i);
                printf("x = %f\n", (double)GL_ptcl_x[i]);
                printf("y = %f\n", (double)GL_ptcl_y[i]);
                printf("z = %f\n", (double)GL_ptcl_z[i]);
                printf("vx = %f\n", (double)GL_ptcl_vx[i]);
                printf("vy = %f\n", (double)GL_ptcl_vy[i]);
                printf("vz = %f\n", (double)GL_ptcl_vz[i]);
                printf("ax = %f\n", (double)GL_ptcl_ax[i]);
                printf("ay = %f\n", (double)GL_ptcl_ay[i]);
                printf("az = %f\n", (double)GL_ptcl_az[i]);
                printf("index = %d\n", i);

                return _FAILURE_;
            }
        }
    }

    return _SUCCESS_;
}

int particle_updating_A(vtype dt)
{

    // Position and velocity are updated first by a "Predictor step"

    //** >> Particle updating A **/

    bool status;    // Boolean value for the updating particles

    // number of parents of the level = GL_tentacles_size[lv]

    status = !GL_ptcl_updating_flag[0];

    for (int lv = GL_tentacles_level_max; lv > -1; lv--)
    {
        //** >> For cycle over parent nodes **/
        for (int i = 0; i < GL_tentacles_size[lv]; i++)
        {

            // ptr_node = GL_tentacles[lv][i];

            if (computing_particles_updating_A(GL_tentacles[lv][i], dt, status) == _FAILURE_)
            {
                printf("Error at function computing_particles_updating_A()\n");
                return _FAILURE_;
            }
        }
    }

    return _SUCCESS_;
}