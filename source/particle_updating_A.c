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

//** >> Local Functions
static int computing_particles_updating_A_HEAD_ONLY(struct node *ptr_node, vtype dt);
static int computing_particles_updating_A_PERIODIC(struct node *ptr_node, vtype dt, bool status);
static int computing_particles_updating_A_REFLEXIVE(struct node *ptr_node, vtype dt, bool status);
static int computing_particles_updating_A_OUTFLOW(struct node *ptr_node, vtype dt, bool status);
static int computing_particles_updating_A_NORMAL(struct node *ptr_node, vtype dt, bool status);

static int computing_particles_updating_A_HEAD_ONLY(struct node *ptr_node, vtype dt)
{

    switch (boundary_type)
    {
        case 0:
        {
            for (int i = 0; i < GL_no_ptcl_final; i++)
            {
                //** >> Updating the new position of the particle **/
                //** >> Velocities **/
                GL_ptcl_vx[i] += GL_ptcl_ax[i] * dt * 0.5;
                GL_ptcl_vy[i] += GL_ptcl_ay[i] * dt * 0.5;
                GL_ptcl_vz[i] += GL_ptcl_az[i] * dt * 0.5;

                //** >> Positions **/
                GL_ptcl_x[i] += GL_ptcl_vx[i] * dt;
                GL_ptcl_y[i] += GL_ptcl_vy[i] * dt;
                GL_ptcl_z[i] += GL_ptcl_vz[i] * dt;

                //** >> Checking if the particle is translated with the boundary of the box simulation **/
                if (GL_ptcl_x[i] < 0. || GL_ptcl_x[i] >= 1.)
                {
                    //** >> Positions **/
                    GL_ptcl_x[i] = GL_ptcl_x[i] >= 1. ? GL_ptcl_x[i] - 1 : 1 + GL_ptcl_x[i];

                    printf("\n Particle %d, is translated at x axis. New position  =  %f :\n", i, GL_ptcl_x[i]);
                }
                if (GL_ptcl_y[i] < 0. || GL_ptcl_y[i] >= 1.)
                {
                    //** >> Positions **/
                    GL_ptcl_y[i] = GL_ptcl_y[i] >= 1. ? GL_ptcl_y[i] - 1 : 1 + GL_ptcl_y[i];

                    printf("\n Particle %d, is translated at y axis. New position  =  %f :\n", i, GL_ptcl_y[i]);
                }

                if (GL_ptcl_z[i] < 0. || GL_ptcl_z[i] >= 1.)
                {
                    //** >> Positions **/
                    GL_ptcl_z[i] = GL_ptcl_z[i] >= 1. ? GL_ptcl_z[i] - 1 : 1 + GL_ptcl_z[i];

                    printf("\n Particle %d, is translated at z axis. New position  =  %f :\n", i, GL_ptcl_z[i]);
                }
            }
            break;
        }
        case 1:
        {
            for (int i = 0; i < GL_no_ptcl_final; i++)
            {
                //** >> Updating the new position of the particle **/
                //** >> Velocities **/
                GL_ptcl_vx[i] += GL_ptcl_ax[i] * dt * 0.5;
                GL_ptcl_vy[i] += GL_ptcl_ay[i] * dt * 0.5;
                GL_ptcl_vz[i] += GL_ptcl_az[i] * dt * 0.5;

                //** >> Positions **/
                GL_ptcl_x[i] += GL_ptcl_vx[i] * dt;
                GL_ptcl_y[i] += GL_ptcl_vy[i] * dt;
                GL_ptcl_z[i] += GL_ptcl_vz[i] * dt;

                //** >> Checking if the particle is reflected with the boundary of the box simulation **/
                if (GL_ptcl_x[i] < 0. || GL_ptcl_x[i] >= 1.)
                {
                    //** >> Positions **/
                    GL_ptcl_x[i] = GL_ptcl_x[i] >= 1. ? 2 - GL_ptcl_x[i] : -GL_ptcl_x[i];

                    //** >> Velocities **/
                    GL_ptcl_vx[i] *= -1;

                    //** >> Accelerations **/
                    GL_ptcl_ax[i] *= -1;

                    printf("\n Particle %d, is reflected at x axis. New position  =  %f :\n", i, GL_ptcl_x[i]);
                }
                if (GL_ptcl_y[i] < 0. || GL_ptcl_y[i] >= 1.)
                {
                    //** >> Positions **/
                    GL_ptcl_y[i] = GL_ptcl_y[i] >= 1. ? 2 - GL_ptcl_y[i] : -GL_ptcl_y[i];

                    //** >> Velocities **/
                    GL_ptcl_vy[i] *= -1;

                    //** >> Accelerations **/
                    GL_ptcl_ay[i] *= -1;

                    printf("\n Particle %d, is reflected at y axis. New position  =  %f :\n", i, GL_ptcl_y[i]);
                }

                if (GL_ptcl_z[i] < 0. || GL_ptcl_z[i] >= 1.)
                {
                    //** >> Positions **/
                    GL_ptcl_z[i] = GL_ptcl_z[i] >= 1. ? 2 - GL_ptcl_z[i] : -GL_ptcl_z[i];

                    //** >> Velocities **/
                    GL_ptcl_vz[i] *= -1;

                    //** >> Accelerations **/
                    GL_ptcl_az[i] *= -1;

                    printf("\n Particle %d, is reflected at z axis. New position  =  %f :\n", i, GL_ptcl_z[i]);
                }
            }
            break;
        }

        case 2:
        {
            for (int i = 0; i < GL_no_ptcl_final; i++)
            {
                //** >> Updating the new position of the particle **/
                //** >> Velocities **/
                GL_ptcl_vx[i] += GL_ptcl_ax[i] * dt * 0.5;
                GL_ptcl_vy[i] += GL_ptcl_ay[i] * dt * 0.5;
                GL_ptcl_vz[i] += GL_ptcl_az[i] * dt * 0.5;

                //** >> Positions **/
                GL_ptcl_x[i] += GL_ptcl_vx[i] * dt;
                GL_ptcl_y[i] += GL_ptcl_vy[i] * dt;
                GL_ptcl_z[i] += GL_ptcl_vz[i] * dt;

                //** >> Checking if the particle exits the simulation **/
                if (GL_ptcl_x[i] < 0. || GL_ptcl_x[i] >= 1. || GL_ptcl_y[i] < 0. || GL_ptcl_y[i] >= 1. || GL_ptcl_z[i] < 0. || GL_ptcl_z[i] >= 1.)
                {
                    printf("\n Particle ID = %d, local id = %d, exits the simulation:\n", GL_ptcl_ID[i], i);
                    GL_ptcl_x[i] = GL_ptcl_x[GL_no_ptcl_final - 1];
                    GL_ptcl_y[i] = GL_ptcl_y[GL_no_ptcl_final - 1];
                    GL_ptcl_z[i] = GL_ptcl_z[GL_no_ptcl_final - 1];
                    GL_ptcl_vx[i] = GL_ptcl_vx[GL_no_ptcl_final - 1];
                    GL_ptcl_vy[i] = GL_ptcl_vy[GL_no_ptcl_final - 1];
                    GL_ptcl_vz[i] = GL_ptcl_vz[GL_no_ptcl_final - 1];
                    GL_ptcl_ax[i] = GL_ptcl_ax[GL_no_ptcl_final - 1];
                    GL_ptcl_ay[i] = GL_ptcl_ay[GL_no_ptcl_final - 1];
                    GL_ptcl_az[i] = GL_ptcl_az[GL_no_ptcl_final - 1];
                    GL_ptcl_ID[i] = GL_ptcl_ID[GL_no_ptcl_final - 1];
                    GL_no_ptcl_final--;

                    if (GL_no_ptcl_final == 0)
                    {
                        printf("Error, There are no particles in the simulation\n");
                        return _FAILURE_;
                    }
                }
            }

            break;
        }
        
        default:
        {
            printf("Error! the boundary_type paramter = %d is different to 0, 1 or 2\n", boundary_type);
            return _FAILURE_;
        }
    }

    return _SUCCESS_;
}

static int computing_particles_updating_A_PERIODIC(struct node *ptr_node, vtype dt, bool status)
{

    int no_ptcl; // Total number of particles in the node

    struct node *ptr_node_ch;  // Child node of the node ptr_node
    struct node *ptr_node_pt;  // parent node of the node ptr_node
    struct node *ptr_node_sib; // sibling node of the node ptr_node

    int ptcl_idx; // Particle grid_idx in the node

    int zone_idx; // Index of the refinement zone

    int box_idx_node_old; // Box index of the current node
    int box_idx_node_new; // Box index of the current node
    int box_idx_pt;       // Box index of the parent node
    int box_idx_sib;      // Box index of the sibling node
    int box_idx_ch;       // Box index of the child node

    if (ptr_node->chn_size == 0)
    {
        for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
        {
            box_idx_node_old = ptr_node->ptr_box_idx[cell_idx];

            no_ptcl = ptr_node->ptr_cell_struct[box_idx_node_old].ptcl_size;
            for (int j = 0; j < no_ptcl; j++)
            {
                ptcl_idx = ptr_node->ptr_cell_struct[box_idx_node_old].ptr_ptcl[j];
                if (GL_ptcl_updating_flag[ptcl_idx] != status)
                {

                    //** >> Updating the new position of the particle **/
                    //** >> Velocities **/
                    GL_ptcl_vx[ptcl_idx] += GL_ptcl_ax[ptcl_idx] * dt * 0.5;
                    GL_ptcl_vy[ptcl_idx] += GL_ptcl_ay[ptcl_idx] * dt * 0.5;
                    GL_ptcl_vz[ptcl_idx] += GL_ptcl_az[ptcl_idx] * dt * 0.5;

                    //** >> Positions **/
                    GL_ptcl_x[ptcl_idx] += GL_ptcl_vx[ptcl_idx] * dt;
                    GL_ptcl_y[ptcl_idx] += GL_ptcl_vy[ptcl_idx] * dt;
                    GL_ptcl_z[ptcl_idx] += GL_ptcl_vz[ptcl_idx] * dt;

                    //** >> Checking if the particle exits the simulation **/
                    if (GL_ptcl_x[ptcl_idx] < 0. || GL_ptcl_x[ptcl_idx] >= 1. || GL_ptcl_y[ptcl_idx] < 0. || GL_ptcl_y[ptcl_idx] >= 1. || GL_ptcl_z[ptcl_idx] < 0. || GL_ptcl_z[ptcl_idx] >= 1.)
                    {
                        printf("Error, Partícula %d, sale de la simulación at positions:\n", ptcl_idx);
                        printf("lv = %d, ID = %d\n", ptr_node->lv, ptr_node->ID);
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
                                if (space_check(&(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_cap), ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptr_ptcl)) == _FAILURE_)
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
                                if (space_check(&(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_cap), ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptr_ptcl)) == _FAILURE_)
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
                            if (space_check(&(ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_cap), ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node->ptr_cell_struct[box_idx_node_new].ptr_ptcl)) == _FAILURE_)
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
            box_idx_node_old = ptr_node->ptr_box_idx[cell_idx];

            if (ptr_node->ptr_box[box_idx_node_old] < 0)
            {
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
                        if (GL_ptcl_x[ptcl_idx] < 0. || GL_ptcl_x[ptcl_idx] >= 1. || GL_ptcl_y[ptcl_idx] < 0. || GL_ptcl_y[ptcl_idx] >= 1. || GL_ptcl_z[ptcl_idx] < 0. || GL_ptcl_z[ptcl_idx] >= 1.)
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
                                if (space_check(&(ptr_node_ch->ptr_cell_struct[box_idx_ch].ptcl_cap), ptr_node_ch->ptr_cell_struct[box_idx_ch].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_ch->ptr_cell_struct[box_idx_ch].ptr_ptcl)) == _FAILURE_)
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
                                    if (space_check(&(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_cap), ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptr_ptcl)) == _FAILURE_)
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
                                    if (space_check(&(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_cap), ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptr_ptcl)) == _FAILURE_)
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
                                if (space_check(&(ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_cap), ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node->ptr_cell_struct[box_idx_node_new].ptr_ptcl)) == _FAILURE_)
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
    }

    return _SUCCESS_;
}

static int computing_particles_updating_A_REFLEXIVE(struct node *ptr_node, vtype dt, bool status)

{

    int no_ptcl; // Total number of particles in the node

    struct node *ptr_node_ch; // Child node of the node ptr_node
    struct node *ptr_node_pt; // parent node of the node ptr_node
    struct node *ptr_node_sib; // sibling node of the node ptr_node

    int ptcl_idx; // Particle grid_idx in the node

    int zone_idx; // Index of the refinement zone

    int box_idx_node_old;   // Box index of the current node
    int box_idx_node_new;   // Box index of the current node
    int box_idx_pt;   // Box index of the parent node
    int box_idx_sib;   // Box index of the sibling node
    int box_idx_ch;   // Box index of the child node

    if (ptr_node->chn_size == 0)
    {
        for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
        {
            box_idx_node_old = ptr_node->ptr_box_idx[cell_idx];

            no_ptcl = ptr_node->ptr_cell_struct[box_idx_node_old].ptcl_size;
            for (int j = 0; j < no_ptcl; j++)
            {
                ptcl_idx = ptr_node->ptr_cell_struct[box_idx_node_old].ptr_ptcl[j];
                if (GL_ptcl_updating_flag[ptcl_idx] != status)
                {

                    //** >> Updating the new position of the particle **/
                    //** >> Velocities **/
                    GL_ptcl_vx[ptcl_idx] += GL_ptcl_ax[ptcl_idx] * dt * 0.5;
                    GL_ptcl_vy[ptcl_idx] += GL_ptcl_ay[ptcl_idx] * dt * 0.5;
                    GL_ptcl_vz[ptcl_idx] += GL_ptcl_az[ptcl_idx] * dt * 0.5;

                    //** >> Positions **/
                    GL_ptcl_x[ptcl_idx] += GL_ptcl_vx[ptcl_idx] * dt;
                    GL_ptcl_y[ptcl_idx] += GL_ptcl_vy[ptcl_idx] * dt;
                    GL_ptcl_z[ptcl_idx] += GL_ptcl_vz[ptcl_idx] * dt;

                    //** >> Checking if the particle is reflected with the boundary of the box simulation **/
                    if (GL_ptcl_x[ptcl_idx] < 0. || GL_ptcl_x[ptcl_idx] >= 1.)
                    {
                        //** >> Positions **/
                        GL_ptcl_x[ptcl_idx] = GL_ptcl_x[ptcl_idx] >= 1. ? 2 - GL_ptcl_x[ptcl_idx] : -GL_ptcl_x[ptcl_idx];

                        //** >> Velocities **/
                        GL_ptcl_vx[ptcl_idx] *= -1;

                        //** >> Accelerations **/
                        GL_ptcl_ax[ptcl_idx] *= -1;

                        printf("\n Particle %d, is reflected at x axis. New position  =  %f :\n", ptcl_idx, GL_ptcl_x[ptcl_idx]);
                    }
                    if (GL_ptcl_y[ptcl_idx] < 0. || GL_ptcl_y[ptcl_idx] >= 1.)
                    {
                        //** >> Positions **/
                        GL_ptcl_y[ptcl_idx] = GL_ptcl_y[ptcl_idx] >= 1. ? 2 - GL_ptcl_y[ptcl_idx] : -GL_ptcl_y[ptcl_idx];

                        //** >> Velocities **/
                        GL_ptcl_vy[ptcl_idx] *= -1;

                        //** >> Accelerations **/
                        GL_ptcl_ay[ptcl_idx] *= -1;

                        printf("\n Particle %d, is reflected at y axis. New position  =  %f :\n", ptcl_idx, GL_ptcl_y[ptcl_idx]);
                    }

                    if (GL_ptcl_z[ptcl_idx] < 0. || GL_ptcl_z[ptcl_idx] >= 1.)
                    {
                        //** >> Positions **/
                        GL_ptcl_z[ptcl_idx] = GL_ptcl_z[ptcl_idx] >= 1. ? 2 - GL_ptcl_z[ptcl_idx] : -GL_ptcl_z[ptcl_idx];

                        //** >> Velocities **/
                        GL_ptcl_vz[ptcl_idx] *= -1;

                        //** >> Accelerations **/
                        GL_ptcl_az[ptcl_idx] *= -1;

                        printf("\n Particle %d, is reflected at z axis. New position  =  %f :\n", ptcl_idx, GL_ptcl_z[ptcl_idx]);
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
                                if (space_check(&(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_cap), ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptr_ptcl)) == _FAILURE_)
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
                                if (space_check(&(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_cap), ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptr_ptcl)) == _FAILURE_)
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
                            if (space_check(&(ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_cap), ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node->ptr_cell_struct[box_idx_node_new].ptr_ptcl)) == _FAILURE_)
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
            box_idx_node_old = ptr_node->ptr_box_idx[cell_idx];

            if (ptr_node->ptr_box[box_idx_node_old] < 0)
            {
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

                        //** >> Checking if the particle is reflected with the boundary of the box simulation **/
                        if (GL_ptcl_x[ptcl_idx] < 0. || GL_ptcl_x[ptcl_idx] >= 1.)
                        {
                            //** >> Positions **/
                            GL_ptcl_x[ptcl_idx] = GL_ptcl_x[ptcl_idx] >= 1. ? 2 - GL_ptcl_x[ptcl_idx] : -GL_ptcl_x[ptcl_idx];

                            //** >> Velocities **/
                            GL_ptcl_vx[ptcl_idx] *= -1;

                            //** >> Accelerations **/
                            GL_ptcl_ax[ptcl_idx] *= -1;

                            printf("\n Particle %d, is reflected at x axis. New position  =  %f :\n", ptcl_idx, GL_ptcl_x[ptcl_idx]);
                        }
                        if (GL_ptcl_y[ptcl_idx] < 0. || GL_ptcl_y[ptcl_idx] >= 1.)
                        {
                            //** >> Positions **/
                            GL_ptcl_y[ptcl_idx] = GL_ptcl_y[ptcl_idx] >= 1. ? 2 - GL_ptcl_y[ptcl_idx] : -GL_ptcl_y[ptcl_idx];

                            //** >> Velocities **/
                            GL_ptcl_vy[ptcl_idx] *= -1;

                            //** >> Accelerations **/
                            GL_ptcl_ay[ptcl_idx] *= -1;

                            printf("\n Particle %d, is reflected at y axis. New position  =  %f :\n", ptcl_idx, GL_ptcl_y[ptcl_idx]);
                        }

                        if (GL_ptcl_z[ptcl_idx] < 0. || GL_ptcl_z[ptcl_idx] >= 1.)
                        {
                            //** >> Positions **/
                            GL_ptcl_z[ptcl_idx] = GL_ptcl_z[ptcl_idx] >= 1. ? 2 - GL_ptcl_z[ptcl_idx] : -GL_ptcl_z[ptcl_idx];

                            //** >> Velocities **/
                            GL_ptcl_vz[ptcl_idx] *= -1;

                            //** >> Accelerations **/
                            GL_ptcl_az[ptcl_idx] *= -1;

                            printf("\n Particle %d, is reflected at z axis. New position  =  %f :\n", ptcl_idx, GL_ptcl_z[ptcl_idx]);
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
                                if (space_check(&(ptr_node_ch->ptr_cell_struct[box_idx_ch].ptcl_cap), ptr_node_ch->ptr_cell_struct[box_idx_ch].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_ch->ptr_cell_struct[box_idx_ch].ptr_ptcl)) == _FAILURE_)
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
                                    if (space_check(&(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_cap), ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptr_ptcl)) == _FAILURE_)
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
                                    if (space_check(&(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_cap), ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptr_ptcl)) == _FAILURE_)
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
                                if (space_check(&(ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_cap), ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node->ptr_cell_struct[box_idx_node_new].ptr_ptcl)) == _FAILURE_)
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
    }

    return _SUCCESS_;
}

static int computing_particles_updating_A_OUTFLOW(struct node *ptr_node, vtype dt, bool status)
{

    int no_ptcl; // Total number of particles in the node

    struct node *ptr_node_ch;  // Child node of the node ptr_node
    struct node *ptr_node_pt;  // parent node of the node ptr_node
    struct node *ptr_node_sib; // sibling node of the node ptr_node

    struct node *ptr_node_aux;
    int box_idx_aux;
    int counter_aux;

    int ptcl_idx; // Particle grid_idx in the node

    int zone_idx; // Index of the refinement zone

    int box_idx_node_old; // Box index of the current node
    int box_idx_node_new; // Box index of the current node
    int box_idx_pt;       // Box index of the parent node
    int box_idx_sib;      // Box index of the sibling node
    int box_idx_ch;       // Box index of the child node

    if (ptr_node->chn_size == 0)
    {
        for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
        {
            box_idx_node_old = ptr_node->ptr_box_idx[cell_idx];

            no_ptcl = ptr_node->ptr_cell_struct[box_idx_node_old].ptcl_size;
            for (int j = 0; j < no_ptcl; j++)
            {
                ptcl_idx = ptr_node->ptr_cell_struct[box_idx_node_old].ptr_ptcl[j];
                if (GL_ptcl_updating_flag[ptcl_idx] != status)
                {
                    
                    //** >> Updating the new position of the particle **/
                    //** >> Velocities **/
                    GL_ptcl_vx[ptcl_idx] += GL_ptcl_ax[ptcl_idx] * dt * 0.5;
                    GL_ptcl_vy[ptcl_idx] += GL_ptcl_ay[ptcl_idx] * dt * 0.5;
                    GL_ptcl_vz[ptcl_idx] += GL_ptcl_az[ptcl_idx] * dt * 0.5;

                    //** >> Positions **/
                    GL_ptcl_x[ptcl_idx] += GL_ptcl_vx[ptcl_idx] * dt;
                    GL_ptcl_y[ptcl_idx] += GL_ptcl_vy[ptcl_idx] * dt;
                    GL_ptcl_z[ptcl_idx] += GL_ptcl_vz[ptcl_idx] * dt;

                    //** >> Checking if the particle exits the simulation **/
                    if (GL_ptcl_x[ptcl_idx] < 0. || GL_ptcl_x[ptcl_idx] >= 1. || GL_ptcl_y[ptcl_idx] < 0. || GL_ptcl_y[ptcl_idx] >= 1. || GL_ptcl_z[ptcl_idx] < 0. || GL_ptcl_z[ptcl_idx] >= 1.)
                    {
                        printf("\n Particle ID = %d, local id = %d, exits the simulation:\n", GL_ptcl_ID[ptcl_idx], ptcl_idx);

                        //** >> Updating Global total mass **/
                        total_mass -= GL_ptcl_mass[ptcl_idx];

                        //** >> Removing the particle from the local cell node **/
                        ptr_node->local_mass -= GL_ptcl_mass[ptcl_idx];
                        ptr_node->ptr_cell_struct[box_idx_node_old].ptr_ptcl[j] = ptr_node->ptr_cell_struct[box_idx_node_old].ptr_ptcl[no_ptcl - 1];
                        no_ptcl--; // The total number of particle decrease
                        j--;       // The last element that was just moved to the current position should also must be analized
                        ptr_node->ptr_cell_struct[box_idx_node_old].cell_mass -= GL_ptcl_mass[ptcl_idx];

                        //** >> Removing the local from the parent nodes **/
                        ptr_node_aux = ptr_node;
                        while (ptr_node_aux != GL_ptr_tree)
                        {
                            ptr_node_aux = ptr_node_aux->ptr_pt;
                            ptr_node_aux->local_mass -= GL_ptcl_mass[ptcl_idx];
                        }

                        if (ptcl_idx != GL_no_ptcl_final - 1)
                        {
                            //** >> Searching the node and cell which contains the particle index "GL_no_ptcl_final - 1"
                            ptr_node_aux = GL_ptr_tree;
                            box_idx_aux = ptcl_idx_to_box_idx(ptr_node_aux, GL_no_ptcl_final - 1);
                            while (ptr_node_aux->ptr_box[box_idx_aux] >= 0)
                            {
                                ptr_node_aux = ptr_node_aux->pptr_chn[ptr_node_aux->ptr_box[box_idx_aux]];
                                box_idx_aux = ptcl_idx_to_box_idx(ptr_node_aux, GL_no_ptcl_final - 1);
                            }

                            //** >> Searching the position in the cell which is equal to the particle index "GL_no_ptcl_final - 1"
                            counter_aux = 0;
                            while (ptr_node_aux->ptr_cell_struct[box_idx_aux].ptr_ptcl[counter_aux] != GL_no_ptcl_final - 1)
                            {
                                counter_aux++;
                            }

                            //** >> Updating the value of the new particle index from "GL_no_ptcl_final - 1" to "ptcl_idx"
                            ptr_node_aux->ptr_cell_struct[box_idx_aux].ptr_ptcl[counter_aux] = ptcl_idx;

                            //** >> Removing the particle from the global array **/
                            GL_ptcl_x[ptcl_idx] = GL_ptcl_x[GL_no_ptcl_final - 1];
                            GL_ptcl_y[ptcl_idx] = GL_ptcl_y[GL_no_ptcl_final - 1];
                            GL_ptcl_z[ptcl_idx] = GL_ptcl_z[GL_no_ptcl_final - 1];
                            GL_ptcl_vx[ptcl_idx] = GL_ptcl_vx[GL_no_ptcl_final - 1];
                            GL_ptcl_vy[ptcl_idx] = GL_ptcl_vy[GL_no_ptcl_final - 1];
                            GL_ptcl_vz[ptcl_idx] = GL_ptcl_vz[GL_no_ptcl_final - 1];
                            GL_ptcl_ax[ptcl_idx] = GL_ptcl_ax[GL_no_ptcl_final - 1];
                            GL_ptcl_ay[ptcl_idx] = GL_ptcl_ay[GL_no_ptcl_final - 1];
                            GL_ptcl_az[ptcl_idx] = GL_ptcl_az[GL_no_ptcl_final - 1];
                            GL_ptcl_updating_flag[ptcl_idx] = GL_ptcl_updating_flag[GL_no_ptcl_final - 1];
                            GL_ptcl_ID[ptcl_idx] = GL_ptcl_ID[GL_no_ptcl_final - 1];
                        }
                        GL_no_ptcl_final--;
                        if (GL_no_ptcl_final == 0)
                        {
                            printf("Error, There are no particles in the simulation\n");
                            return _FAILURE_;
                        }
                        
                    }
                    else
                    {
                        
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
                                    if (space_check(&(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_cap), ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptr_ptcl)) == _FAILURE_)
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
                                    if (space_check(&(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_cap), ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptr_ptcl)) == _FAILURE_)
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
                                if (space_check(&(ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_cap), ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node->ptr_cell_struct[box_idx_node_new].ptr_ptcl)) == _FAILURE_)
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
            }
            ptr_node->ptr_cell_struct[box_idx_node_old].ptcl_size = no_ptcl;
        }
    } // End cycle over particles in the node
    else
    {
        for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
        {
            box_idx_node_old = ptr_node->ptr_box_idx[cell_idx];

            if (ptr_node->ptr_box[box_idx_node_old] < 0)
            {
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
                        if (GL_ptcl_x[ptcl_idx] < 0. || GL_ptcl_x[ptcl_idx] >= 1. || GL_ptcl_y[ptcl_idx] < 0. || GL_ptcl_y[ptcl_idx] >= 1. || GL_ptcl_z[ptcl_idx] < 0. || GL_ptcl_z[ptcl_idx] >= 1.)
                        {
                            printf("\n Particle ID = %d, local id = %d, exits the simulation:\n", GL_ptcl_ID[ptcl_idx], ptcl_idx);

                            //** >> Updating Global total mass **/
                            total_mass -= GL_ptcl_mass[ptcl_idx];

                            //** >> Removing the particle from the local cell node **/
                            ptr_node->local_mass -= GL_ptcl_mass[ptcl_idx];
                            ptr_node->ptr_cell_struct[box_idx_node_old].ptr_ptcl[j] = ptr_node->ptr_cell_struct[box_idx_node_old].ptr_ptcl[no_ptcl - 1];
                            no_ptcl--; // The total number of particle decrease
                            j--;       // The last element that was just moved to the current position should also must be analized
                            ptr_node->ptr_cell_struct[box_idx_node_old].cell_mass -= GL_ptcl_mass[ptcl_idx];

                            //** >> Removing the local from the parent nodes **/
                            ptr_node_aux = ptr_node;
                            while (ptr_node_aux != GL_ptr_tree)
                            {
                                ptr_node_aux = ptr_node_aux->ptr_pt;
                                ptr_node_aux->local_mass -= GL_ptcl_mass[ptcl_idx];
                            }

                            if (ptcl_idx != GL_no_ptcl_final - 1)
                            {
                                //** >> Searching the node and cell which contains the particle index "GL_no_ptcl_final - 1"
                                ptr_node_aux = GL_ptr_tree;
                                box_idx_aux = ptcl_idx_to_box_idx(ptr_node_aux, GL_no_ptcl_final - 1);
                                while (ptr_node_aux->ptr_box[box_idx_aux] >= 0)
                                {
                                    ptr_node_aux = ptr_node_aux->pptr_chn[ptr_node_aux->ptr_box[box_idx_aux]];
                                    box_idx_aux = ptcl_idx_to_box_idx(ptr_node_aux, GL_no_ptcl_final - 1);
                                }

                                //** >> Searching the position in the cell which is equal to the particle index "GL_no_ptcl_final - 1"
                                counter_aux = 0;
                                while (ptr_node_aux->ptr_cell_struct[box_idx_aux].ptr_ptcl[counter_aux] != GL_no_ptcl_final - 1)
                                {
                                    counter_aux++;
                                }

                                //** >> Updating the value of the new particle index from "GL_no_ptcl_final - 1" to "ptcl_idx"
                                ptr_node_aux->ptr_cell_struct[box_idx_aux].ptr_ptcl[counter_aux] = ptcl_idx;

                                //** >> Removing the particle from the global array **/
                                GL_ptcl_x[ptcl_idx] = GL_ptcl_x[GL_no_ptcl_final - 1];
                                GL_ptcl_y[ptcl_idx] = GL_ptcl_y[GL_no_ptcl_final - 1];
                                GL_ptcl_z[ptcl_idx] = GL_ptcl_z[GL_no_ptcl_final - 1];
                                GL_ptcl_vx[ptcl_idx] = GL_ptcl_vx[GL_no_ptcl_final - 1];
                                GL_ptcl_vy[ptcl_idx] = GL_ptcl_vy[GL_no_ptcl_final - 1];
                                GL_ptcl_vz[ptcl_idx] = GL_ptcl_vz[GL_no_ptcl_final - 1];
                                GL_ptcl_ax[ptcl_idx] = GL_ptcl_ax[GL_no_ptcl_final - 1];
                                GL_ptcl_ay[ptcl_idx] = GL_ptcl_ay[GL_no_ptcl_final - 1];
                                GL_ptcl_az[ptcl_idx] = GL_ptcl_az[GL_no_ptcl_final - 1];
                                GL_ptcl_updating_flag[ptcl_idx] = GL_ptcl_updating_flag[GL_no_ptcl_final - 1];
                                GL_ptcl_ID[ptcl_idx] = GL_ptcl_ID[GL_no_ptcl_final - 1];
                            }
                            GL_no_ptcl_final--;
                            if (GL_no_ptcl_final == 0)
                            {
                                printf("Error, There are no particles in the simulation\n");
                                return _FAILURE_;
                            }
                        }
                        else
                        {
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
                                    if (space_check(&(ptr_node_ch->ptr_cell_struct[box_idx_ch].ptcl_cap), ptr_node_ch->ptr_cell_struct[box_idx_ch].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_ch->ptr_cell_struct[box_idx_ch].ptr_ptcl)) == _FAILURE_)
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
                                        if (space_check(&(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_cap), ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptr_ptcl)) == _FAILURE_)
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
                                        if (space_check(&(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_cap), ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptr_ptcl)) == _FAILURE_)
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
                                    if (space_check(&(ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_cap), ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node->ptr_cell_struct[box_idx_node_new].ptr_ptcl)) == _FAILURE_)
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
                }
                ptr_node->ptr_cell_struct[box_idx_node_old].ptcl_size = no_ptcl;
            }
        }
    }

    return _SUCCESS_;
}

static int computing_particles_updating_A_NORMAL(struct node *ptr_node, vtype dt, bool status)
{

    int no_ptcl; // Total number of particles in the node

    struct node *ptr_node_ch;  // Child node of the node ptr_node
    struct node *ptr_node_pt;  // parent node of the node ptr_node
    struct node *ptr_node_sib; // sibling node of the node ptr_node

    int ptcl_idx; // Particle grid_idx in the node

    int zone_idx; // Index of the refinement zone

    int box_idx_node_old; // Box index of the current node
    int box_idx_node_new; // Box index of the current node
    int box_idx_pt;       // Box index of the parent node
    int box_idx_sib;      // Box index of the sibling node
    int box_idx_ch;       // Box index of the child node

    if (ptr_node->chn_size == 0)
    {
        for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
        {
            box_idx_node_old = ptr_node->ptr_box_idx[cell_idx];

            no_ptcl = ptr_node->ptr_cell_struct[box_idx_node_old].ptcl_size;
            for (int j = 0; j < no_ptcl; j++)
            {
                ptcl_idx = ptr_node->ptr_cell_struct[box_idx_node_old].ptr_ptcl[j];
                if (GL_ptcl_updating_flag[ptcl_idx] != status)
                {

                    //** >> Updating the new position of the particle **/
                    //** >> Velocities **/
                    GL_ptcl_vx[ptcl_idx] += GL_ptcl_ax[ptcl_idx] * dt * 0.5;
                    GL_ptcl_vy[ptcl_idx] += GL_ptcl_ay[ptcl_idx] * dt * 0.5;
                    GL_ptcl_vz[ptcl_idx] += GL_ptcl_az[ptcl_idx] * dt * 0.5;

                    //** >> Positions **/
                    GL_ptcl_x[ptcl_idx] += GL_ptcl_vx[ptcl_idx] * dt;
                    GL_ptcl_y[ptcl_idx] += GL_ptcl_vy[ptcl_idx] * dt;
                    GL_ptcl_z[ptcl_idx] += GL_ptcl_vz[ptcl_idx] * dt;

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
                                if (space_check(&(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_cap), ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptr_ptcl)) == _FAILURE_)
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
                                if (space_check(&(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_cap), ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptr_ptcl)) == _FAILURE_)
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
                            if (space_check(&(ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_cap), ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node->ptr_cell_struct[box_idx_node_new].ptr_ptcl)) == _FAILURE_)
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
            box_idx_node_old = ptr_node->ptr_box_idx[cell_idx];

            if (ptr_node->ptr_box[box_idx_node_old] < 0)
            {
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
                        if (GL_ptcl_x[ptcl_idx] < 0. || GL_ptcl_x[ptcl_idx] >= 1. || GL_ptcl_y[ptcl_idx] < 0. || GL_ptcl_y[ptcl_idx] >= 1. || GL_ptcl_z[ptcl_idx] < 0. || GL_ptcl_z[ptcl_idx] >= 1.)
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
                                if (space_check(&(ptr_node_ch->ptr_cell_struct[box_idx_ch].ptcl_cap), ptr_node_ch->ptr_cell_struct[box_idx_ch].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_ch->ptr_cell_struct[box_idx_ch].ptr_ptcl)) == _FAILURE_)
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
                                    if (space_check(&(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_cap), ptr_node_sib->ptr_cell_struct[box_idx_sib].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_sib->ptr_cell_struct[box_idx_sib].ptr_ptcl)) == _FAILURE_)
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
                                    if (space_check(&(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_cap), ptr_node_pt->ptr_cell_struct[box_idx_pt].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node_pt->ptr_cell_struct[box_idx_pt].ptr_ptcl)) == _FAILURE_)
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
                                if (space_check(&(ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_cap), ptr_node->ptr_cell_struct[box_idx_node_new].ptcl_size + 1, 1.0f, "p1i1", &(ptr_node->ptr_cell_struct[box_idx_node_new].ptr_ptcl)) == _FAILURE_)
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
    }

    return _SUCCESS_;
}

int particle_updating_A(vtype dt)
{

    // Position and velocity are updated first by a "Predictor step"

    //** >> Particle updating A **/

    if(lmin < lmax)
    {
        struct node *ptr_node;
        bool status; // Boolean value for the updating particles

        // number of parents of the level = GL_tentacles_size[lv]

        status = !GL_ptcl_updating_flag[0];

        switch (boundary_type)
        {
            case 0:
            {
                for (int lv = GL_tentacles_level_max; lv > -1; lv--)
                {
                    //** >> For cycle over parent nodes **/
                    for (int i = 0; i < GL_tentacles_size[lv]; i++)
                    {

                        ptr_node = GL_tentacles[lv][i];

                        if (ptr_node->boundary_simulation_contact == false)
                        {
                            if (computing_particles_updating_A_NORMAL(GL_tentacles[lv][i], dt, status) == _FAILURE_)
                            {
                                printf("Error at function computing_particles_updating_A()\n");
                                return _FAILURE_;
                            }
                        }
                        else
                        {
                            if (computing_particles_updating_A_PERIODIC(GL_tentacles[lv][i], dt, status) == _FAILURE_)
                            {
                                printf("Error at function computing_particles_updating_A()\n");
                                return _FAILURE_;
                            }
                        }
                    }
                }
                break;
            }

            case 1:
            {
                for (int lv = GL_tentacles_level_max; lv > -1; lv--)
                {
                    //** >> For cycle over parent nodes **/
                    for (int i = 0; i < GL_tentacles_size[lv]; i++)
                    {

                        ptr_node = GL_tentacles[lv][i];

                        if (ptr_node->boundary_simulation_contact == false)
                        {
                            if (computing_particles_updating_A_NORMAL(GL_tentacles[lv][i], dt, status) == _FAILURE_)
                            {
                                printf("Error at function computing_particles_updating_A()\n");
                                return _FAILURE_;
                            }
                        }
                        else
                        {
                            if (computing_particles_updating_A_REFLEXIVE(GL_tentacles[lv][i], dt, status) == _FAILURE_)
                            {
                                printf("Error at function computing_particles_updating_A()\n");
                                return _FAILURE_;
                            }
                        }
                    }
                }
            }
            case 2:
            {
                for (int lv = GL_tentacles_level_max; lv > -1; lv--)
                {
                    //** >> For cycle over parent nodes **/
                    for (int i = 0; i < GL_tentacles_size[lv]; i++)
                    {

                        ptr_node = GL_tentacles[lv][i];

                        if (ptr_node->boundary_simulation_contact == false)
                        {
                            if (computing_particles_updating_A_NORMAL(GL_tentacles[lv][i], dt, status) == _FAILURE_)
                            {
                                printf("Error at function computing_particles_updating_A()\n");
                                return _FAILURE_;
                            }
                        }
                        else
                        {
                            if (computing_particles_updating_A_OUTFLOW(GL_tentacles[lv][i], dt, status) == _FAILURE_)
                            {
                                printf("Error at function computing_particles_updating_A()\n");
                                return _FAILURE_;
                            }
                        }
                    }
                }
                break;
            }

            default:
            {
                printf("Error! the boundary_type paramter = %d is different to 0, 1 or 2\n", boundary_type);
                return _FAILURE_;
            }
        }
    }
    else
    {
        if (computing_particles_updating_A_HEAD_ONLY(GL_tentacles[0][0], dt) == _FAILURE_)
        {
            printf("Error at function computing_particles_updating_A_HEAD_ONLY()\n");
            return _FAILURE_;
        }
    }

    return _SUCCESS_;
}