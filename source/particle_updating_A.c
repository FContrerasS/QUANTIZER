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
    int lv;      // Level of refinement

    struct node *ptr_node_ch; // Child node of the node ptr_node
    struct node *ptr_node_pt; // parent node of the node ptr_node
    struct node *ptr_node_sib; // sibling node of the node ptr_node

    int ptcl_idx; // Particle grid_idx in the node

    int zone_idx; // Index of the refinement zone

    vtype pos_x; // Particle position in the grid level
    vtype pos_y;
    vtype pos_z;

    int pos_x_floor; // floor of the particle position in the grid level
    int pos_y_floor;
    int pos_z_floor;

    int box_idx_x; // Box index in X direcction
    int box_idx_y; // Box index in Y direcction
    int box_idx_z; // Box index in Z direcction
    int box_idx;   // Box index

    no_ptcl = ptr_node->ptcl_size;
    lv = ptr_node->l;

    

    if (ptr_node->chn_size == 0) 
    {
        
        for (int i = 0; i < no_ptcl; i++)
        {
            
            ptcl_idx = ptr_node->ptr_ptcl[i];
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
               printf("x = %f\n", GL_ptcl_x[ptcl_idx]);
               printf("y = %f\n", GL_ptcl_y[ptcl_idx]);
               printf("z = %f\n", GL_ptcl_z[ptcl_idx]);
               printf("vx = %f\n", GL_ptcl_vx[ptcl_idx]);
               printf("vy = %f\n", GL_ptcl_vy[ptcl_idx]);
               printf("vz = %f\n", GL_ptcl_vz[ptcl_idx]);
               printf("ax = %f\n", GL_ptcl_ax[ptcl_idx]);
               printf("ay = %f\n", GL_ptcl_ay[ptcl_idx]);
               printf("az = %f\n", GL_ptcl_az[ptcl_idx]);
               printf("index = %d\n", ptcl_idx);

                return _FAILURE_;
            }

            //** >> Moving the particle to the new node if it is necessary **/

            //** >> Firstly we ask if the particle stays on the node **/

            //** >> Position of the particles in the grid level of the current node **/
            pos_x = GL_ptcl_x[ptcl_idx] * (1 << lv);
            pos_y = GL_ptcl_y[ptcl_idx] * (1 << lv);
            pos_z = GL_ptcl_z[ptcl_idx] * (1 << lv);

            //** >> Floor of the particles positions in the grid level of the current node **/
            pos_x_floor = (int)pos_x;
            pos_y_floor = (int)pos_y;
            pos_z_floor = (int)pos_z;

            //** >> Box index in the current node **/
            box_idx_x = pos_x_floor - ptr_node->box_ts_x;
            box_idx_y = pos_y_floor - ptr_node->box_ts_y;
            box_idx_z = pos_z_floor - ptr_node->box_ts_z;
            box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;


            //** Secondly we ask if the particle leaves the node
            //** >> The particle moves towards its parent node or towards some sibling node  **/
            if (ptr_node->ptr_box_old[box_idx] < -3)
            {
                ptr_node_pt = ptr_node->ptr_pt;

                //** >> Position of the particles in the grid level of the parent node **/
                pos_x = GL_ptcl_x[ptcl_idx] * (1 << (lv - 1));
                pos_y = GL_ptcl_y[ptcl_idx] * (1 << (lv - 1));
                pos_z = GL_ptcl_z[ptcl_idx] * (1 << (lv - 1));

                //** >> Floor of the particles positions in the grid level of the parent node **/
                pos_x_floor = (int)pos_x;
                pos_y_floor = (int)pos_y;
                pos_z_floor = (int)pos_z;

                //** >> Box index in the parent node **/
                box_idx_x = pos_x_floor - ptr_node_pt->box_ts_x;
                box_idx_y = pos_y_floor - ptr_node_pt->box_ts_y;
                box_idx_z = pos_z_floor - ptr_node_pt->box_ts_z;
                box_idx = box_idx_x + box_idx_y * ptr_node_pt->box_real_dim_x + box_idx_z * ptr_node_pt->box_real_dim_x * ptr_node_pt->box_real_dim_y;

                //** >> If the particle moves towards a sibling node **/
                if (ptr_node_pt->ptr_box_old[box_idx] >= 0)
                {
                    zone_idx = ptr_node_pt->ptr_box_old[box_idx];
                    ptr_node_sib = ptr_node_pt->pptr_chn[zone_idx];

                    //** >> Space checking of the particle capacity in the sibling node **/
                    if (space_check(&(ptr_node_sib->ptcl_cap), ptr_node_sib->ptcl_size + 1, "p1i1", &(ptr_node_sib->ptr_ptcl)) == _FAILURE_)
                    {
                        printf("Error, in space_check function\n");
                        return _FAILURE_;
                    }

                    //** >> Adding the particle in the sibling node **/
                    ptr_node_sib->ptr_ptcl[ptr_node_sib->ptcl_size] = ptcl_idx;
                    ptr_node_sib->ptcl_size += 1; // +1 to the total particles in the sibling node
                }

                // Whether the particle stays at the parent node or moves to a sibling node, it must be removed from the current node
                //** >> Removing the particle index of the current node ptr_node **/
                // We move the last element of the array to the current position
                ptr_node->ptr_ptcl[i] = ptr_node->ptr_ptcl[no_ptcl - 1];
                no_ptcl--; // The total number of particle decrease
                i--;       // The last element that was just moved to the current position should also must be analized
            }

            //** >> The status of the particle is changed from not updated to updated **/
            GL_ptcl_updating_flag[ptcl_idx] = status;
        }
    } // End cycle over particles in the node
    else
    {
        for (int i = 0; i < no_ptcl; i++)
        {
            ptcl_idx = ptr_node->ptr_ptcl[i];
            //** >> Particle has not been updated yet
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
                    printf("x = %f\n", GL_ptcl_x[ptcl_idx]);
                    printf("y = %f\n", GL_ptcl_y[ptcl_idx]);
                    printf("z = %f\n", GL_ptcl_z[ptcl_idx]);
                    printf("vx = %f\n", GL_ptcl_vx[ptcl_idx]);
                    printf("vy = %f\n", GL_ptcl_vy[ptcl_idx]);
                    printf("vz = %f\n", GL_ptcl_vz[ptcl_idx]);
                    printf("ax = %f\n", GL_ptcl_ax[ptcl_idx]);
                    printf("ay = %f\n", GL_ptcl_ay[ptcl_idx]);
                    printf("az = %f\n", GL_ptcl_az[ptcl_idx]);
                    printf("index = %d\n", ptcl_idx);

                    return _FAILURE_;
                }

                //** >> Moving the particle to the new node if it is necessary **/

                //** >> Firstly we ask if the particle stays on the node **/

                //** >> Position of the particles in the grid level of the current node **/
                pos_x = GL_ptcl_x[ptcl_idx] * (1 << lv);
                pos_y = GL_ptcl_y[ptcl_idx] * (1 << lv);
                pos_z = GL_ptcl_z[ptcl_idx] * (1 << lv);

                //** >> Floor of the particles positions in the grid level of the current node **/
                pos_x_floor = (int)pos_x;
                pos_y_floor = (int)pos_y;
                pos_z_floor = (int)pos_z;

                //** >> Box index in the current node **/
                box_idx_x = pos_x_floor - ptr_node->box_ts_x;
                box_idx_y = pos_y_floor - ptr_node->box_ts_y;
                box_idx_z = pos_z_floor - ptr_node->box_ts_z;
                box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

                //** >> The particle moves towards one of its child nodes **/
                if (ptr_node->ptr_box_old[box_idx] >= 0)
                {
                    zone_idx = ptr_node->ptr_box_old[box_idx];
                    ptr_node_ch = ptr_node->pptr_chn[zone_idx];

                    //** >> Space checking of the particle capacity in the child node **/
                    if (space_check(&(ptr_node_ch->ptcl_cap), ptr_node_ch->ptcl_size + 1, "p1i1", &(ptr_node_ch->ptr_ptcl)) == _FAILURE_)
                    {
                        printf("Error, in space_check function\n");
                        return _FAILURE_;
                    }

                    //** >> Adding the particle in the sibling node **/
                    ptr_node_ch->ptr_ptcl[ptr_node_ch->ptcl_size] = ptcl_idx;
                    ptr_node_ch->ptcl_size += 1; // +1 to the total particles in the child node
                }

                //** Secondly we ask if the particle leaves the node
                //** >> The particle moves towards its parent node or towards some sibling node  **/
                if (ptr_node->ptr_box_old[box_idx] < -3)
                {
                    ptr_node_pt = ptr_node->ptr_pt;

                    //** >> Position of the particles in the grid level of the parent node **/
                    pos_x = GL_ptcl_x[ptcl_idx] * (1 << (lv - 1));
                    pos_y = GL_ptcl_y[ptcl_idx] * (1 << (lv - 1));
                    pos_z = GL_ptcl_z[ptcl_idx] * (1 << (lv - 1));

                    //** >> Floor of the particles positions in the grid level of the parent node **/
                    pos_x_floor = (int)pos_x;
                    pos_y_floor = (int)pos_y;
                    pos_z_floor = (int)pos_z;

                    //** >> Box index in the parent node **/
                    box_idx_x = pos_x_floor - ptr_node_pt->box_ts_x;
                    box_idx_y = pos_y_floor - ptr_node_pt->box_ts_y;
                    box_idx_z = pos_z_floor - ptr_node_pt->box_ts_z;
                    box_idx = box_idx_x + box_idx_y * ptr_node_pt->box_real_dim_x + box_idx_z * ptr_node_pt->box_real_dim_x * ptr_node_pt->box_real_dim_y;

                    

                    //** >> If the particle moves towards a sibling node **/
                    if (ptr_node_pt->ptr_box_old[box_idx] >= 0)
                    {
                        printf("\n\nSe entra en un heramno\n\n");
                        zone_idx = ptr_node_pt->ptr_box_old[box_idx];
                        ptr_node_sib = ptr_node_pt->pptr_chn[zone_idx];

                        //** >> Space checking of the particle capacity in the sibling node **/
                        if (space_check(&(ptr_node_sib->ptcl_cap), ptr_node_sib->ptcl_size + 1, "p1i1", &(ptr_node_sib->ptr_ptcl)) == _FAILURE_)
                        {
                            printf("Error, in space_check function\n");
                            return _FAILURE_;
                        }

                        //** >> Adding the particle in the sibling node **/
                        ptr_node_sib->ptr_ptcl[ptr_node_sib->ptcl_size] = ptcl_idx;
                        ptr_node_sib->ptcl_size += 1; // +1 to the total particles in the sibling node
                    }

                    // Whether the particle stays at the parent node or moves to a sibling node, it must be removed from the current node
                    //** >> Removing the particle index of the current node ptr_node **/
                    // We move the last element of the array to the current position
                    ptr_node->ptr_ptcl[i] = ptr_node->ptr_ptcl[no_ptcl - 1];
                    no_ptcl--; // The total number of particle decrease
                    i--;       // The last element that was just moved to the current position should also must be analized
                }

                //** >> The status of the particle is changed from not updated to updated **/
                GL_ptcl_updating_flag[ptcl_idx] = status;
            }
        } // End cycle over particles in the node
    }

    //** >> Updating the number of particles in the current node **/
    ptr_node->ptcl_size = no_ptcl;

    return _SUCCESS_;
}


int particle_updating_A(vtype dt)
{

    // Position and velocity are updated first by a "Predictor step"

    //** >> Particle updating A **/

    struct node *ptr_node;
    bool status;    // Boolean value for the updating particles

    int no_pts; // Number of parents in the cycle
    status = !GL_ptcl_updating_flag[0];

    

    for (int lv = GL_tentacles_level_max; lv > -1; lv--)
    {
        no_pts = GL_tentacles_size[lv];
        //** >> For cycle over parent nodes **/
        for (int i = 0; i < no_pts; i++)
        {
            ptr_node = GL_tentacles_old[lv][i];

            if (computing_particles_updating_A(ptr_node, dt, status) == _FAILURE_)
            {
                printf("Error at function computing_particles_updating_A()\n");
                return _FAILURE_;
            }
        }
    }

    return _SUCCESS_;
}