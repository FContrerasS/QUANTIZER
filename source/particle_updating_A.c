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

static int ptcl_idx_to_box_idx(struct node *ptr_node, int ptcl_idx)
{

    int lv;      // Level of refinement

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

    lv = ptr_node->lv;

    //** >> Position of the particles in the grid level of the current node before updating**/
    pos_x = GL_ptcl_x[ptcl_idx] * (1 << lv);
    pos_y = GL_ptcl_y[ptcl_idx] * (1 << lv);
    pos_z = GL_ptcl_z[ptcl_idx] * (1 << lv);

    //** >> Floor of the particles positions in the grid level of the current node before updating**/
    pos_x_floor = (int)pos_x;
    pos_y_floor = (int)pos_y;
    pos_z_floor = (int)pos_z;

    //** >> Box index in the current node before updating **/
    box_idx_x = pos_x_floor - ptr_node->box_ts_x;
    box_idx_y = pos_y_floor - ptr_node->box_ts_y;
    box_idx_z = pos_z_floor - ptr_node->box_ts_z;
    box_idx = box_idx_x + box_idx_y * ptr_node->box_real_dim_x + box_idx_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;

    return box_idx;
}

static int computing_particles_updating_A(struct node *ptr_node, vtype dt, bool status)
{

    int no_ptcl; // Total number of particles in the node

    struct node *ptr_node_ch = NULL; // Child node of the node ptr_node
    struct node *ptr_node_pt = NULL; // parent node of the node ptr_node
    struct node *ptr_node_sib = NULL; // sibling node of the node ptr_node

    int ptcl_idx; // Particle grid_idx in the node

    int zone_idx; // Index of the refinement zone

    int box_idx_node;   // Box index of the current node
    int box_idx_pt;   // Box index of the parent node
    int box_idx_sib;   // Box index of the sibling node
    int box_idx_ch;   // Box index of the child node

    no_ptcl = ptr_node->ptcl_size;

    vtype aux_mass1 = 0;
    vtype aux_mass2 = 0;

    int box_idx_aux_x;
    int box_idx_aux_y;
    int box_idx_aux_z;
    int box_idx_aux;

    for (int i = 0; i < ptr_node->box_real_dim_x * ptr_node->box_real_dim_y * ptr_node->box_real_dim_z; i++)
    {
        aux_mass1 += ptr_node->ptr_box_mass[i];
    }
    aux_mass2 = 0;
    for (int i = 0; i < ptr_node->cell_size; i++)
    {
        box_idx_aux_x = ptr_node->ptr_cell_idx_x[i] - ptr_node->box_ts_x;
        box_idx_aux_y = ptr_node->ptr_cell_idx_y[i] - ptr_node->box_ts_y;
        box_idx_aux_z = ptr_node->ptr_cell_idx_z[i] - ptr_node->box_ts_z;
        box_idx_aux = box_idx_aux_x + box_idx_aux_y * ptr_node->box_real_dim_x + box_idx_aux_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
        aux_mass2 += ptr_node->ptr_box_mass[box_idx_aux];
    }

    printf("\n\n Error, particle updating Antes\n\n");
    printf("lv = %d, aux mass per box = %f\n", ptr_node->lv, aux_mass1);
    printf("lv = %d, aux mass per cell = %f\n",ptr_node->lv, aux_mass2);
    printf("lv = %d, local mass = %f\n",ptr_node->lv, ptr_node->local_mass);
    printf("lv = %d, no ptcl = %d\n",ptr_node->lv,ptr_node->ptcl_size);
    printf("\n\n");
    



    if(ptr_node->lv == 5)
    {
        struct node *ptr_ch;
        ptr_ch = ptr_node->pptr_chn[0];
        aux_mass1 = 0;
        for (int i = 0; i < ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y * ptr_ch->box_real_dim_z; i++)
        {
            aux_mass1 += ptr_ch->ptr_box_mass[i];
        }
        aux_mass2 = 0;
        for (int i = 0; i < ptr_ch->cell_size; i++)
        {
            box_idx_aux_x = ptr_ch->ptr_cell_idx_x[i] - ptr_ch->box_ts_x;
            box_idx_aux_y = ptr_ch->ptr_cell_idx_y[i] - ptr_ch->box_ts_y;
            box_idx_aux_z = ptr_ch->ptr_cell_idx_z[i] - ptr_ch->box_ts_z;
            box_idx_aux = box_idx_aux_x + box_idx_aux_y * ptr_ch->box_real_dim_x + box_idx_aux_z * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
            aux_mass2 += ptr_ch->ptr_box_mass[box_idx_aux];
        }

        printf("\n\n Error, particle updating Antes\n\n");
        printf("lv = %d, aux mass per box = %f\n", ptr_ch->lv, aux_mass1);
        printf("lv = %d, aux mass per cell = %f\n",ptr_ch->lv, aux_mass2);
        printf("lv = %d, local mass = %f\n",ptr_ch->lv, ptr_ch->local_mass);
        printf("lv = %d, no ptcl = %d\n",ptr_ch->lv,ptr_ch->ptcl_size);
        printf("\n\n");
    
    }










    if(lmin < lmax)
    {
        if (ptr_node->chn_size == 0)
        {

            for (int i = 0; i < no_ptcl; i++)
            {

                ptcl_idx = ptr_node->ptr_ptcl[i]; // Particle index

                //** >> First we remove the mass of the box_mass array of the node **/
                box_idx_node = ptcl_idx_to_box_idx(ptr_node, ptcl_idx);

                
                // if(i == no_ptcl -1)
                // {
                //     vtype aux_mass1 = 0;
                //     vtype aux_mass2 = 0;
                //     for (int i = 0; i < ptr_node->box_real_dim_x * ptr_node->box_real_dim_y * ptr_node->box_real_dim_z; i++)
                //     {
                //         aux_mass1 += ptr_node->ptr_box_mass[i];
                //     }
                //     aux_mass2 = 0;
                //     int box_idx_aux_x;
                //     int box_idx_aux_y;
                //     int box_idx_aux_z;
                //     int box_idx_aux;
                //     for (int i = 0; i < ptr_node->cell_size; i++)
                //     {
                //         box_idx_aux_x = ptr_node->ptr_cell_idx_x[i] - ptr_node->box_ts_x;
                //         box_idx_aux_y = ptr_node->ptr_cell_idx_y[i] - ptr_node->box_ts_y;
                //         box_idx_aux_z = ptr_node->ptr_cell_idx_z[i] - ptr_node->box_ts_z;
                //         box_idx_aux = box_idx_aux_x + box_idx_aux_y * ptr_node->box_real_dim_x + box_idx_aux_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                //         aux_mass2 += ptr_node->ptr_box_mass[box_idx_aux];
                //     }

                //     if (aux_mass1 != aux_mass2 || aux_mass1 != ptr_node->local_mass)
                //     {
                //         printf("\n\n Error \n\n");
                //         printf("aux mass per box = %f, ", aux_mass1);
                //         printf("aux mass per cell = %f, ", aux_mass2);
                //         printf("local mass = %f\n", ptr_node->local_mass);
                //         printf("\n\n");
                //     }
                // }




                //** >> Removing the mass of the mass box array of the node**/
                ptr_node->ptr_box_mass[box_idx_node] -= GL_ptcl_mass[ptcl_idx];



                if(ptr_node->ptr_box_mass[box_idx_node] < 0)
                {
                    printf("\nError (if), %d, ptcl = %d box cell mass = %f, local mass = %f\n", i, ptcl_idx, ptr_node->ptr_box_mass[box_idx_node], ptr_node->local_mass);
                }
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
                box_idx_node = ptcl_idx_to_box_idx(ptr_node, ptcl_idx);

                //** We ask if the particle leaves the node
                //** >> The particle moves towards its parent node or towards some sibling node  **/
                if (ptr_node->ptr_box[box_idx_node] < -3)
                {
                    //** The local mass is reduced **/
                    ptr_node->local_mass -= GL_ptcl_mass[ptcl_idx];

                    ptr_node_pt = ptr_node->ptr_pt;
                    //** >> Box index in the parent node **/
                    box_idx_pt = ptcl_idx_to_box_idx(ptr_node_pt, ptcl_idx);

                    //** >> If the particle moves towards a sibling node **/
                    if (ptr_node_pt->ptr_box[box_idx_pt] >= 0)
                    {
                        zone_idx = ptr_node_pt->ptr_box[box_idx_pt];
                        ptr_node_sib = ptr_node_pt->pptr_chn[zone_idx];

                        //** >> Space checking of the particle capacity in the sibling node **/
                        if (space_check(&(ptr_node_sib->ptcl_cap), ptr_node_sib->ptcl_size + 1, 2.0f, "p1i1", &(ptr_node_sib->ptr_ptcl)) == _FAILURE_)
                        {
                            printf("Error, in space_check function\n");
                            return _FAILURE_;
                        }

                        //** >> Adding the particle in the sibling node **/
                        ptr_node_sib->ptr_ptcl[ptr_node_sib->ptcl_size] = ptcl_idx;
                        ptr_node_sib->ptcl_size += 1; // +1 to the total particles in the sibling node

                        //** >> Adding the mass of the mass box array of the sibling node in the new cell position**/
                        box_idx_sib = ptcl_idx_to_box_idx(ptr_node_sib, ptcl_idx);
                        ptr_node_sib->ptr_box_mass[box_idx_sib] += GL_ptcl_mass[ptcl_idx];
                        //** The local mass is increased **/
                        ptr_node_sib->local_mass += GL_ptcl_mass[ptcl_idx];
                    }
                    //** If the particle is only in the parent node **/
                    else
                    {
                        //** >> Adding the mass of the mass box array of the parent node in the new cell position**/
                        ptr_node_pt->ptr_box_mass[box_idx_pt] += GL_ptcl_mass[ptcl_idx];
                    }

                    // Whether the particle stays at the parent node or moves to a sibling node, it must be removed from the current node
                    //** >> Removing the particle index of the current node ptr_node **/
                    // We move the last element of the array to the current position
                    ptr_node->ptr_ptcl[i] = ptr_node->ptr_ptcl[no_ptcl - 1];
                    no_ptcl--; // The total number of particle decrease
                    i--;       // The last element that was just moved to the current position should also must be analized
                }
                //** >> The particle stay in the node **/
                else
                {
                    //** >> Adding the mass of the mass box array of the node in the new cell position**/
                    ptr_node->ptr_box_mass[box_idx_node] += GL_ptcl_mass[ptcl_idx];
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

                    // if(i < no_ptcl -9800)
                    // {
                    //     vtype aux_mass1 = 0;
                    //     vtype aux_mass2 = 0;
                    //     for (int i = 0; i < ptr_node->box_real_dim_x * ptr_node->box_real_dim_y * ptr_node->box_real_dim_z; i++)
                    //     {
                    //         aux_mass1 += ptr_node->ptr_box_mass[i];
                    //     }
                    //     aux_mass2 = 0;
                    //     int box_idx_aux_x;
                    //     int box_idx_aux_y;
                    //     int box_idx_aux_z;
                    //     int box_idx_aux;
                    //     for (int i = 0; i < ptr_node->cell_size; i++)
                    //     {
                    //         box_idx_aux_x = ptr_node->ptr_cell_idx_x[i] - ptr_node->box_ts_x;
                    //         box_idx_aux_y = ptr_node->ptr_cell_idx_y[i] - ptr_node->box_ts_y;
                    //         box_idx_aux_z = ptr_node->ptr_cell_idx_z[i] - ptr_node->box_ts_z;
                    //         box_idx_aux = box_idx_aux_x + box_idx_aux_y * ptr_node->box_real_dim_x + box_idx_aux_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
                    //         aux_mass2 += ptr_node->ptr_box_mass[box_idx_aux];
                    //     }

                    //     if (aux_mass1 != aux_mass2 || aux_mass1 != ptr_node->local_mass)
                    //     {
                    //         printf("\n\n Error \n\n");
                    //         printf("aux mass per box = %f, ", aux_mass1);
                    //         printf("aux mass per cell = %f, ", aux_mass2);
                    //         printf("local mass = %f\n", ptr_node->local_mass);
                    //         printf("\n\n");
                    //     }
                    // }

                    //** >> Removing the mass of the mass box array of the node**/
                    box_idx_node = ptcl_idx_to_box_idx(ptr_node, ptcl_idx);
                    ptr_node->ptr_box_mass[box_idx_node] -= GL_ptcl_mass[ptcl_idx];

                    if (ptr_node->ptr_box_mass[box_idx_node] < 0)
                    {
                        printf("\nError (else), %d, ptcl = %d mass = %f\n", i, ptcl_idx, ptr_node->ptr_box_mass[box_idx_node]);
                    }

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

                    //** >> We ask if the particle stays on the node **/
                    box_idx_node = ptcl_idx_to_box_idx(ptr_node, ptcl_idx);

                    //** >> The particle moves towards one of its child nodes **/
                    if (ptr_node->ptr_box[box_idx_node] >= 0)
                    {
                        zone_idx = ptr_node->ptr_box[box_idx_node];
                        ptr_node_ch = ptr_node->pptr_chn[zone_idx];

                        //** >> Space checking of the particle capacity in the child node **/
                        if (space_check(&(ptr_node_ch->ptcl_cap), ptr_node_ch->ptcl_size + 1, 2.0f, "p1i1", &(ptr_node_ch->ptr_ptcl)) == _FAILURE_)
                        {
                            printf("Error, in space_check function\n");
                            return _FAILURE_;
                        }

                        //** >> Adding the particle in the sibling node **/
                        ptr_node_ch->ptr_ptcl[ptr_node_ch->ptcl_size] = ptcl_idx;
                        ptr_node_ch->ptcl_size += 1; // +1 to the total particles in the child node

                        //** >> Adding the mass of the mass box array of the chile node in the new cell position**/
                        box_idx_ch = ptcl_idx_to_box_idx(ptr_node_ch, ptcl_idx);
                        ptr_node_ch->ptr_box_mass[box_idx_ch] += GL_ptcl_mass[ptcl_idx];
                        //** The local mass is increased **/
                        ptr_node_ch->local_mass += GL_ptcl_mass[ptcl_idx];
                    }
                    //** >> The particle moves towards its parent node or towards some sibling node  **/
                    else if (ptr_node->ptr_box[box_idx_node] < -3)
                    {
                        //** The local mass is reduced **/
                        ptr_node->local_mass -= GL_ptcl_mass[ptcl_idx];

                        ptr_node_pt = ptr_node->ptr_pt;
                        //** >> Box index in the parent node **/
                        box_idx_pt = ptcl_idx_to_box_idx(ptr_node_pt, ptcl_idx);

                        //** >> If the particle moves towards a sibling node **/
                        if (ptr_node_pt->ptr_box[box_idx_pt] >= 0)
                        {
                            zone_idx = ptr_node_pt->ptr_box[box_idx_pt];
                            ptr_node_sib = ptr_node_pt->pptr_chn[zone_idx];

                            //** >> Space checking of the particle capacity in the sibling node **/
                            if (space_check(&(ptr_node_sib->ptcl_cap), ptr_node_sib->ptcl_size + 1, 2.0f, "p1i1", &(ptr_node_sib->ptr_ptcl)) == _FAILURE_)
                            {
                                printf("Error, in space_check function\n");
                                return _FAILURE_;
                            }

                            //** >> Adding the particle in the sibling node **/
                            ptr_node_sib->ptr_ptcl[ptr_node_sib->ptcl_size] = ptcl_idx;
                            ptr_node_sib->ptcl_size += 1; // +1 to the total particles in the sibling node

                            //** >> Adding the mass of the mass box array of the sibling node in the new cell position**/
                            box_idx_sib = ptcl_idx_to_box_idx(ptr_node_sib, ptcl_idx);
                            ptr_node_sib->ptr_box_mass[box_idx_sib] += GL_ptcl_mass[ptcl_idx];
                            //** The local mass is increased **/
                            ptr_node_sib->local_mass += GL_ptcl_mass[ptcl_idx];
                        }
                        //** >> The particle stay in the parent node node **/
                        else
                        {
                            //** >> Adding the mass of the mass box array of the parent node in the new cell position**/
                            ptr_node_pt->ptr_box_mass[box_idx_pt] += GL_ptcl_mass[ptcl_idx];
                        }

                        // Whether the particle stays at the parent node or moves to a sibling node, it must be removed from the current node
                        //** >> Removing the particle index of the current node ptr_node **/
                        // We move the last element of the array to the current position
                        ptr_node->ptr_ptcl[i] = ptr_node->ptr_ptcl[no_ptcl - 1];
                        no_ptcl--; // The total number of particle decrease
                        i--;       // The last element that was just moved to the current position should also must be analized
                    }
                    //** >> The particle stay in the node but not inside of a child **/
                    else
                    {
                        //** >> Adding the mass of the mass box array of the node in the new cell position**/
                        ptr_node->ptr_box_mass[box_idx_node] += GL_ptcl_mass[ptcl_idx];
                    }

                    //** >> The status of the particle is changed from not updated to updated **/
                    GL_ptcl_updating_flag[ptcl_idx] = status;
                }
            } // End cycle over particles in the node
        }

        //** >> Updating the number of particles in the current node **/
        ptr_node->ptcl_size = no_ptcl;
    }   //End if(lmin < lmax)
    else // Case only 1 level, i.e. lmin = lmax
    {
        for (int i = 0; i < no_ptcl; i++)
        {
            ptcl_idx = ptr_node->ptr_ptcl[i]; // Particle index

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
        }
    }




    ptr_node_ch = NULL;
    ptr_node_pt = NULL; 
    ptr_node_sib = NULL;







    aux_mass1 = 0;
    for (int i = 0; i < ptr_node->box_real_dim_x * ptr_node->box_real_dim_y * ptr_node->box_real_dim_z; i++)
    {
        aux_mass1 += ptr_node->ptr_box_mass[i];
    }
    aux_mass2 = 0;
    for (int i = 0; i < ptr_node->cell_size; i++)
    {
        box_idx_aux_x = ptr_node->ptr_cell_idx_x[i] - ptr_node->box_ts_x;
        box_idx_aux_y = ptr_node->ptr_cell_idx_y[i] - ptr_node->box_ts_y;
        box_idx_aux_z = ptr_node->ptr_cell_idx_z[i] - ptr_node->box_ts_z;
        box_idx_aux = box_idx_aux_x + box_idx_aux_y * ptr_node->box_real_dim_x + box_idx_aux_z * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
        aux_mass2 += ptr_node->ptr_box_mass[box_idx_aux];
    }

    printf("\n\n Error, particle updating Despues\n\n");
    printf("lv = %d, aux mass per box = %f\n", ptr_node->lv, aux_mass1);
    printf("lv = %d, aux mass per cell = %f\n",ptr_node->lv, aux_mass2);
    printf("lv = %d, local mass = %f\n",ptr_node->lv, ptr_node->local_mass);
    printf("lv = %d, no ptcl = %d\n",ptr_node->lv,ptr_node->ptcl_size);
    printf("\n\n");



    if(ptr_node->lv == 5)
    {
        struct node *ptr_ch;
        ptr_ch = ptr_node->pptr_chn[0];
        aux_mass1 = 0;
        for (int i = 0; i < ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y * ptr_ch->box_real_dim_z; i++)
        {
            aux_mass1 += ptr_ch->ptr_box_mass[i];
        }
        aux_mass2 = 0;
        for (int i = 0; i < ptr_ch->cell_size; i++)
        {
            box_idx_aux_x = ptr_ch->ptr_cell_idx_x[i] - ptr_ch->box_ts_x;
            box_idx_aux_y = ptr_ch->ptr_cell_idx_y[i] - ptr_ch->box_ts_y;
            box_idx_aux_z = ptr_ch->ptr_cell_idx_z[i] - ptr_ch->box_ts_z;
            box_idx_aux = box_idx_aux_x + box_idx_aux_y * ptr_ch->box_real_dim_x + box_idx_aux_z * ptr_ch->box_real_dim_x * ptr_ch->box_real_dim_y;
            aux_mass2 += ptr_ch->ptr_box_mass[box_idx_aux];
        }

        printf("\n\n Error, particle updating Despues\n\n");
        printf("lv = %d, aux mass per box = %f\n", ptr_ch->lv, aux_mass1);
        printf("lv = %d, aux mass per cell = %f\n",ptr_ch->lv, aux_mass2);
        printf("lv = %d, local mass = %f\n",ptr_ch->lv, ptr_ch->local_mass);
        printf("lv = %d, no ptcl = %d\n",ptr_ch->lv,ptr_ch->ptcl_size);
        printf("\n\n");
    }

    return _SUCCESS_;
}

int particle_updating_A(vtype dt)
{

    // Position and velocity are updated first by a "Predictor step"

    //** >> Particle updating A **/

    struct node *ptr_node = NULL;
    bool status;    // Boolean value for the updating particles

    int no_pts; // Number of parents in the cycle
    status = !GL_ptcl_updating_flag[0];

    

    for (int lv = GL_tentacles_level_max; lv > -1; lv--)
    {
        no_pts = GL_tentacles_size[lv];
        //** >> For cycle over parent nodes **/
        for (int i = 0; i < no_pts; i++)
        {
            ptr_node = GL_tentacles[lv][i];

            if (computing_particles_updating_A(ptr_node, dt, status) == _FAILURE_)
            {
                printf("Error at function computing_particles_updating_A()\n");
                return _FAILURE_;
            }
        }
    }

    ptr_node = NULL;

    return _SUCCESS_;
}