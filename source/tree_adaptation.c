/*
 * tree_adaptation.c
 *
 * Adapting the tree to the next time-step iteration
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

#include "tree_adaptation.h"

static int updating_box_mass(struct node *ptr_node)
{

    int no_chn; // Number of child nodes
    struct node *ptr_node_ch; //child node
    vtype aux_box_mass; //Total mass in the 8 child cells

    int box_idx_x_ch; // Box index in X direcction of the child cell 
    int box_idx_y_ch; // Box index in Y direcction of the child cell
    int box_idx_z_ch; // Box index in Z direcction of the child cell
    int box_idx_ch;  // Box index of the child cell

    int box_idx_x_node;   // Box index in X direcction of the node cell
    int box_idx_y_node;   // Box index in Y direcction of the node cell
    int box_idx_z_node;   // Box index in Z direcction of the node cell
    int box_idx_node;  // Box index of the node cell

    no_chn = ptr_node->chn_size;

    for (int i = 0; i < no_chn; i++)    //Cycle over children
    {
        ptr_node_ch = ptr_node->pptr_chn[i];
        for (int j = 0; j < ptr_node_ch->cell_size; j += 8) // Cycle over packeges of 8 cells
        {
            //** >> Computing the mass of the 8 child cells **/
            box_idx_x_ch = ptr_node_ch->ptr_cell_idx_x[j] - ptr_node_ch->box_ts_x;
            box_idx_y_ch = ptr_node_ch->ptr_cell_idx_y[j] - ptr_node_ch->box_ts_y;
            box_idx_z_ch = ptr_node_ch->ptr_cell_idx_z[j] - ptr_node_ch->box_ts_z;

            aux_box_mass = 0;

            for (int kk = 0; kk < 2; kk++)
            {
                for (int jj = 0; jj < 2; jj++)
                {
                    for (int ii = 0; ii < 2; ii++)
                    {
                        box_idx_ch = (box_idx_x_ch + ii) + (box_idx_y_ch + jj) * ptr_node_ch->box_real_dim_x + (box_idx_z_ch + kk) * ptr_node_ch->box_real_dim_x * ptr_node_ch->box_real_dim_y;
                        aux_box_mass += ptr_node_ch->ptr_box_mass[box_idx_ch];
                    }
                }
            }
            //** Re-defining the mass of the node cell corresponding to those 8 child cells
            box_idx_x_node = (ptr_node_ch->ptr_cell_idx_x[j] >> 1) - ptr_node->box_ts_x;
            box_idx_y_node = (ptr_node_ch->ptr_cell_idx_y[j] >> 1) - ptr_node->box_ts_y;
            box_idx_z_node = (ptr_node_ch->ptr_cell_idx_z[j] >> 1) - ptr_node->box_ts_z;
            box_idx_node = box_idx_x_node + box_idx_y_node * ptr_node->box_real_dim_x + box_idx_z_node * ptr_node->box_real_dim_x * ptr_node->box_real_dim_y;
            ptr_node->ptr_box_mass[box_idx_node] = aux_box_mass;
            total_mass += aux_box_mass;
        }
    }



    return _SUCCESS_;
}

static int node_adaptation(struct node *ptr_node)
{

    return _SUCCESS_;
}

int tree_adaptation()
{

    //** >> Working in the refinement zones **/
    if (lmin < lmax)
    {
        struct node *ptr_node;
        int no_pts; // Number of parents in the cycle
        int no_lvs; // Number of level of refinement to adapt

        no_lvs = GL_tentacles_level_max < lmax - lmin ? GL_tentacles_level_max : GL_tentacles_level_max - 1;

        for (int lv = no_lvs; lv > -1; lv--)
        {
            no_pts = GL_tentacles_size[lv];
            //** >> For cycle over parent nodes **/
            for (int i = 0; i < no_pts; i++)
            {
                ptr_node = GL_tentacles_old[lv][i];

                //** Updating the box mass information **/
                if (updating_box_mass(ptr_node) == _FAILURE_)
                {
                    printf("Error at function node_adaptation()\n");
                    return _FAILURE_;
                }

                if (node_adaptation(ptr_node) == _FAILURE_)
                {
                    printf("Error at function node_adaptation()\n");
                    return _FAILURE_;
                }
            }
        }
    }

        return _SUCCESS_;
}