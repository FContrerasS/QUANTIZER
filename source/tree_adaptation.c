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

    no_chn = ptr_node->chn_size;

    if (no_chn > 0)
    {

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