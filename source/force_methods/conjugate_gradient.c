/*
 * conjugate_gradient.c
 *
 * Solve the Discrete Poisson equation using the conjugate gradient method
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

#include "conjugate_gradient.h"

int conjugate_gradient(struct node *ptr_node)
{

    clock_t aux_clock;

    int box_grid_idx;

    vtype H = 1.0L / (1 << ptr_node->lv);
    vtype one_over_H_pow_2 = 1.0L / (H * H);
    int size = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1) * (ptr_node->box_real_dim_z + 1);

    vtype *r, *p, *Ap;
    r = (vtype *)calloc(size, sizeof(vtype));
    p = (vtype *)malloc(size * sizeof(vtype));
    Ap = (vtype *)malloc(size * sizeof(vtype));

    vtype alpha, rsold, rsnew;

    int grid_box_real_dim_X = ptr_node->box_real_dim_x + 1;
    int grid_box_real_dim_X_times_Y = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

    // Computing r
    for (int i = 0; i < ptr_node->grid_intr_size; i++)
    {
        box_grid_idx = ptr_node->ptr_intr_grid_idx[i];
        r[box_grid_idx] = ptr_node->ptr_d[box_grid_idx] + one_over_H_pow_2 *
                                                              (6.0 * ptr_node->ptr_pot[box_grid_idx] - ptr_node->ptr_pot[box_grid_idx + 1] - ptr_node->ptr_pot[box_grid_idx - 1] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X] - ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X] - ptr_node->ptr_pot[box_grid_idx + grid_box_real_dim_X_times_Y] - ptr_node->ptr_pot[box_grid_idx - grid_box_real_dim_X_times_Y]);
    }

    // Copy r to p
    memcpy(p, r, size * sizeof(vtype));

    // Computing rsold
    rsold = 0;
    for (int i = 0; i < ptr_node->grid_intr_size; i++)
    {
        box_grid_idx = ptr_node->ptr_intr_grid_idx[i];
        rsold += r[box_grid_idx] * r[box_grid_idx];
    }

    // //** >> CHEKING ERROR SOLUTION CONDITION **/
    // aux_clock = clock();
    // if (poisson_error(ptr_node, r, rsold, 1) == true)
    // {
    //     GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
    //     free(r);
    //     free(p);
    //     free(Ap);
    //     return _SUCCESS_;
    // }
    // GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

    int cycle_size = size < _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_ ? size : _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_;
    // For cycle over mutually conjugate vectors p
    for (int i = 0; i < cycle_size; i++)
    {

        // Computing Ap
        for (int j = 0; j < ptr_node->grid_intr_size; j++)
        {
            box_grid_idx = ptr_node->ptr_intr_grid_idx[j];
            Ap[box_grid_idx] = -one_over_H_pow_2 * (6.0 * p[box_grid_idx] - p[box_grid_idx + 1] - p[box_grid_idx - 1] - p[box_grid_idx + grid_box_real_dim_X] - p[box_grid_idx - grid_box_real_dim_X] - p[box_grid_idx + grid_box_real_dim_X_times_Y] - p[box_grid_idx - grid_box_real_dim_X_times_Y]);
            //printf("Ap = %1.3Le\n", Ap[box_grid_idx]);
        }

        // Computing alpha
        alpha = 0;
        for (int j = 0; j < ptr_node->grid_intr_size; j++)
        {
            box_grid_idx = ptr_node->ptr_intr_grid_idx[j];
            alpha += p[box_grid_idx] * Ap[box_grid_idx];
        }

        alpha = rsold / alpha;

        // New solution
        for (int j = 0; j < ptr_node->grid_intr_size; j++)
        {
            box_grid_idx = ptr_node->ptr_intr_grid_idx[j];
            ptr_node->ptr_pot[box_grid_idx] += alpha * p[box_grid_idx];
        }

        // New rest
        for (int j = 0; j < ptr_node->grid_intr_size; j++)
        {
            box_grid_idx = ptr_node->ptr_intr_grid_idx[j];
            r[box_grid_idx] -= alpha * Ap[box_grid_idx];
        }

        // Computing rsnew
        rsnew = 0;
        for (int j = 0; j < ptr_node->grid_intr_size; j++)
        {
            box_grid_idx = ptr_node->ptr_intr_grid_idx[j];
            rsnew += r[box_grid_idx] * r[box_grid_idx];
        }

        //** >> CHEKING ERROR SOLUTION CONDITION **/
        if (i % 25 == 0)
        {
            aux_clock = clock();
            if (poisson_error(ptr_node, r, rsnew, 1) == true)
            {
                //printf("\ni = %d\n", i);
                GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
                free(r);
                free(p);
                free(Ap);
                return _SUCCESS_;
            }
            GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
        }
        //printf("\ni = %d, rsnew = %1.3Le, alpha = %1.3Le\n", i, rsnew,alpha);
        // Updating p
        for (int j = 0; j < ptr_node->grid_intr_size; j++)
        {
            box_grid_idx = ptr_node->ptr_intr_grid_idx[j];
            p[box_grid_idx] = r[box_grid_idx] + rsnew / rsold * p[box_grid_idx];
        }
        // Updating rsold
        rsold = rsnew;
    }

    

    free(r);
    free(p);
    free(Ap);
    return _FAILURE_;
}
