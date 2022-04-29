/*
 * potential.c
 *
 * Potential computation in each node
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

#include "potential.h"


static int fill_red_and_black(int **pptr_red_black, int *ptr_red_black_cap, int *ptr_red_black_size, const struct node *ptr_node)
{
    int red_size;   // Number of red points
    int black_size; // Number of black points
    int red_or_black;   // Desition if the point is red or black

    int box_idx_x; // Box index in X direcction
    int box_idx_y; // Box index in Y direcction
    int box_idx_z; // Box index in Z direcction

    int box_grid_idx; // box grid index

    red_size = 0;
    black_size = 0;

    //** >> Filling Red and Black arrays with the correspondig grid points **/
    for (int i = 0; i < ptr_node->grid_intr_size; i++)
    {
        box_grid_idx = ptr_node->ptr_grid_intr[i];
        box_idx_z = box_grid_idx / ((ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1));
        box_idx_y = (box_grid_idx - box_idx_z * (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1)) / (ptr_node->box_real_dim_x + 1);
        box_idx_x = box_grid_idx - box_idx_z * (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1) - box_idx_y * (ptr_node->box_real_dim_x + 1);
        red_or_black = box_idx_x + box_idx_y + box_idx_z;

        //** >> Space checking of red and black pointers **/
        if (red_or_black % 2 == 0)   // If it is red
        {
            red_size++;

            //** >> Space checking of the particle capacity of the red array **/
            if (space_check(&(ptr_red_black_cap[0]), red_size, 2.0f, "p1i1", &(pptr_red_black[0])) == _FAILURE_)
            {
                printf("Error, in space_check function\n");
                return _FAILURE_;
            }

            //** >> Adding element to the red array **/
            pptr_red_black[0][red_size - 1] = box_grid_idx;
        }   
        else // If it is black
        {
            black_size++;

            //** >> Space checking of the particle capacity of the black array **/
            if (space_check(&(ptr_red_black_cap[1]), black_size, 2.0f, "p1i1", &(pptr_red_black[1])) == _FAILURE_)
            {
                printf("Error, in space_check function\n");
                return _FAILURE_;
            }

            //** >> Adding element to the black array **/
            pptr_red_black[1][black_size - 1] = box_grid_idx;
        }
    }

    //** >> Size of the red and black arrays **/
    ptr_red_black_size[0] = red_size;
    ptr_red_black_size[1] = black_size;

    return _SUCCESS_;
}

static void initial_potential(const struct node *ptr_node_pt, struct node *ptr_node_ch)
{
    int ch_box_grid_idx; // Child box grid index
    int ch_box_idx_x;   // Child box index
    int ch_box_idx_y;
    int ch_box_idx_z;


    int ch_cell_idx_x;  // Chilld cell index at x direcction
    int ch_cell_idx_y;
    int ch_cell_idx_z;

    int pt_box_idx_x0; // Parent box index at X direcction
    int pt_box_idx_x1;
    int pt_box_idx_y0;
    int pt_box_idx_y1;
    int pt_box_idx_z0;
    int pt_box_idx_z1;

    int pt_box_grid_i0_j0_k0; // Parent grid index at i0, j0, k0 position
    int pt_box_grid_i1_j0_k0;
    int pt_box_grid_i0_j1_k0;
    int pt_box_grid_i1_j1_k0;
    int pt_box_grid_i0_j0_k1;
    int pt_box_grid_i1_j0_k1;
    int pt_box_grid_i0_j1_k1;
    int pt_box_grid_i1_j1_k1;

    vtype aux_pot; // Auxiliar Potential in the grid point

    //** >> Passing the potential from coarse parent to fine child **/
    // Border grid points
    for (int i = 0; i < ptr_node_ch->grid_bder_size; i++)
    {
        //** >> Child box indexes **/
        ch_box_grid_idx = ptr_node_ch->ptr_grid_bder[i];
        ch_box_idx_z = ch_box_grid_idx / ((ptr_node_ch->box_real_dim_x + 1) * (ptr_node_ch->box_real_dim_y + 1));
        ch_box_idx_y = (ch_box_grid_idx - ch_box_idx_z * (ptr_node_ch->box_real_dim_x + 1) * (ptr_node_ch->box_real_dim_y + 1)) / (ptr_node_ch->box_real_dim_x + 1);
        ch_box_idx_x = ch_box_grid_idx - ch_box_idx_z * (ptr_node_ch->box_real_dim_x + 1) * (ptr_node_ch->box_real_dim_y + 1) - ch_box_idx_y * (ptr_node_ch->box_real_dim_x + 1);

        //** >> Child cell indexes **/
        ch_cell_idx_x = ch_box_idx_x + ptr_node_ch->box_ts_x;
        ch_cell_idx_y = ch_box_idx_y + ptr_node_ch->box_ts_y;
        ch_cell_idx_z = ch_box_idx_z + ptr_node_ch->box_ts_z;

        //** >> Parent box indexes **/
        pt_box_idx_x0 = (ch_cell_idx_x >> 1) - ptr_node_pt->box_ts_x;
        pt_box_idx_x1 = ((ch_cell_idx_x + 1) >> 1) - ptr_node_pt->box_ts_x;
        pt_box_idx_y0 = (ch_cell_idx_y >> 1) - ptr_node_pt->box_ts_y;
        pt_box_idx_y1 = ((ch_cell_idx_y + 1) >> 1) - ptr_node_pt->box_ts_y;
        pt_box_idx_z0 = (ch_cell_idx_z >> 1) - ptr_node_pt->box_ts_z;
        pt_box_idx_z1 = ((ch_cell_idx_z + 1) >> 1) - ptr_node_pt->box_ts_z;

        //** >> Parent grid indexes **/
        pt_box_grid_i0_j0_k0 = pt_box_idx_x0 + pt_box_idx_y0 * (ptr_node_pt->box_real_dim_x + 1) + pt_box_idx_z0 * (ptr_node_pt->box_real_dim_x + 1) * (ptr_node_pt->box_real_dim_y + 1);
        pt_box_grid_i1_j0_k0 = pt_box_idx_x1 + pt_box_idx_y0 * (ptr_node_pt->box_real_dim_x + 1) + pt_box_idx_z0 * (ptr_node_pt->box_real_dim_x + 1) * (ptr_node_pt->box_real_dim_y + 1);
        pt_box_grid_i0_j1_k0 = pt_box_idx_x0 + pt_box_idx_y1 * (ptr_node_pt->box_real_dim_x + 1) + pt_box_idx_z0 * (ptr_node_pt->box_real_dim_x + 1) * (ptr_node_pt->box_real_dim_y + 1);
        pt_box_grid_i1_j1_k0 = pt_box_idx_x1 + pt_box_idx_y1 * (ptr_node_pt->box_real_dim_x + 1) + pt_box_idx_z0 * (ptr_node_pt->box_real_dim_x + 1) * (ptr_node_pt->box_real_dim_y + 1);
        pt_box_grid_i0_j0_k1 = pt_box_idx_x0 + pt_box_idx_y0 * (ptr_node_pt->box_real_dim_x + 1) + pt_box_idx_z1 * (ptr_node_pt->box_real_dim_x + 1) * (ptr_node_pt->box_real_dim_y + 1);
        pt_box_grid_i1_j0_k1 = pt_box_idx_x1 + pt_box_idx_y0 * (ptr_node_pt->box_real_dim_x + 1) + pt_box_idx_z1 * (ptr_node_pt->box_real_dim_x + 1) * (ptr_node_pt->box_real_dim_y + 1);
        pt_box_grid_i0_j1_k1 = pt_box_idx_x0 + pt_box_idx_y1 * (ptr_node_pt->box_real_dim_x + 1) + pt_box_idx_z1 * (ptr_node_pt->box_real_dim_x + 1) * (ptr_node_pt->box_real_dim_y + 1);
        pt_box_grid_i1_j1_k1 = pt_box_idx_x1 + pt_box_idx_y1 * (ptr_node_pt->box_real_dim_x + 1) + pt_box_idx_z1 * (ptr_node_pt->box_real_dim_x + 1) * (ptr_node_pt->box_real_dim_y + 1);

        aux_pot = 0.125 * (ptr_node_pt->ptr_pot[pt_box_grid_i0_j0_k0] +
                           ptr_node_pt->ptr_pot[pt_box_grid_i1_j0_k0] +
                           ptr_node_pt->ptr_pot[pt_box_grid_i0_j1_k0] +
                           ptr_node_pt->ptr_pot[pt_box_grid_i1_j1_k0] +
                           ptr_node_pt->ptr_pot[pt_box_grid_i0_j0_k1] +
                           ptr_node_pt->ptr_pot[pt_box_grid_i1_j0_k1] +
                           ptr_node_pt->ptr_pot[pt_box_grid_i0_j1_k1] +
                           ptr_node_pt->ptr_pot[pt_box_grid_i1_j1_k1]);

        ptr_node_ch->ptr_pot[ch_box_grid_idx] = aux_pot;
    }

    // Interior grid points
    for (int i = 0; i < ptr_node_ch->grid_intr_size; i++)
    {
        //** >> Child box indexes **/
        ch_box_grid_idx = ptr_node_ch->ptr_grid_intr[i];
        ch_box_idx_z = ch_box_grid_idx / ((ptr_node_ch->box_real_dim_x + 1) * (ptr_node_ch->box_real_dim_y + 1));
        ch_box_idx_y = (ch_box_grid_idx - ch_box_idx_z * (ptr_node_ch->box_real_dim_x + 1) * (ptr_node_ch->box_real_dim_y + 1)) / (ptr_node_ch->box_real_dim_x + 1);
        ch_box_idx_x = ch_box_grid_idx - ch_box_idx_z * (ptr_node_ch->box_real_dim_x + 1) * (ptr_node_ch->box_real_dim_y + 1) - ch_box_idx_y * (ptr_node_ch->box_real_dim_x + 1);

        //** >> Child cell indexes **/
        ch_cell_idx_x = ch_box_idx_x + ptr_node_ch->box_ts_x;
        ch_cell_idx_y = ch_box_idx_y + ptr_node_ch->box_ts_y;
        ch_cell_idx_z = ch_box_idx_z + ptr_node_ch->box_ts_z;

        //** >> Parent box indexes **/
        pt_box_idx_x0 = (ch_cell_idx_x >> 1) - ptr_node_pt->box_ts_x;
        pt_box_idx_x1 = ((ch_cell_idx_x + 1) >> 1) - ptr_node_pt->box_ts_x;
        pt_box_idx_y0 = (ch_cell_idx_y >> 1) - ptr_node_pt->box_ts_y;
        pt_box_idx_y1 = ((ch_cell_idx_y + 1) >> 1) - ptr_node_pt->box_ts_y;
        pt_box_idx_z0 = (ch_cell_idx_z >> 1) - ptr_node_pt->box_ts_z;
        pt_box_idx_z1 = ((ch_cell_idx_z + 1) >> 1) - ptr_node_pt->box_ts_z;

        //** >> Parent grid indexes **/
        pt_box_grid_i0_j0_k0 = pt_box_idx_x0 + pt_box_idx_y0 * (ptr_node_pt->box_real_dim_x + 1) + pt_box_idx_z0 * (ptr_node_pt->box_real_dim_x + 1) * (ptr_node_pt->box_real_dim_y + 1);
        pt_box_grid_i1_j0_k0 = pt_box_idx_x1 + pt_box_idx_y0 * (ptr_node_pt->box_real_dim_x + 1) + pt_box_idx_z0 * (ptr_node_pt->box_real_dim_x + 1) * (ptr_node_pt->box_real_dim_y + 1);
        pt_box_grid_i0_j1_k0 = pt_box_idx_x0 + pt_box_idx_y1 * (ptr_node_pt->box_real_dim_x + 1) + pt_box_idx_z0 * (ptr_node_pt->box_real_dim_x + 1) * (ptr_node_pt->box_real_dim_y + 1);
        pt_box_grid_i1_j1_k0 = pt_box_idx_x1 + pt_box_idx_y1 * (ptr_node_pt->box_real_dim_x + 1) + pt_box_idx_z0 * (ptr_node_pt->box_real_dim_x + 1) * (ptr_node_pt->box_real_dim_y + 1);
        pt_box_grid_i0_j0_k1 = pt_box_idx_x0 + pt_box_idx_y0 * (ptr_node_pt->box_real_dim_x + 1) + pt_box_idx_z1 * (ptr_node_pt->box_real_dim_x + 1) * (ptr_node_pt->box_real_dim_y + 1);
        pt_box_grid_i1_j0_k1 = pt_box_idx_x1 + pt_box_idx_y0 * (ptr_node_pt->box_real_dim_x + 1) + pt_box_idx_z1 * (ptr_node_pt->box_real_dim_x + 1) * (ptr_node_pt->box_real_dim_y + 1);
        pt_box_grid_i0_j1_k1 = pt_box_idx_x0 + pt_box_idx_y1 * (ptr_node_pt->box_real_dim_x + 1) + pt_box_idx_z1 * (ptr_node_pt->box_real_dim_x + 1) * (ptr_node_pt->box_real_dim_y + 1);
        pt_box_grid_i1_j1_k1 = pt_box_idx_x1 + pt_box_idx_y1 * (ptr_node_pt->box_real_dim_x + 1) + pt_box_idx_z1 * (ptr_node_pt->box_real_dim_x + 1) * (ptr_node_pt->box_real_dim_y + 1);

        aux_pot = 0.125 * (ptr_node_pt->ptr_pot[pt_box_grid_i0_j0_k0] +
                           ptr_node_pt->ptr_pot[pt_box_grid_i1_j0_k0] +
                           ptr_node_pt->ptr_pot[pt_box_grid_i0_j1_k0] +
                           ptr_node_pt->ptr_pot[pt_box_grid_i1_j1_k0] +
                           ptr_node_pt->ptr_pot[pt_box_grid_i0_j0_k1] +
                           ptr_node_pt->ptr_pot[pt_box_grid_i1_j0_k1] +
                           ptr_node_pt->ptr_pot[pt_box_grid_i0_j1_k1] +
                           ptr_node_pt->ptr_pot[pt_box_grid_i1_j1_k1]);

        ptr_node_ch->ptr_pot[ch_box_grid_idx] = aux_pot;
    }
}

static int potential_branch_node(int **pptr_red_black, const int *ptr_red_black_size, struct node *ptr_node)
{
    vtype H;
    vtype H_pow2;

    int red_size;
    int black_size;
    vtype aux_pot;

    red_size = ptr_red_black_size[0];
    black_size = ptr_red_black_size[1];

    int box_grid_idx;
    int box_grid_idxNbr_x_plus;  // Box grid index in the neigborhood on the right
    int box_grid_idxNbr_x_minus; // Box grid index in the neigborhood on the left
    int box_grid_idxNbr_y_plus;  // Box grid index in the neigborhood behind
    int box_grid_idxNbr_y_minus; // Box grid index in the neigborhood in front
    int box_grid_idxNbr_z_plus;  // Box grid index in the neigborhood up
    int box_grid_idxNbr_z_minus; // Box grid index in the neigborhood down

    H = 1.0L / (1 << ptr_node->lv);
    H_pow2 = H * H;

    //** >> Cycle over the Successive over-relaxation **/
    for (int iter = 0; iter < _Iter_branches_SOR_; iter++)
    {
        //** >> Cycle over red points **/
        for (int i = 0; i < red_size; i++)
        {
            box_grid_idx = pptr_red_black[0][i];
            box_grid_idxNbr_x_plus = box_grid_idx + 1;
            box_grid_idxNbr_x_minus = box_grid_idx - 1;
            box_grid_idxNbr_y_plus = box_grid_idx + (ptr_node->box_real_dim_x + 1);
            box_grid_idxNbr_y_minus = box_grid_idx - (ptr_node->box_real_dim_x + 1);
            box_grid_idxNbr_z_plus = box_grid_idx + (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);
            box_grid_idxNbr_z_minus = box_grid_idx - (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

            aux_pot = 1.0L / 6.0L * (-H_pow2 * ptr_node->ptr_d[box_grid_idx] + ptr_node->ptr_pot[box_grid_idxNbr_x_plus] + ptr_node->ptr_pot[box_grid_idxNbr_x_minus] + ptr_node->ptr_pot[box_grid_idxNbr_y_plus] + ptr_node->ptr_pot[box_grid_idxNbr_y_minus] + ptr_node->ptr_pot[box_grid_idxNbr_z_plus] + ptr_node->ptr_pot[box_grid_idxNbr_z_minus]);
            ptr_node->ptr_pot[box_grid_idx] = _w_SOR_ * aux_pot + (1 - _w_SOR_) * ptr_node->ptr_pot[box_grid_idx];
        }

        //** >> Cycle over black points **/
        for (int i = 0; i < black_size; i++)
        {
            box_grid_idx = pptr_red_black[1][i];
            box_grid_idxNbr_x_plus = box_grid_idx + 1;
            box_grid_idxNbr_x_minus = box_grid_idx - 1;
            box_grid_idxNbr_y_plus = box_grid_idx + (ptr_node->box_real_dim_x + 1);
            box_grid_idxNbr_y_minus = box_grid_idx - (ptr_node->box_real_dim_x + 1);
            box_grid_idxNbr_z_plus = box_grid_idx + (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);
            box_grid_idxNbr_z_minus = box_grid_idx - (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

            aux_pot = 1.0L / 6.0L * (-H_pow2 * ptr_node->ptr_d[box_grid_idx] + ptr_node->ptr_pot[box_grid_idxNbr_x_plus] + ptr_node->ptr_pot[box_grid_idxNbr_x_minus] + ptr_node->ptr_pot[box_grid_idxNbr_y_plus] + ptr_node->ptr_pot[box_grid_idxNbr_y_minus] + ptr_node->ptr_pot[box_grid_idxNbr_z_plus] + ptr_node->ptr_pot[box_grid_idxNbr_z_minus]);
            ptr_node->ptr_pot[box_grid_idx] = _w_SOR_ * aux_pot + (1 - _w_SOR_) * ptr_node->ptr_pot[box_grid_idx];
        }
    }

    return _SUCCESS_;
}

int potential()
{
    int iter;
    bool check;

    clock_t aux_clock;


    //**>> POTENTIAL HEAD NODE **/

    //** >> CHEKING ERROR SOLUTION CONDITION **/
    aux_clock = clock();
    check = poisson_error(GL_ptr_tree);
    GL_times[21] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

    //** >> SOLVING POISSON EQUATION **/
    //printf("potential Head\n");
    iter = 0;
    while (iter < _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_ && check == false)
    {
        iter = iter + 1;

        //** >> INITIAL POTENTIAL COMPUTATION **/
        aux_clock = clock();
        if (potential_head_node() == _FAILURE_)
        {
            printf("\n\n Error running potential_head_node() function \n\n ");
            return _FAILURE_;
        }
        GL_times[20] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

        //** >> CHEKING ERROR SOLUTION CONDITION **/
        aux_clock = clock();
        check = poisson_error(GL_ptr_tree);
        GL_times[21] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
    }

    if (iter == _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_)
    {
        printf("\nERROR: The precision was not reached in the parent node. Too many Multigrid Iterations, plz choose a lower precision than %1.3e\n", _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_);
        return _FAILURE_;
    }

    //printf("Potential in Branches\n");
    //** >> POTENTIAL BRANCH NODES **/
    if(lmin < lmax)
    {
        struct node *ptr_node_pt = NULL;
        struct node *ptr_node_ch = NULL;

        int no_pts;   // Number of parents in the cycle

        int **pptr_red_black = NULL;
        int *ptr_red_black_cap = NULL;
        int *ptr_red_black_size = NULL;

        // Red and Black arrays contain the index element of the box grid of the potential
        pptr_red_black = (int **)malloc(2 * sizeof(int *));
        pptr_red_black[0] = NULL;
        pptr_red_black[1] = NULL;
        ptr_red_black_cap = (int *)malloc(2 * sizeof(int));
        ptr_red_black_cap[0] = 0;
        ptr_red_black_cap[1] = 0;
        ptr_red_black_size = (int *)malloc(2 * sizeof(int));
        ptr_red_black_size[0] = 0;
        ptr_red_black_size[1] = 0;

        for (int lv = 0; lv < lmax - lmin + 1; lv++)
        {
            no_pts = GL_tentacles_size[lv];

            //** >> For cycle over parent nodes **/
            for (int i = 0; i < no_pts; i++)
            {
                ptr_node_pt = GL_tentacles[lv][i];
                for (int j = 0; j < ptr_node_pt->chn_size; j++)
                {
                    ptr_node_ch = ptr_node_pt->pptr_chn[j];
                    //printf("Parent lv = %d, Parent ID = %d, child ID = %d\n",ptr_node_pt->lv, ptr_node_pt->ID,ptr_node_ch->ID);
                    //** >> Transfer the potential from parent to child nodes **/
                    aux_clock = clock();
                    initial_potential(ptr_node_pt, ptr_node_ch);
                    GL_times[22] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                    //** >> Computing the potential in the child node **/
                    aux_clock = clock();
                    if (fill_red_and_black(pptr_red_black, ptr_red_black_cap, ptr_red_black_size, ptr_node_ch) == _FAILURE_)
                    {
                        printf("Error at function fill_red_and_black()\n");
                        return _FAILURE_;
                    }
                    GL_times[23] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                    //** >> CHEKING ERROR SOLUTION CONDITION **/
                    aux_clock = clock();
                    check = poisson_error(ptr_node_ch);
                    GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                    iter = 0;
                    while (iter < _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_ && check == false)
                    {
                        iter = iter + 1;

                        //** >> Computing the potential in the child node **/
                        aux_clock = clock();
                        if (potential_branch_node(pptr_red_black, ptr_red_black_size, ptr_node_ch) == _FAILURE_)
                        {
                            printf("Error at function potential_branch_node()\n");
                            return _FAILURE_;
                        }
                        GL_times[24] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

                        //** >> CHEKING ERROR SOLUTION CONDITION **/
                        aux_clock = clock();
                        check = poisson_error(ptr_node_ch);
                        GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
                    }
                    if (iter == _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_)
                    {
                        printf("\nERROR: The precision was not reached in the branch node. Too many SOR Iterations, plz choose a lower precision than %1.3e\n", _ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_);
                        return _FAILURE_;
                    }
                }   // End cycle over number of children
            } // End cycle over number of parents

            if (GL_tentacles_size[lv + 1] == 0)
            {
                lv = lmax - lmin + 1;
            }
        }

        //** >> Free pointers **/
        free(pptr_red_black[0]); // Free red 
        pptr_red_black[0] = NULL;
        free(pptr_red_black[1]); // Free black
        pptr_red_black[1] = NULL;
        free(pptr_red_black); // Free red and black
        pptr_red_black = NULL;
        free(ptr_red_black_cap); // Free capacity of red and black
        ptr_red_black_cap = NULL;
        free(ptr_red_black_size); // Free size of red and black
        ptr_red_black_size = NULL;
        
    }

        return _SUCCESS_;
}