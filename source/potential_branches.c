/*
 * potential_branches.c
 *
 * Compute the potential in the branch nodes of the level of refinement
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

#include "potential_branches.h"

static int fill_red_and_black_branch_node(int **pptr_red_black, int *ptr_red_black_cap, int *ptr_red_black_size, const struct node *ptr_node)
{
	int red_size;	  // Number of red points
	int black_size;	  // Number of black points
	int red_or_black; // Desition if the point is red or black

	// int box_idx_x; // Box index in X direcction
	// int box_idx_y; // Box index in Y direcction
	// int box_idx_z; // Box index in Z direcction

	int box_grid_idx; // box grid index

	red_size = 0;
	black_size = 0;

	//** >> Filling Red and Black arrays with the correspondig grid points **/
	for (int i = 0; i < ptr_node->grid_intr_size; i++)
	{
		box_grid_idx = ptr_node->ptr_intr_grid_idx[i];
		// box_idx_z = box_grid_idx / ((ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1));
		// box_idx_y = (box_grid_idx - box_idx_z * (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1)) / (ptr_node->box_real_dim_x + 1);
		// box_idx_x = box_grid_idx - box_idx_z * (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1) - box_idx_y * (ptr_node->box_real_dim_x + 1);
		// box_idx_x = ptr_node->ptr_intr_grid_cell_idx_x[i];
		// box_idx_y = ptr_node->ptr_intr_grid_cell_idx_y[i];
		// box_idx_z = ptr_node->ptr_intr_grid_cell_idx_z[i];
		// red_or_black = box_idx_x + box_idx_y + box_idx_z;

		red_or_black = ptr_node->ptr_intr_grid_cell_idx_x[i] +
					   ptr_node->ptr_intr_grid_cell_idx_y[i] +
					   ptr_node->ptr_intr_grid_cell_idx_z[i];

		//** >> Space checking of red and black pointers **/
		if (red_or_black % 2 == 0) // If it is red
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

static void initial_potential_branch_node(const struct node *ptr_node, struct node *ptr_ch)
{
	int ch_box_grid_idx; // Child box grid index
	// int ch_box_idx_x;   // Child box index
	// int ch_box_idx_y;
	// int ch_box_idx_z;

	// int ch_cell_idx_x; // Chilld cell index at x direcction
	// int ch_cell_idx_y;
	// int ch_cell_idx_z;

	int box_idx_x0_node; // Parent box index at X direcction
	int box_idx_x1_node;
	int box_idx_y0_node;
	int box_idx_y1_node;
	int box_idx_z0_node;
	int box_idx_z1_node;

	int box_grid_i0_j0_k0_node; // Parent grid index at i0, j0, k0 position
	int box_grid_i1_j0_k0_node;
	int box_grid_i0_j1_k0_node;
	int box_grid_i1_j1_k0_node;
	int box_grid_i0_j0_k1_node;
	int box_grid_i1_j0_k1_node;
	int box_grid_i0_j1_k1_node;
	int box_grid_i1_j1_k1_node;

	vtype aux_pot; // Auxiliar Potential in the grid point

	int grid_box_real_dim_X_node = (ptr_node->box_real_dim_x + 1);
	int grid_box_real_dim_X_times_Y_node = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

	//** >> Passing the potential from coarse parent to fine child **/
	// Border grid points
	for (int i = 0; i < ptr_ch->grid_bder_size; i++)
	{
		//** >> Child box indexes **/
		ch_box_grid_idx = ptr_ch->ptr_bder_grid_idx[i];
		// ch_box_idx_z = ch_box_grid_idx / ((ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1));
		// ch_box_idx_y = (ch_box_grid_idx - ch_box_idx_z * (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1)) / (ptr_ch->box_real_dim_x + 1);
		// ch_box_idx_x = ch_box_grid_idx - ch_box_idx_z * (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1) - ch_box_idx_y * (ptr_ch->box_real_dim_x + 1);

		//** >> Child cell indexes **/
		// ch_cell_idx_x = ch_box_idx_x + ptr_ch->box_ts_x;
		// ch_cell_idx_y = ch_box_idx_y + ptr_ch->box_ts_y;
		// ch_cell_idx_z = ch_box_idx_z + ptr_ch->box_ts_z;

		// ch_cell_idx_x = ptr_ch->ptr_bder_grid_cell_idx_x[i];
		// ch_cell_idx_y = ptr_ch->ptr_bder_grid_cell_idx_y[i];
		// ch_cell_idx_z = ptr_ch->ptr_bder_grid_cell_idx_z[i];

		//** >> Parent box indexes **/
		box_idx_x0_node = (ptr_ch->ptr_bder_grid_cell_idx_x[i] >> 1) - ptr_node->box_ts_x;
		box_idx_x1_node = ((ptr_ch->ptr_bder_grid_cell_idx_x[i] + 1) >> 1) - ptr_node->box_ts_x;
		box_idx_y0_node = (ptr_ch->ptr_bder_grid_cell_idx_y[i] >> 1) - ptr_node->box_ts_y;
		box_idx_y1_node = ((ptr_ch->ptr_bder_grid_cell_idx_y[i] + 1) >> 1) - ptr_node->box_ts_y;
		box_idx_z0_node = (ptr_ch->ptr_bder_grid_cell_idx_z[i] >> 1) - ptr_node->box_ts_z;
		box_idx_z1_node = ((ptr_ch->ptr_bder_grid_cell_idx_z[i] + 1) >> 1) - ptr_node->box_ts_z;

		//** >> Parent grid indexes **/
		box_grid_i0_j0_k0_node = box_idx_x0_node + box_idx_y0_node * grid_box_real_dim_X_node + box_idx_z0_node * grid_box_real_dim_X_times_Y_node;
		box_grid_i1_j0_k0_node = box_idx_x1_node + box_idx_y0_node * grid_box_real_dim_X_node + box_idx_z0_node * grid_box_real_dim_X_times_Y_node;
		box_grid_i0_j1_k0_node = box_idx_x0_node + box_idx_y1_node * grid_box_real_dim_X_node + box_idx_z0_node * grid_box_real_dim_X_times_Y_node;
		box_grid_i1_j1_k0_node = box_idx_x1_node + box_idx_y1_node * grid_box_real_dim_X_node + box_idx_z0_node * grid_box_real_dim_X_times_Y_node;
		box_grid_i0_j0_k1_node = box_idx_x0_node + box_idx_y0_node * grid_box_real_dim_X_node + box_idx_z1_node * grid_box_real_dim_X_times_Y_node;
		box_grid_i1_j0_k1_node = box_idx_x1_node + box_idx_y0_node * grid_box_real_dim_X_node + box_idx_z1_node * grid_box_real_dim_X_times_Y_node;
		box_grid_i0_j1_k1_node = box_idx_x0_node + box_idx_y1_node * grid_box_real_dim_X_node + box_idx_z1_node * grid_box_real_dim_X_times_Y_node;
		box_grid_i1_j1_k1_node = box_idx_x1_node + box_idx_y1_node * grid_box_real_dim_X_node + box_idx_z1_node * grid_box_real_dim_X_times_Y_node;

		aux_pot = 0.125 * (ptr_node->ptr_pot[box_grid_i0_j0_k0_node] +
						   ptr_node->ptr_pot[box_grid_i1_j0_k0_node] +
						   ptr_node->ptr_pot[box_grid_i0_j1_k0_node] +
						   ptr_node->ptr_pot[box_grid_i1_j1_k0_node] +
						   ptr_node->ptr_pot[box_grid_i0_j0_k1_node] +
						   ptr_node->ptr_pot[box_grid_i1_j0_k1_node] +
						   ptr_node->ptr_pot[box_grid_i0_j1_k1_node] +
						   ptr_node->ptr_pot[box_grid_i1_j1_k1_node]);

		ptr_ch->ptr_pot[ch_box_grid_idx] = aux_pot;
	}

	// Interior grid points
	for (int i = 0; i < ptr_ch->grid_intr_size; i++)
	{
		//** >> Child box indexes **/
		ch_box_grid_idx = ptr_ch->ptr_intr_grid_idx[i];

		// ch_box_idx_z = ch_box_grid_idx / ((ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1));
		// ch_box_idx_y = (ch_box_grid_idx - ch_box_idx_z * (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1)) / (ptr_ch->box_real_dim_x + 1);
		// ch_box_idx_x = ch_box_grid_idx - ch_box_idx_z * (ptr_ch->box_real_dim_x + 1) * (ptr_ch->box_real_dim_y + 1) - ch_box_idx_y * (ptr_ch->box_real_dim_x + 1);

		//** >> Child cell indexes **/
		// ch_cell_idx_x = ch_box_idx_x + ptr_ch->box_ts_x;
		// ch_cell_idx_y = ch_box_idx_y + ptr_ch->box_ts_y;
		// ch_cell_idx_z = ch_box_idx_z + ptr_ch->box_ts_z;

		// ch_cell_idx_x = ptr_ch->ptr_intr_grid_cell_idx_x[i];
		// ch_cell_idx_y = ptr_ch->ptr_intr_grid_cell_idx_y[i];
		// ch_cell_idx_z = ptr_ch->ptr_intr_grid_cell_idx_z[i];

		//** >> Parent box indexes **/
		box_idx_x0_node = (ptr_ch->ptr_intr_grid_cell_idx_x[i] >> 1) - ptr_node->box_ts_x;
		box_idx_x1_node = ((ptr_ch->ptr_intr_grid_cell_idx_x[i] + 1) >> 1) - ptr_node->box_ts_x;
		box_idx_y0_node = (ptr_ch->ptr_intr_grid_cell_idx_y[i] >> 1) - ptr_node->box_ts_y;
		box_idx_y1_node = ((ptr_ch->ptr_intr_grid_cell_idx_y[i] + 1) >> 1) - ptr_node->box_ts_y;
		box_idx_z0_node = (ptr_ch->ptr_intr_grid_cell_idx_z[i] >> 1) - ptr_node->box_ts_z;
		box_idx_z1_node = ((ptr_ch->ptr_intr_grid_cell_idx_z[i] + 1) >> 1) - ptr_node->box_ts_z;

		//** >> Parent grid indexes **/
		box_grid_i0_j0_k0_node = box_idx_x0_node + box_idx_y0_node * grid_box_real_dim_X_node + box_idx_z0_node * grid_box_real_dim_X_times_Y_node;
		box_grid_i1_j0_k0_node = box_idx_x1_node + box_idx_y0_node * grid_box_real_dim_X_node + box_idx_z0_node * grid_box_real_dim_X_times_Y_node;
		box_grid_i0_j1_k0_node = box_idx_x0_node + box_idx_y1_node * grid_box_real_dim_X_node + box_idx_z0_node * grid_box_real_dim_X_times_Y_node;
		box_grid_i1_j1_k0_node = box_idx_x1_node + box_idx_y1_node * grid_box_real_dim_X_node + box_idx_z0_node * grid_box_real_dim_X_times_Y_node;
		box_grid_i0_j0_k1_node = box_idx_x0_node + box_idx_y0_node * grid_box_real_dim_X_node + box_idx_z1_node * grid_box_real_dim_X_times_Y_node;
		box_grid_i1_j0_k1_node = box_idx_x1_node + box_idx_y0_node * grid_box_real_dim_X_node + box_idx_z1_node * grid_box_real_dim_X_times_Y_node;
		box_grid_i0_j1_k1_node = box_idx_x0_node + box_idx_y1_node * grid_box_real_dim_X_node + box_idx_z1_node * grid_box_real_dim_X_times_Y_node;
		box_grid_i1_j1_k1_node = box_idx_x1_node + box_idx_y1_node * grid_box_real_dim_X_node + box_idx_z1_node * grid_box_real_dim_X_times_Y_node;

		aux_pot = 0.125 * (ptr_node->ptr_pot[box_grid_i0_j0_k0_node] +
						   ptr_node->ptr_pot[box_grid_i1_j0_k0_node] +
						   ptr_node->ptr_pot[box_grid_i0_j1_k0_node] +
						   ptr_node->ptr_pot[box_grid_i1_j1_k0_node] +
						   ptr_node->ptr_pot[box_grid_i0_j0_k1_node] +
						   ptr_node->ptr_pot[box_grid_i1_j0_k1_node] +
						   ptr_node->ptr_pot[box_grid_i0_j1_k1_node] +
						   ptr_node->ptr_pot[box_grid_i1_j1_k1_node]);

		ptr_ch->ptr_pot[ch_box_grid_idx] = aux_pot;
	}
}

static int compute_potential_branch_node(int **pptr_red_black, const int *ptr_red_black_size, struct node *ptr_node)
{
	int red_size;
	int black_size;
	vtype aux_pot;

	red_size = ptr_red_black_size[0];
	black_size = ptr_red_black_size[1];

	int box_grid_idx;
	int box_grid_idxNbr_x_plus;	 // Box grid index in the neigborhood on the right
	int box_grid_idxNbr_x_minus; // Box grid index in the neigborhood on the left
	int box_grid_idxNbr_y_plus;	 // Box grid index in the neigborhood behind
	int box_grid_idxNbr_y_minus; // Box grid index in the neigborhood in front
	int box_grid_idxNbr_z_plus;	 // Box grid index in the neigborhood up
	int box_grid_idxNbr_z_minus; // Box grid index in the neigborhood down

	int grid_box_real_dim_X_node = (ptr_node->box_real_dim_x + 1);
	int grid_box_real_dim_X_times_Y_node = (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1);

	vtype H = 1.0L / (1 << ptr_node->lv);
	vtype H_pow2 = H * H;

	//** >> Cycle over the Successive over-relaxation **/
	for (int iter = 0; iter < _Iter_branches_SOR_; iter++)
	{
		//** >> Cycle over red points **/
		for (int i = 0; i < red_size; i++)
		{
			box_grid_idx = pptr_red_black[0][i];
			box_grid_idxNbr_x_plus = box_grid_idx + 1;
			box_grid_idxNbr_x_minus = box_grid_idx - 1;
			box_grid_idxNbr_y_plus = box_grid_idx + grid_box_real_dim_X_node;
			box_grid_idxNbr_y_minus = box_grid_idx - grid_box_real_dim_X_node;
			box_grid_idxNbr_z_plus = box_grid_idx + grid_box_real_dim_X_times_Y_node;
			box_grid_idxNbr_z_minus = box_grid_idx - grid_box_real_dim_X_times_Y_node;

			aux_pot = 1.0L / 6.0L * (-H_pow2 * ptr_node->ptr_d[box_grid_idx] + ptr_node->ptr_pot[box_grid_idxNbr_x_plus] + ptr_node->ptr_pot[box_grid_idxNbr_x_minus] + ptr_node->ptr_pot[box_grid_idxNbr_y_plus] + ptr_node->ptr_pot[box_grid_idxNbr_y_minus] + ptr_node->ptr_pot[box_grid_idxNbr_z_plus] + ptr_node->ptr_pot[box_grid_idxNbr_z_minus]);
			ptr_node->ptr_pot[box_grid_idx] = _w_SOR_ * aux_pot + (1.0 - _w_SOR_) * ptr_node->ptr_pot[box_grid_idx];
		}

		//** >> Cycle over black points **/
		for (int i = 0; i < black_size; i++)
		{
			box_grid_idx = pptr_red_black[1][i];
			box_grid_idxNbr_x_plus = box_grid_idx + 1;
			box_grid_idxNbr_x_minus = box_grid_idx - 1;
			box_grid_idxNbr_y_plus = box_grid_idx + grid_box_real_dim_X_node;
			box_grid_idxNbr_y_minus = box_grid_idx - grid_box_real_dim_X_node;
			box_grid_idxNbr_z_plus = box_grid_idx + grid_box_real_dim_X_times_Y_node;
			box_grid_idxNbr_z_minus = box_grid_idx - grid_box_real_dim_X_times_Y_node;

			aux_pot = 1.0L / 6.0L * (-H_pow2 * ptr_node->ptr_d[box_grid_idx] + ptr_node->ptr_pot[box_grid_idxNbr_x_plus] + ptr_node->ptr_pot[box_grid_idxNbr_x_minus] + ptr_node->ptr_pot[box_grid_idxNbr_y_plus] + ptr_node->ptr_pot[box_grid_idxNbr_y_minus] + ptr_node->ptr_pot[box_grid_idxNbr_z_plus] + ptr_node->ptr_pot[box_grid_idxNbr_z_minus]);
			ptr_node->ptr_pot[box_grid_idx] = _w_SOR_ * aux_pot + (1.0 - _w_SOR_) * ptr_node->ptr_pot[box_grid_idx];
		}
	}

	return _SUCCESS_;
}

const int conjugate_gradient_branch_node(struct node *ptr_node)
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

	//** >> CHEKING ERROR SOLUTION CONDITION **/
	aux_clock = clock();
	if (poisson_error(ptr_node, r, rsold, 1) == true)
	{
		GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
		free(r);
		free(p);
		free(Ap);
		return _SUCCESS_;
	}
	GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

	int cycle_size = size < _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_ ? size : _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_;
	// For cycle over mutually conjugate vectors p
	for (int i = 0; i < cycle_size; i++)
	{

		// Computing Ap
		for (int j = 0; j < ptr_node->grid_intr_size; j++)
		{
			box_grid_idx = ptr_node->ptr_intr_grid_idx[j];
			Ap[box_grid_idx] = -one_over_H_pow_2 * (6.0 * p[box_grid_idx] - p[box_grid_idx + 1] - p[box_grid_idx - 1] - p[box_grid_idx + grid_box_real_dim_X] - p[box_grid_idx - grid_box_real_dim_X] - p[box_grid_idx + grid_box_real_dim_X_times_Y] - p[box_grid_idx - grid_box_real_dim_X_times_Y]);
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
		aux_clock = clock();
		if (poisson_error(ptr_node, r, rsnew, 1) == true)
		{
			GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
			free(r);
			free(p);
			free(Ap);
			return _SUCCESS_;
		}
		GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

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

int potential_branches()
{
	struct node *ptr_node;
	struct node *ptr_ch;

	int no_pts; // Number of parents in the cycle

	int **pptr_red_black;
	int *ptr_red_black_cap;
	int *ptr_red_black_size;

	int iter;
	bool check;

	clock_t aux_clock;

	vtype *dummy_pvtype = NULL;
	vtype dummy_vtype = 0;

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

	for (int lv = 0; lv < GL_tentacles_level_max; lv++)
	{
		no_pts = GL_tentacles_size[lv];

		//** >> For cycle over parent nodes **/
		for (int i = 0; i < no_pts; i++)
		{
			ptr_node = GL_tentacles[lv][i];
			for (int j = 0; j < ptr_node->chn_size; j++)
			{
				ptr_ch = ptr_node->pptr_chn[j];
				//** >> Transfer the potential from parent to child nodes **/
				aux_clock = clock();
				initial_potential_branch_node(ptr_node, ptr_ch);
				GL_times[22] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

				if (ptr_ch->cell_size < 513)	//216 = node with minimum size of 1 (+1 n_exp) size, (1 + 2*n_exp)^3 * 8 
				{
					//** >> Computing the potential in the child node **/
					aux_clock = clock();
					if (conjugate_gradient_branch_node(ptr_ch) == _FAILURE_)
					{
						printf("Error at function conjugate_gradient_branch()\n");
						return _FAILURE_;
					}
					GL_times[24] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
				}
				else
				{
					//** >> Computing the potential in the child node **/
					aux_clock = clock();
					if (fill_red_and_black_branch_node(pptr_red_black, ptr_red_black_cap, ptr_red_black_size, ptr_ch) == _FAILURE_)
					{
						printf("Error at function fill_red_and_black()\n");
						return _FAILURE_;
					}
					GL_times[23] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

					//** >> CHEKING ERROR SOLUTION CONDITION **/
					aux_clock = clock();
					check = poisson_error(ptr_ch,dummy_pvtype,dummy_vtype,0);
					GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

					iter = 0;
					while (iter < _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_ && check == false)
					{
						iter = iter + 1;

						//** >> Computing the potential in the child node **/
						aux_clock = clock();
						if (compute_potential_branch_node(pptr_red_black, ptr_red_black_size, ptr_ch) == _FAILURE_)
						{
							printf("Error at function potential_branch_node()\n");
							return _FAILURE_;
						}
						GL_times[24] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;

						//** >> CHEKING ERROR SOLUTION CONDITION **/
						aux_clock = clock();
						check = poisson_error(ptr_ch,dummy_pvtype,dummy_vtype,0);
						GL_times[25] += (double)(clock() - aux_clock) / CLOCKS_PER_SEC;
					}
					if (iter == _MAX_NUMBER_OF_ITERATIONS_IN_POISSON_EQUATION_)
					{
						printf("\nERROR: The precision was not reached in the branch node. Too many SOR Iterations, plz choose a lower precision than %1.3e\n", (double)_ERROR_THRESHOLD_IN_THE_POISSON_EQUATION_);
						return _FAILURE_;
					}
				}

			} // End cycle over number of children
		}	  // End cycle over number of parents
	}

	//** >> Free pointers **/
	if (pptr_red_black[0] != NULL)
	{
		free(pptr_red_black[0]); // Free red
	}
	if (pptr_red_black[1] != NULL)
	{
		free(pptr_red_black[1]); // Free black
	}
	free(pptr_red_black);	  // Free red and black
	free(ptr_red_black_cap);  // Free capacity of red and black
	free(ptr_red_black_size); // Free size of red and black

	return _SUCCESS_;
}
