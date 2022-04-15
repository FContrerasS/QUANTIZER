/*
 * initialization.c
 *
 * Define basic parameters and structures of the simulation
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

#include "initialization.h"

static void initializing_particle_flag_updating()
{
	GL_ptcl_updating_flag = (bool *)malloc(GL_no_ptcl * sizeof(bool));

	for (int i = 0;i< GL_no_ptcl; i++)
	{
		GL_ptcl_updating_flag[i] = false;
	}
}

static void initializing_head_node()
{
	struct node *ptr_head;
	ptr_head = GL_ptr_tree;

	//** >> Basic values of the node **/
	initialize_node(ptr_head);

	//** >> Global node properties **/
	ptr_head->ID = 0;
	ptr_head->lv = lmin;

	//** >> Cells in the node **/
	ptr_head->ptr_cell_idx_x = (int *)malloc(no_lmin_cell_pow3 * sizeof(int));
	ptr_head->ptr_cell_idx_y = (int *)malloc(no_lmin_cell_pow3 * sizeof(int));
	ptr_head->ptr_cell_idx_z = (int *)malloc(no_lmin_cell_pow3 * sizeof(int));
	int idx;
	for (int k = 0; k < no_lmin_cell; k++)
	{
		for (int j = 0; j < no_lmin_cell; j++)
		{
			for (int i = 0; i < no_lmin_cell; i++)
			{
				idx = i + j * no_lmin_cell + k * no_lmin_cell_pow2;
				ptr_head->ptr_cell_idx_x[idx] = i;
				ptr_head->ptr_cell_idx_y[idx] = j;
				ptr_head->ptr_cell_idx_z[idx] = k;
			}
		}
	}
	ptr_head->cell_cap = no_lmin_cell_pow3;
	ptr_head->cell_size = no_lmin_cell_pow3;

	//** >> Particles in the node **/
	ptr_head->ptcl_cap = GL_no_ptcl;
	ptr_head->ptcl_size = GL_no_ptcl;
	ptr_head->ptr_ptcl = (int *)malloc(ptr_head->ptcl_cap * sizeof(int));
	for (int i = 0; i < ptr_head->ptcl_size; i++)
	{
		ptr_head->ptr_ptcl[i] = i;
	}

	//** >> Boxes **/
	// Including border when particles goes out of the simulation and also because the n_exp parameter checking
	ptr_head->ptr_box = (int *)malloc(box_side_lmin_pow3 * sizeof(int));
	ptr_head->ptr_box_aux = (int *)malloc(box_side_lmin_pow3 * sizeof(int));
	int box_idx; // Index in the box
	for (int k = 0; k < box_side_lmin; k++)
	{
		for (int j = 0; j < box_side_lmin; j++)
		{
			for (int i = 0; i < box_side_lmin; i++)
			{
				box_idx = i + j * box_side_lmin + k * box_side_lmin_pow2;
				if (i <= bder_os_sim - 1 || i >= box_side_lmin - bder_os_sim || j <= bder_os_sim - 1 || j >= box_side_lmin - bder_os_sim || k <= bder_os_sim - 1 || k >= box_side_lmin - bder_os_sim)
				{
					ptr_head->ptr_box[box_idx] = -5; //-5 corresponds to cell out of the box simulation
				}
				else if (i == bder_os_sim || i == box_side_lmin - bder_os_sim - 1 || j == bder_os_sim || j == box_side_lmin - bder_os_sim - 1 || k == bder_os_sim || k == box_side_lmin - bder_os_sim - 1)
				{
					ptr_head->ptr_box[box_idx] = -2; //-2 corresponds to cell in the border of the block
				}
				else
				{
					ptr_head->ptr_box[box_idx] = -3; //-3 corresponds to cell belonging to the block
				}
			}
		}
	}

	ptr_head->box_dim_x = no_lmin_cell;
	ptr_head->box_dim_y = no_lmin_cell;
	ptr_head->box_dim_z = no_lmin_cell;
	ptr_head->box_real_dim_x = box_side_lmin;
	ptr_head->box_real_dim_y = box_side_lmin;
	ptr_head->box_real_dim_z = box_side_lmin;
	ptr_head->box_ts_x = -bder_os_sim;
	ptr_head->box_ts_y = -bder_os_sim;
	ptr_head->box_ts_z = -bder_os_sim;

	//** >> Grid points **/
	ptr_head->grid_intr_cap = (no_lmin_cell - 1) * (no_lmin_cell - 1) * (no_lmin_cell - 1);
	ptr_head->grid_intr_size = 0;
	ptr_head->grid_bder_cap = (no_lmin_cell + 1) * (no_lmin_cell + 1) * (no_lmin_cell + 1) - ptr_head->grid_intr_cap;
	ptr_head->grid_bder_size = 0;
	ptr_head->ptr_grid_intr = (int *)malloc(ptr_head->grid_intr_cap * sizeof(int));
	ptr_head->ptr_grid_bder = (int *)malloc(ptr_head->grid_bder_cap * sizeof(int));
	int box_grid_idx; // Index in the box grid point
	//** Filling grid points indexes **/
	for (int k = bder_os_sim; k < box_side_lmin + 1 - bder_os_sim; k++)
	{
		for (int j = bder_os_sim; j < box_side_lmin + 1 - bder_os_sim; j++)
		{
			for (int i = bder_os_sim; i < box_side_lmin + 1 - bder_os_sim; i++)
			{
				box_grid_idx = i + j * (box_side_lmin + 1) + k * (box_side_lmin + 1) * (box_side_lmin + 1);
				//** >> Border grid point **/
				if (i == bder_os_sim || i == box_side_lmin - bder_os_sim || j == bder_os_sim || j == box_side_lmin - bder_os_sim || k == bder_os_sim || k == box_side_lmin - bder_os_sim)
				{
					ptr_head->ptr_grid_bder[ptr_head->grid_bder_size] = box_grid_idx;
					ptr_head->grid_bder_size += 1; // Increasing the border grid points
				}
				//** >> Interior grid point **/
				else
				{
					ptr_head->ptr_grid_intr[ptr_head->grid_intr_size] = box_grid_idx;
					ptr_head->grid_intr_size += 1; // Increasing the interior grid points
				}
			}
		}
	}

	//** >> Refinement Criterion **/
	// The default mass is 0
	if (lmin < lmax)
	{
		ptr_head->ptr_box_mass = (vtype *)calloc(box_side_lmin_pow3, sizeof(vtype));
	}
	ptr_head->local_mass = total_mass;

	//** >> Potential, Acceleration and density of the grid **/
	int cap = (ptr_head->box_real_dim_x + 1) * (ptr_head->box_real_dim_y + 1) * (ptr_head->box_real_dim_z + 1);
	ptr_head->ptr_pot = (vtype *)calloc(cap, sizeof(vtype)); // Potential
	ptr_head->ptr_ax = (vtype *)calloc(cap, sizeof(vtype));	 // Acceleration
	ptr_head->ptr_ay = (vtype *)calloc(cap, sizeof(vtype));
	ptr_head->ptr_az = (vtype *)calloc(cap, sizeof(vtype));
	ptr_head->ptr_d = (vtype *)calloc(cap, sizeof(vtype));	 // Density

	//** >> Initial values for the potential and acceleration **/
	initial_potential_and_acceleration_head(ptr_head);

} // end function initializing_main_node

static void initializing_tentacles()
{
	GL_tentacles_old = (struct node ***)malloc((lmax - lmin + 1) * sizeof(struct node **));
	GL_tentacles_new = (struct node ***)malloc((lmax - lmin + 1) * sizeof(struct node **));
	GL_tentacles_cap = (int *)calloc((lmax - lmin + 1), sizeof(int)); // Capacity of pointer of the tentacles in each level
	GL_tentacles_size = (int *)calloc((lmax - lmin + 1), sizeof(int));

	//** >> Head node **/
	GL_tentacles_old[0] = (struct node **)malloc(1 * sizeof(struct node *));	
	GL_tentacles_new[0] = (struct node **)malloc(1 * sizeof(struct node *));
	GL_tentacles_cap[0] = 1;
	GL_tentacles_size[0] = 1;

	for (int i = 1; i < (lmax - lmin + 1); i++)
	{
		GL_tentacles_old[i] = NULL;
		GL_tentacles_new[i] = NULL;
	}

	//** >> Filling Head node **/
	GL_tentacles_old[0][0] = GL_ptr_tree;
	GL_tentacles_new[0][0] = GL_ptr_tree;

	GL_tentacles_level_max = 0; // Maximum depth of the tentacles

	//** >> Memory computing **/
	TOTAL_MEMORY_TENTACLES += 2 * (lmax - lmin + 1) * ( sizeof(struct node **) + sizeof(int)) + 2 * sizeof(struct node *);
}

int initialization()
{

	//** >> Particle updating flag initialization **/
	initializing_particle_flag_updating();

	// struct node *GL_ptr_tree; // Head or Main node
	//** >> Initializing the Head of the Tree **/
	initializing_head_node();

	//** >>  Initializing tentacles struct pointer **/
	initializing_tentacles();




	return _SUCCESS_;
}
