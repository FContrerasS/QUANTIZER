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

//** >> Local Functions
static void initializing_global_parameters(void);
static void initializing_particle_updating_flag_and_ID(void);
static int initializing_head_node(void);
static void initializing_tentacles(void);

static void initializing_global_parameters(void)
{

	GL_total_mass_initial = 0;

	for(int i = 0; i< GL_no_ptcl_initial; i++)
	{
		GL_total_mass_initial += GL_ptcl_mass[i];
	}

	GL_total_mass_final = GL_total_mass_initial;

	meanmass = GL_total_mass_initial / GL_no_ptcl_initial;

	
}

static void initializing_particle_updating_flag_and_ID(void)
{
	//** >> Updating flag **/
	GL_ptcl_updating_flag = (bool *)malloc(GL_no_ptcl_initial * sizeof(bool));

	for (int i = 0;i< GL_no_ptcl_initial; i++)
	{
		GL_ptcl_updating_flag[i] = false;
	}

	//** >> Particle Id **/
	GL_ptcl_ID = (int *)malloc(GL_no_ptcl_initial * sizeof(int));

	for (int i = 0; i < GL_no_ptcl_initial; i++)
	{
		GL_ptcl_ID[i] = i;
	}
}

static int initializing_head_node(void)
{

	struct node *ptr_head;
	ptr_head = GL_ptr_tree;

	int box_idx_ptcl;

	int box_idx_x;
	int box_idx_y;
	int box_idx_z;

	int box_grid_idx; // Index in the box grid point

	//** >> Basic values of the node **/
	initialize_node(ptr_head);

	//** >> Global node properties **/
	ptr_head->ID = 0;
	ptr_head->lv = lmin;

	//** >> Boxes **/
	// Including border when particles goes out of the simulation and also because the n_exp parameter checking
	ptr_head->ptr_box = (int *)malloc(box_side_lmin_pow3 * sizeof(int));
	ptr_head->ptr_box_old = (int *)malloc(box_side_lmin_pow3 * sizeof(int));
	int box_idx; // Index in the box
	//-6 = teleport cell (Only for boundary_type 0) ; -5 = cell outside of the box simulation (Only for boundary_type 1 or 2); -4 = cell out of the node, the cell is in the parent node
	for (int k = 0; k < box_side_lmin; k++)
	{
		for (int j = 0; j < box_side_lmin; j++)
		{
			for (int i = 0; i < box_side_lmin; i++)
			{
				box_idx = i + j * box_side_lmin + k * box_side_lmin_pow2;
				if (i <= bder_os_sim - 1 || i >= box_side_lmin - bder_os_sim || j <= bder_os_sim - 1 || j >= box_side_lmin - bder_os_sim || k <= bder_os_sim - 1 || k >= box_side_lmin - bder_os_sim)
				{
					ptr_head->ptr_box[box_idx] = boundary_type == 0 ? -6 : -5; 
				}
				else
				{
					ptr_head->ptr_box[box_idx] = -3; //-3 corresponds to cell belonging to the block
				}
			}
		}
	}

	ptr_head->box_min_x = 0;
	ptr_head->box_min_y = 0;
	ptr_head->box_min_z = 0;
	ptr_head->box_max_x = (1 << lmin) - 1;
	ptr_head->box_max_y = (1 << lmin) - 1;
	ptr_head->box_max_z = (1 << lmin) - 1;

	ptr_head->box_dim_x = no_lmin_cell;
	ptr_head->box_dim_y = no_lmin_cell;
	ptr_head->box_dim_z = no_lmin_cell;
	ptr_head->box_real_dim_x = box_side_lmin;
	ptr_head->box_real_dim_y = box_side_lmin;
	ptr_head->box_real_dim_z = box_side_lmin;
	ptr_head->box_cap = box_side_lmin_pow3;
	ptr_head->box_ts_x = -bder_os_sim;
	ptr_head->box_ts_y = -bder_os_sim;
	ptr_head->box_ts_z = -bder_os_sim;

	//** >> Cells in the node **/
	ptr_head->ptr_cell_idx_x = (int *)malloc(no_lmin_cell_pow3 * sizeof(int));
	ptr_head->ptr_cell_idx_y = (int *)malloc(no_lmin_cell_pow3 * sizeof(int));
	ptr_head->ptr_cell_idx_z = (int *)malloc(no_lmin_cell_pow3 * sizeof(int));
	ptr_head->ptr_box_idx = (int *)malloc(no_lmin_cell_pow3 * sizeof(int));

	

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
				box_idx_x = i - ptr_head->box_ts_x;
				box_idx_y = j - ptr_head->box_ts_y;
				box_idx_z = k - ptr_head->box_ts_z;
				ptr_head->ptr_box_idx[idx] = box_idx_x + box_idx_y * ptr_head->box_real_dim_x + box_idx_z * ptr_head->box_real_dim_x * ptr_head->box_real_dim_y;
			}
		}
	}
	ptr_head->cell_cap = no_lmin_cell_pow3;
	ptr_head->cell_size = no_lmin_cell_pow3;

	//** >> Struct of cells (Particles and cell mass)
	if(lmin < lmax)
	{
		ptr_head->ptr_cell_struct = (struct cell_struct *)malloc(ptr_head->box_cap * sizeof(struct cell_struct));
		for (int j = 0; j < ptr_head->box_cap; j++)
		{
			initialize_cell_struct(&(ptr_head->ptr_cell_struct[j]));
		}

		for (int j = 0; j < GL_no_ptcl_initial; j++)
		{
			box_idx_ptcl = ptcl_idx_to_box_idx(ptr_head, j);

			//** >> Space checking of number of particles in the cell **/

			if (space_check(&(ptr_head->ptr_cell_struct[box_idx_ptcl].ptcl_cap), ptr_head->ptr_cell_struct[box_idx_ptcl].ptcl_size + 1, 1.0f, "p1i1", &(ptr_head->ptr_cell_struct[box_idx_ptcl].ptr_ptcl)) == _FAILURE_)
			{
				printf("Error, in space_check function\n");
				return _FAILURE_;
			}

			ptr_head->ptr_cell_struct[box_idx_ptcl].ptr_ptcl[ptr_head->ptr_cell_struct[box_idx_ptcl].ptcl_size] = j;

			ptr_head->ptr_cell_struct[box_idx_ptcl].ptcl_size += 1;

			ptr_head->ptr_cell_struct[box_idx_ptcl].cell_mass += GL_ptcl_mass[j];
		}
	}

	//** >> Total mass in the node
	ptr_head->local_mass = GL_total_mass_initial;
	//** >> Total number of particles in the node
	ptr_head->local_no_ptcl_full_node = GL_no_ptcl_initial;
	ptr_head->local_no_ptcl_to_use_outside_refinement_zones = GL_no_ptcl_initial;

	//** >> Grid points **/
	ptr_head->grid_intr_cap = (no_lmin_cell - 1) * (no_lmin_cell - 1) * (no_lmin_cell - 1);
	ptr_head->grid_intr_size = 0;
	ptr_head->ptr_intr_grid_cell_idx_x = (int *)malloc(ptr_head->grid_intr_cap * sizeof(int));
	ptr_head->ptr_intr_grid_cell_idx_y = (int *)malloc(ptr_head->grid_intr_cap * sizeof(int));
	ptr_head->ptr_intr_grid_cell_idx_z = (int *)malloc(ptr_head->grid_intr_cap * sizeof(int));
	ptr_head->ptr_intr_box_grid_idx = (int *)malloc(ptr_head->grid_intr_cap * sizeof(int));

	ptr_head->grid_SIMULATION_BOUNDARY_cap = (no_lmin_cell + 1) * (no_lmin_cell + 1) * (no_lmin_cell + 1) - ptr_head->grid_intr_cap;
	ptr_head->grid_SIMULATION_BOUNDARY_size = 0;
	ptr_head->ptr_SIMULATION_BOUNDARY_grid_cell_idx_x = (int *)malloc(ptr_head->grid_SIMULATION_BOUNDARY_cap * sizeof(int));
	ptr_head->ptr_SIMULATION_BOUNDARY_grid_cell_idx_y = (int *)malloc(ptr_head->grid_SIMULATION_BOUNDARY_cap * sizeof(int));
	ptr_head->ptr_SIMULATION_BOUNDARY_grid_cell_idx_z = (int *)malloc(ptr_head->grid_SIMULATION_BOUNDARY_cap * sizeof(int));
	ptr_head->ptr_SIMULATION_BOUNDARY_box_grid_idx = (int *)malloc(ptr_head->grid_SIMULATION_BOUNDARY_cap * sizeof(int));

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
					ptr_head->ptr_SIMULATION_BOUNDARY_grid_cell_idx_x[ptr_head->grid_SIMULATION_BOUNDARY_size] = i + ptr_head->box_ts_x;
					ptr_head->ptr_SIMULATION_BOUNDARY_grid_cell_idx_y[ptr_head->grid_SIMULATION_BOUNDARY_size] = j + ptr_head->box_ts_y;
					ptr_head->ptr_SIMULATION_BOUNDARY_grid_cell_idx_z[ptr_head->grid_SIMULATION_BOUNDARY_size] = k + ptr_head->box_ts_z;
					ptr_head->ptr_SIMULATION_BOUNDARY_box_grid_idx[ptr_head->grid_SIMULATION_BOUNDARY_size] = box_grid_idx;
					ptr_head->grid_SIMULATION_BOUNDARY_size += 1; // Increasing the border grid points
				}
				//** >> Interior grid point **/
				else
				{
					ptr_head->ptr_intr_grid_cell_idx_x[ptr_head->grid_intr_size] = i + ptr_head->box_ts_x;
					ptr_head->ptr_intr_grid_cell_idx_y[ptr_head->grid_intr_size] = j + ptr_head->box_ts_y;
					ptr_head->ptr_intr_grid_cell_idx_z[ptr_head->grid_intr_size] = k + ptr_head->box_ts_z;
					ptr_head->ptr_intr_box_grid_idx[ptr_head->grid_intr_size] = box_grid_idx;
					ptr_head->grid_intr_size += 1; // Increasing the interior grid points
				}
			}
		}
	}

	//** >> Potential, Acceleration and density of the grid **/
	int cap = (ptr_head->box_real_dim_x + 1) * (ptr_head->box_real_dim_y + 1) * (ptr_head->box_real_dim_z + 1);
	ptr_head->grid_properties_cap = cap;
	ptr_head->ptr_pot = (vtype *)calloc(cap, sizeof(vtype)); // Potential
	ptr_head->ptr_pot_old = (vtype *)calloc(cap, sizeof(vtype)); // Potential
	ptr_head->ptr_ax = (vtype *)calloc(cap, sizeof(vtype));	 // Acceleration
	ptr_head->ptr_ay = (vtype *)calloc(cap, sizeof(vtype));
	ptr_head->ptr_az = (vtype *)calloc(cap, sizeof(vtype));
	ptr_head->ptr_d = (vtype *)calloc(cap, sizeof(vtype));	 // Density

	//** >> Initial values for the potential and acceleration **/
	initial_potential_and_acceleration_head(ptr_head);

	//** >> Boundary of the simulation box **/
	ptr_head->boundary_simulation_contact = true;
	ptr_head->boundary_simulation_contact_x = true;
	ptr_head->boundary_simulation_contact_y = true;
	ptr_head->boundary_simulation_contact_z = true;

	if(boundary_type == 0)
	{
		ptr_head->pbc_crosses_the_boundary_simulation_box = true;
		ptr_head->pbc_crosses_the_boundary_simulation_box_x = true;
		ptr_head->pbc_crosses_the_boundary_simulation_box_y = true;
		ptr_head->pbc_crosses_the_boundary_simulation_box_z = true;

		ptr_head->pbc_crosses_the_whole_simulation_box = true;
		ptr_head->pbc_crosses_the_whole_simulation_box_x = true;
		ptr_head->pbc_crosses_the_whole_simulation_box_y = true;
		ptr_head->pbc_crosses_the_whole_simulation_box_z = true;
	}


	return _SUCCESS_;

} // end function initializing_main_node

static void initializing_tentacles(void)
{
	GL_tentacles = (struct node ***)malloc((lmax - lmin + 1) * sizeof(struct node **));
	GL_tentacles_cap = (int *)calloc((lmax - lmin + 1), sizeof(int)); // Capacity of pointer of the tentacles in each level
	GL_tentacles_size = (int *)calloc((lmax - lmin + 1), sizeof(int));

	for (int i = 0; i < (lmax - lmin + 1); i++)
	{
		GL_tentacles[i] = NULL;
	}

	//** >> Head node **/
	GL_tentacles[0] = (struct node **)malloc(1 * sizeof(struct node *));	
	GL_tentacles_cap[0] = 1;
	GL_tentacles_size[0] = 1;


	//** >> Filling Head node **/
	GL_tentacles[0][0] = GL_ptr_tree;

	GL_tentacles_level_max = 0; // Maximum depth of the tentacles
}

int initialization(void)
{

	//** >> Initializing global parameters **/
	initializing_global_parameters();

	//** >> Particle updating flag and ID initialization **/
	initializing_particle_updating_flag_and_ID();

	// struct node *GL_ptr_tree; // Head or Main node
	//** >> Initializing the Head of the Tree **/
	if (initializing_head_node() == _FAILURE_)
	{
		printf("Error at function initializing_head_node()\n");
		return _FAILURE_;
	}

	//** >> Initializing tentacles struct pointer **/
	initializing_tentacles();

	return _SUCCESS_;
}
