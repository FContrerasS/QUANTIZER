/*
 * reset.c
 *
 * Reset parameters before the new potential computation
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

#include "reset.h"

// static void cleaning()
// {
// 	//** >> Cleaning nodes**/

// 	struct node *ptr_node;

// 	int no_pts; // Number of parents in the cycle

// 	for (int lv = GL_tentacles_level_max; lv > 0; lv--)
// 	{
// 		no_pts = GL_tentacles_size[lv];
// 		//** >> For cycle over parent nodes **/
// 		for (int i = 0; i < no_pts; i++)
// 		{
// 			ptr_node = GL_tentacles[lv][i];

// 			free(ptr_node->ptr_cell_idx_x);
// 			free(ptr_node->ptr_cell_idx_y);
// 			free(ptr_node->ptr_cell_idx_z);
// 			free(ptr_node->ptr_ptcl);
// 			free(ptr_node->ptr_box);
// 			free(ptr_node->ptr_box_aux);
// 			free(ptr_node->ptr_grid_intr);
// 			free(ptr_node->ptr_grid_bder);
// 			free(ptr_node->ptr_box_mass);
// 			free(ptr_node->ptr_pot);
// 			free(ptr_node->ptr_ax);
// 			free(ptr_node->ptr_ay);
// 			free(ptr_node->ptr_az);
// 			free(ptr_node->ptr_d);
// 			ptr_node->ptr_pt = NULL;
// 			for (int j = 0; j < ptr_node->zones_size; j++)
// 			{
// 				ptr_node->pptr_chn[j] = NULL;
// 				free(ptr_node->pptr_zones[j]);
// 			}
// 			free(ptr_node->pptr_chn);
// 			free(ptr_node->ptr_cell_ref);
// 			free(ptr_node->pptr_zones);
// 			free(ptr_node->ptr_zone_cap);
// 			free(ptr_node->ptr_zone_size);
// 			free(ptr_node->ptr_aux_idx);
// 			free(ptr_node);
// 		}
// 	}
// }

int reset()
{

	TOTAL_MEMORY_NODES = 0;
	TOTAL_MEMORY_CELDAS = 0;
	TOTAL_MEMORY_PARTICULAS = 0;
	TOTAL_MEMORY_CAJAS = 0;
	TOTAL_MEMORY_OTROS = 0;
	TOTAL_MEMORY_TENTACLES = 0;
	TOTAL_MEMORY_AUX = 0;

	//Global Particles properties
	TOTAL_MEMORY_OTROS += GL_no_ptcl * 10 * sizeof(vtype);	

	//Head and brances nodes
	struct node *ptr_node = NULL;
	int no_pts;
	for (int lv = GL_tentacles_level_max; lv > -1; lv--)
	{
		no_pts = GL_tentacles_size[lv];
		//** >> For cycle over parent nodes **/
		for (int i = 0; i < no_pts; i++)
		{
			ptr_node = GL_tentacles[lv][i];

			TOTAL_MEMORY_NODES += sizeof(struct node);
			TOTAL_MEMORY_CELDAS += 3 * ptr_node->cell_cap * sizeof(int);
			TOTAL_MEMORY_PARTICULAS += ptr_node->ptcl_cap * sizeof(int);
			TOTAL_MEMORY_CAJAS += ptr_node->box_real_dim_x * ptr_node->box_real_dim_y * ptr_node->box_real_dim_z * (2 * sizeof(int) + sizeof(vtype));	  // old, new and mass density
			TOTAL_MEMORY_CAJAS += (ptr_node->box_real_dim_x + 1) * (ptr_node->box_real_dim_y + 1) * (ptr_node->box_real_dim_z + 1) * (5 * sizeof(vtype)); // Pot, 3 acc and dens
			TOTAL_MEMORY_AUX += ptr_node->cell_ref_cap * sizeof(int);
			TOTAL_MEMORY_AUX += ptr_node->zones_cap * sizeof(int *);
			for (int j = 0; j < ptr_node->zones_size;j++)
			{
				TOTAL_MEMORY_AUX += ptr_node->ptr_zone_size[j] * sizeof(int);
			}
			TOTAL_MEMORY_AUX += 2 * ptr_node->zones_cap * sizeof(int);
			TOTAL_MEMORY_AUX += 2 * ptr_node->aux_idx_cap * sizeof(int);
		}
	}



	//** >> Global particle accelerations **/
	for (int i = 0; i < GL_no_ptcl; i++)
	{
		GL_ptcl_ax[i] = 0;
		GL_ptcl_ay[i] = 0;
		GL_ptcl_az[i] = 0;
	}

	//** >> Head node **/
	struct node *ptr_head = NULL;
	ptr_head = GL_ptr_tree;

	// int box_idx_x; // Box index in X direcction
	// int box_idx_y; // Box index in Y direcction
	// int box_idx_z; // Box index in Z direcction
	// int box_idx;   // Box index
	// int cell_idx;  // The cell index is simply i of the for loop

	// if(lmin < lmax)
	// {
	// 	//** >> Reset of the boxes **/
	// 	for (int i = 0; i < ptr_head->cell_ref_size; i++)
	// 	{
	// 		cell_idx = ptr_head->ptr_cell_ref[i];
	// 		box_idx_x = ptr_head->ptr_cell_idx_x[cell_idx] - ptr_head->box_ts_x;
	// 		box_idx_y = ptr_head->ptr_cell_idx_y[cell_idx] - ptr_head->box_ts_y;
	// 		box_idx_z = ptr_head->ptr_cell_idx_z[cell_idx] - ptr_head->box_ts_z;
	// 		box_idx = box_idx_x + box_idx_y * ptr_head->box_real_dim_x + box_idx_z * ptr_head->box_real_dim_x * ptr_head->box_real_dim_y;
	// 		//** >> Restarting the box values
	// 		ptr_head->ptr_box[box_idx] = -3;
	// 		ptr_head->ptr_box_aux[box_idx] = -3;
	// 	}

	// 	// //** Reset of the mass box **/
	// 	for (int i = 0; i < box_side_lmin_pow3; i++)
	// 	{
	// 		ptr_head->ptr_box_mass[i] = 0;
	// 	}

	// 	// //** >> Tree structure reset **/
	// 	ptr_head->chn_size = 0;

	// 	// //** >> Reset of the refinement zones auxiliary arrays **/
	// 	ptr_head->cell_ref_size = 0;
	// 	for (int i = 0; i < ptr_head->zones_size;i++)
	// 	{
	// 		ptr_head->ptr_zone_size[i] = 0;
	// 	}
	// 	ptr_head->zones_size = 0;
	// }

	//** >> Potential, Acceleration and density of the grid **/
	int cap = (ptr_head->box_real_dim_x + 1) * (ptr_head->box_real_dim_y + 1) * (ptr_head->box_real_dim_z + 1);

	for (int i = 0; i < cap; i++)
	{
		ptr_head->ptr_d[i] = 0;
	}

	//** >> Initial values for the potential and acceleration **/
	//** Here we are modified the boundary values of the potential and acceleration
	initial_potential_and_acceleration_head(ptr_head);

	//** Free branch and leaf nodes
	//cleaning();

	//** Reset Tentacles **/
	// for (int i = 1; i < GL_tentacles_level_max + 1; i++)
	// {
	// 	GL_tentacles_size[i] = 0;
	// }
	// GL_tentacles_level_max = 0;

	ptr_node = NULL;
	ptr_head = NULL;

	return _SUCCESS_;
}