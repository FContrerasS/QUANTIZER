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

	//** >> Potential, Acceleration and density of the grid **/
	int cap = (ptr_head->box_real_dim_x + 1) * (ptr_head->box_real_dim_y + 1) * (ptr_head->box_real_dim_z + 1);

	for (int i = 0; i < cap; i++)
	{
		ptr_head->ptr_d[i] = 0;
	}

	//** >> Initial values for the potential and acceleration **/
	//** Here we are modified the boundary values of the potential and acceleration
	initial_potential_and_acceleration_head(ptr_head);


	ptr_node = NULL;
	ptr_head = NULL;

	return _SUCCESS_;
}