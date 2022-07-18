/*
 * garbage_collector.c
 *
 * free memory not used
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

#include "garbage_collector.h"

//** >> Local Functions
static int computing_memory(void);
static void free_memory_pool(void);
static int free_nodes_voids(void);

static int computing_memory(void)
{

	int no_pts;
	struct node *ptr_node;

	//** >> Tentacles **/
	TOTAL_MEMORY_TENTACLES += (lmax - lmin + 1) * sizeof(struct node **);

	for (int lv = 0; lv < lmax - lmin + 1; lv++)
	{
		TOTAL_MEMORY_TENTACLES += GL_tentacles_cap[lv] * sizeof(struct node *);
	}

	//** >> Global particles **/
	TOTAL_MEMORY_PARTICLES += GL_no_ptcl_initial * (sizeof(bool) + sizeof(int) + 10 * sizeof(vtype)); // Global particles

	//** >> Node properties **/
	for (int lv = 0; lv < GL_tentacles_level_max + 1; lv++)
	{
		no_pts = GL_tentacles_size[lv];

		for (int i = 0; i < no_pts; i++)
		{
			ptr_node = GL_tentacles[lv][i];
			TOTAL_MEMORY_CELDAS += 4 * ptr_node->cell_cap * sizeof(int);

			if (lmin < lmax)
			{
				TOTAL_MEMORY_CELL_STRUCT += ptr_node->box_cap * sizeof(struct cell_struct);
				TOTAL_MEMORY_CELL_STRUCT += ptr_node->cell_struct_old_cap * sizeof(struct cell_struct);
				if (ptr_node->ptr_cell_struct != NULL)
				{
					for (int j = 0; j < ptr_node->box_cap; j++)
					{
						TOTAL_MEMORY_CELL_STRUCT += ptr_node->ptr_cell_struct[j].ptcl_cap * sizeof(int);
					}
				}

				if (ptr_node->ptr_cell_struct_old != NULL)
				{
					for (int j = 0; j < ptr_node->cell_struct_old_cap; j++)
					{
						TOTAL_MEMORY_CELL_STRUCT += ptr_node->ptr_cell_struct_old[j].ptcl_cap * sizeof(int);
					}
				}
			}
			TOTAL_MEMORY_CAJAS += 2 * ptr_node->box_cap * (sizeof(int) + sizeof(vtype)); // Boxes and mass boxes
			if (lv != 0)
			{
				TOTAL_MEMORY_GRID_POINTS += 4 * (ptr_node->grid_bder_cap + ptr_node->grid_intr_cap) * sizeof(int); // Grid interior and border points
			}
			else
			{
				TOTAL_MEMORY_GRID_POINTS += (ptr_node->grid_bder_cap + ptr_node->grid_intr_cap) * sizeof(int); // Grid interior and border points
			}

			TOTAL_MEMORY_GRID_PROPERTIES += 6 * ptr_node->grid_properties_cap * sizeof(vtype); // Grid properties, accelerations, potential and density
			TOTAL_MEMORY_AUX += ptr_node->zones_cap * sizeof(int *) + ptr_node->cell_ref_cap * sizeof(int);
			for (int j = 0; j < ptr_node->zones_cap; j++)
			{
				TOTAL_MEMORY_AUX += ptr_node->ptr_zone_cap[j] * sizeof(int);
			}
			TOTAL_MEMORY_AUX += ptr_node->aux_idx_cap * sizeof(int);
			TOTAL_MEMORY_AUX += 4 * ptr_node->links_cap * sizeof(int);
		}
		TOTAL_MEMORY_NODES += no_pts * sizeof(struct node);
	}

	// Stack of memory pool
	int cntr_nodes_memory_pool = 0;
	ptr_node = GL_pool_node_start;
	while (ptr_node != NULL)
	{
		cntr_nodes_memory_pool++;
		TOTAL_MEMORY_STACK += 4 * ptr_node->cell_cap * sizeof(int);
		TOTAL_MEMORY_STACK += ptr_node->box_cap * sizeof(struct cell_struct);

		TOTAL_MEMORY_STACK += ptr_node->box_cap * sizeof(struct cell_struct);
		TOTAL_MEMORY_STACK += ptr_node->cell_struct_old_cap * sizeof(struct cell_struct);
		if (ptr_node->ptr_cell_struct != NULL)
		{
			for (int j = 0; j < ptr_node->box_cap; j++)
			{
				TOTAL_MEMORY_STACK += ptr_node->ptr_cell_struct[j].ptcl_cap * sizeof(int);
			}
		}

		if (ptr_node->ptr_cell_struct_old != NULL)
		{
			for (int j = 0; j < ptr_node->cell_struct_old_cap; j++)
			{
				TOTAL_MEMORY_STACK += ptr_node->ptr_cell_struct_old[j].ptcl_cap * sizeof(int);
			}
		}

		TOTAL_MEMORY_STACK += 2 * ptr_node->box_cap * (sizeof(int) + sizeof(vtype));				 // Boxes and mass boxes
		TOTAL_MEMORY_STACK += 4 * (ptr_node->grid_bder_cap + ptr_node->grid_intr_cap) * sizeof(int); // Grid interior and border points
		TOTAL_MEMORY_STACK += 6 * ptr_node->grid_properties_cap * sizeof(vtype);					 // Grid properties, accelerations, potential and density
		TOTAL_MEMORY_STACK += ptr_node->zones_cap * sizeof(int *) + ptr_node->cell_ref_cap * sizeof(int);
		for (int j = 0; j < ptr_node->zones_cap; j++)
		{
			TOTAL_MEMORY_STACK += ptr_node->ptr_zone_cap[j] * sizeof(int);
		}
		TOTAL_MEMORY_STACK += ptr_node->aux_idx_cap * sizeof(int);
		TOTAL_MEMORY_STACK += sizeof(struct node);
		if (ptr_node == GL_pool_node_end)
		{
			ptr_node = NULL;
		}
		else
		{
			ptr_node = ptr_node->ptr_pt;
		}
	}

	return cntr_nodes_memory_pool;
}

static void free_memory_pool(void)
{
	// Stack of memory pool
	struct node *ptr_node = GL_pool_node_start;
	while (ptr_node != NULL)
	{
		//** >> Boxes **/
		free(ptr_node->ptr_box);
		free(ptr_node->ptr_box_old);

		//** >> Cells in the node **/
		free(ptr_node->ptr_cell_idx_x);
		free(ptr_node->ptr_cell_idx_y);
		free(ptr_node->ptr_cell_idx_z);
		free(ptr_node->ptr_box_idx);

		// //** >> Struct of cells (Particles and cell mass)
		for (int i = 0; i < ptr_node->box_cap; i++)
		{
			free(ptr_node->ptr_cell_struct[i].ptr_ptcl);
		}
		free(ptr_node->ptr_cell_struct);
		if (ptr_node->ptr_cell_struct_old !=NULL)
		{
			for (int i = 0; i < ptr_node->cell_struct_old_cap; i++)
			{
				if (ptr_node->ptr_cell_struct_old[i].ptr_ptcl != NULL)
				{
					free(ptr_node->ptr_cell_struct_old[i].ptr_ptcl);
				}
			}
			free(ptr_node->ptr_cell_struct_old);
		}

		//** >> Grid points **/
		free(ptr_node->ptr_intr_grid_cell_idx_x);
		free(ptr_node->ptr_intr_grid_cell_idx_y);
		free(ptr_node->ptr_intr_grid_cell_idx_z);
		free(ptr_node->ptr_intr_grid_idx);

		free(ptr_node->ptr_bder_grid_cell_idx_x);
		free(ptr_node->ptr_bder_grid_cell_idx_y);
		free(ptr_node->ptr_bder_grid_cell_idx_z);
		free(ptr_node->ptr_bder_grid_idx);

		//** >> Potential, acceleration and density of the grid **/
		free(ptr_node->ptr_pot);
		free(ptr_node->ptr_pot_old);
		free(ptr_node->ptr_ax);
		free(ptr_node->ptr_ay);
		free(ptr_node->ptr_az);
		free(ptr_node->ptr_d);

		//** >> Auxililary arrays to go from old box to new box **/
		free(ptr_node->ptr_cell_ref);
		for(int i = 0; i< ptr_node->zones_cap;i++)
		{
			if (ptr_node->pptr_zones[i] !=NULL)
			{
				free(ptr_node->pptr_zones[i]);
			}
		}
		free(ptr_node->ptr_zone_cap);
		free(ptr_node->ptr_zone_size);
		free(ptr_node->pptr_zones);

		free(ptr_node->ptr_aux_idx);

		//** >> Links in Tree adaptation **/
		free(ptr_node->ptr_links_old_ord_old);
		free(ptr_node->ptr_links_new_ord_old);
		free(ptr_node->ptr_links_old_ord_new);
		free(ptr_node->ptr_links_new_ord_new);

		if (ptr_node == GL_pool_node_end)
		{
			ptr_node = NULL;
		}
		else
		{
			ptr_node = ptr_node->ptr_pt;
		}
		free(GL_pool_node_start);
		GL_pool_node_start = ptr_node;
	}

	GL_pool_node_start = NULL;
	GL_pool_node_end = NULL;
}

static int free_nodes_voids(void)
{

	struct node *ptr_node;
	int no_pts;

	int *ptr_aux;

	for (int lv = 0; lv < GL_tentacles_level_max + 1; lv++)
	{
		no_pts = GL_tentacles_size[lv];

		for (int i = 0; i < no_pts; i++)
		{
			ptr_node = GL_tentacles[lv][i];

			for (int j = 0; j < ptr_node->box_cap; j++)
			{
				if (ptr_node->ptr_cell_struct[j].ptcl_cap > ptr_node->ptr_cell_struct[j].ptcl_size)
				{
					
					if (ptr_node->ptr_cell_struct[j].ptcl_size == 0)
					{
						free(ptr_node->ptr_cell_struct[j].ptr_ptcl);
						ptr_node->ptr_cell_struct[j].ptr_ptcl = NULL;
					}
					else
					{
						ptr_aux = NULL;
						ptr_aux = (int *)realloc(ptr_node->ptr_cell_struct[j].ptr_ptcl, ptr_node->ptr_cell_struct[j].ptcl_size * sizeof(int));
						if (ptr_aux == NULL)
						{
							printf("\nError in the realocation of ptr_int\n");
							return _FAILURE_;
						}
						else
						{
							ptr_node->ptr_cell_struct[j].ptr_ptcl = ptr_aux;
						}
					}
					
					ptr_node->ptr_cell_struct[j].ptcl_cap = ptr_node->ptr_cell_struct[j].ptcl_size;
				}
			}

			if (ptr_node->ptr_cell_struct_old != NULL)
			{
				for (int j = 0; j < ptr_node->cell_struct_old_cap; j++)
				{
					if (ptr_node->ptr_cell_struct_old[j].ptr_ptcl != NULL)
					{
						free(ptr_node->ptr_cell_struct_old[j].ptr_ptcl);
						ptr_node->ptr_cell_struct_old[j].ptr_ptcl = NULL;
					}
				}
				free(ptr_node->ptr_cell_struct_old);
				ptr_node->ptr_cell_struct_old = NULL;
				ptr_node->cell_struct_old_cap = 0;
			}
		}
	}

	return _SUCCESS_;
}

int garbage_collector(void)
{

	TOTAL_MEMORY_NODES = 0;
	TOTAL_MEMORY_CELDAS = 0;
	TOTAL_MEMORY_PARTICLES = 0;
	TOTAL_MEMORY_CELL_STRUCT = 0;
	TOTAL_MEMORY_CAJAS = 0;
	TOTAL_MEMORY_GRID_POINTS = 0;
	TOTAL_MEMORY_GRID_PROPERTIES = 0;
	TOTAL_MEMORY_AUX = 0;
	TOTAL_MEMORY_TENTACLES = 0;
	TOTAL_MEMORY_STACK = 0;

	int cntr_nodes_memory_pool;
	cntr_nodes_memory_pool = computing_memory();

	double sum = TOTAL_MEMORY_NODES + TOTAL_MEMORY_CELDAS + TOTAL_MEMORY_PARTICLES + TOTAL_MEMORY_CELL_STRUCT + TOTAL_MEMORY_CAJAS + TOTAL_MEMORY_GRID_POINTS + TOTAL_MEMORY_GRID_PROPERTIES + TOTAL_MEMORY_AUX + TOTAL_MEMORY_TENTACLES + TOTAL_MEMORY_STACK;

	printf("\n\n%sMEMORY [MB]:%s\n\n", KRED, KNRM);

	double Total_memory[10] = {TOTAL_MEMORY_NODES,
							   TOTAL_MEMORY_CELDAS,
							   TOTAL_MEMORY_PARTICLES,
							   TOTAL_MEMORY_CELL_STRUCT,
							   TOTAL_MEMORY_CAJAS,
							   TOTAL_MEMORY_GRID_POINTS,
							   TOTAL_MEMORY_GRID_PROPERTIES,
							   TOTAL_MEMORY_AUX,
							   TOTAL_MEMORY_TENTACLES,
							   TOTAL_MEMORY_STACK};

	char Memory_names[50][100] = {
		"Nodes",
		"Celdas",
		"Particulas",
		"Cell struct",
		"Cajas",
		"Grid Points",
		"Grid Properties",
		"Auxiliary",
		"Tentacles",
		"Memory Pool",
	};

	for (int i = 0; i < 10; i++)
	{
		printf("%s = %f ~ %.1f %%\n", Memory_names[i], Total_memory[i] / 1000000, Total_memory[i] * 100 / sum);
		
	}
	printf("Nodes in memory pool = %d\n", cntr_nodes_memory_pool);

	printf("\n%sTOTAL = %f MB%s\n\n", KMAG, sum / 1000000, KNRM);

	free_memory_pool();


	if (free_nodes_voids() == _FAILURE_)
	{
		printf("\n\n Error running free_nodes_voids() function\n\n");
		return _FAILURE_;
	}

	return _SUCCESS_;
}