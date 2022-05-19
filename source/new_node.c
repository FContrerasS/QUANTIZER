/*
 * new_node.c
 *
 * Return a pointer to a node using the stack of the memory pool 
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

#include "new_node.h"


void re_initialize_node(struct node *ptr_node)
{
  	int box_idx_node;

	//** >> Boxes **/
	ptr_node->box_min_x = INT_MAX; // Already minimal box value index in the real local space at the dimension X
	ptr_node->box_min_y = INT_MAX; // Already minimal box value index in the real local space at the dimension Y
	ptr_node->box_min_z = INT_MAX; // Already minimal box value index in the real local space at the dimension Z
	ptr_node->box_max_x = 0;	   // Already maximum box value index in the real local space at the dimension X
	ptr_node->box_max_y = 0;	   // Already maximum box value index in the real local space at the dimension Y
	ptr_node->box_max_z = 0;	   // Already maximum box value index in the real local space at the dimension Z
	ptr_node->box_check_fit = true; // Check if the new box will fit in the old one

	//** >> Struct of cells (Particles and cell mass)
	for (int i = 0; i < ptr_node->cell_size; i++)
	{
		box_idx_node = ptr_node->ptr_box_idx[i];
		ptr_node->ptr_cell_struct[box_idx_node].cell_mass = 0;
		ptr_node->ptr_cell_struct[box_idx_node].ptcl_size = 0;
	}

	//** >> Cells in the node **/
	ptr_node->cell_size = 0; // Number of existing cells in the node

	//** >> Total mass in the node **/
	ptr_node->local_mass = 0;

	//** >> Grid points **/
    ptr_node->grid_intr_size = 0; // Number of existing grid interior points in the block
    ptr_node->grid_bder_size = 0; // Number of existing grid border points in the block

	//** >> Potential, acceleration and density of the grid **/

  	for(int i = 0 ; i< ptr_node->grid_properties_cap; i++)
  	{
  	  ptr_node->ptr_d[i] = 0;
  	}

	//** >> Tree structure **/
	ptr_node->chn_size = 0;	   // Number of children of the node

	//** >> Auxililary arrays to go from old box to new box **/
	ptr_node->cell_ref_size = 0;	// Number of cells to be refined

  	for (int j = 0; j < ptr_node->zones_size; j++)
  	{
  	  ptr_node->ptr_zone_size[j] = 0;
  	}
	ptr_node->zones_size = 0;		// Number of refined zones in the node

}

struct node* new_node()
{
	struct node *ptr_node;

	if (GL_pool_node_start == NULL)
	{
		ptr_node = (struct node *)malloc(sizeof(struct node));
		initialize_node(ptr_node);
	}
	else
	{
		ptr_node = GL_pool_node_start;
		re_initialize_node(ptr_node);
		if (GL_pool_node_start == GL_pool_node_end)
		{
			GL_pool_node_start = NULL;
			GL_pool_node_end = NULL;
		}
		else
		{
			GL_pool_node_start = ptr_node->ptr_pt;
		}
	}

	return ptr_node;
}
