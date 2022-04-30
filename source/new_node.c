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

struct node *new_node()
{
	struct node *ptr_node = NULL;

	if (GL_pool_node_start == NULL)
	{
		ptr_node = (struct node *)malloc(sizeof(struct node));
		initialize_node(ptr_node);
	}
	else
	{
		ptr_node = GL_pool_node_start;
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