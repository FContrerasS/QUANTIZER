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

int reset(void)
{
  // Global Particles properties

  // Head and brances nodes

  //* >> Global particle accelerations *//
  for (int i = 0; i < GL_no_ptcl_final; i++)
  {
    GL_ptcl_ax[i] = 0.0;
    GL_ptcl_ay[i] = 0.0;
    GL_ptcl_az[i] = 0.0;
  }

  //* >> Head node *//
  struct node *ptr_head = GL_ptr_tree;

  //* >> Potential, Acceleration and density of the grid *//
  int cap = (ptr_head->box_real_dim_x + 1) * (ptr_head->box_real_dim_y + 1) * (ptr_head->box_real_dim_z + 1);

  for (int i = 0; i < cap; i++)
  {
    ptr_head->ptr_d[i] = 0.0;
  }

  // for (int i = 0; i < cap; i++)
  // {
  // 	ptr_head->ptr_pot[i] = 0;
  // }

  //* >> Initial values for the potential and acceleration *//
  //* Here we are modified the boundary values of the potential and acceleration
  // initial_potential_and_acceleration_head(ptr_head);

  return _SUCCESS_;
}
