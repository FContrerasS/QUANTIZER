/*
 * ptcl_idx_to_box_idx.h
 *
 * Header file of the ptcl_idx_to_box_idx.c source file
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

#ifndef __PTCLIDXTOBOXIDX__
#define __PTCLIDXTOBOXIDX__

#include "common.h"

struct node;

int ptcl_idx_to_box_idx(struct node *ptr_node, int ptcl_idx);

#endif

