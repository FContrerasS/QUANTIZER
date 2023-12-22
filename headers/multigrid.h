/*
 * multigrid.h
 *
 * Header file with the most multigrid parameters declared
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

#ifndef __MULTIGRID__
#define __MULTIGRID__

#include "common.h"
#include "jacobi.h"
#include "gauss_saidel.h"
#include "conjugate_gradient_multigrid.h"

void multigrid(vtype *phi, const vtype *rho, int gridsize, int type_multigrid);

#endif
