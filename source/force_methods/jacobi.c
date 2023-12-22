/*
 * jacobi.c
 *
 * Solve the Discrete Poisson equation using the Jacobi method
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

#include "jacobi.h"

static void compute_jacobi(vtype *phi_old, vtype *phi_new, const vtype *rho,
                           int gridsize, int iter_max)
{

  vtype *pointer_auxiliary;
  vtype aux_pot;
  vtype C = -(1.0 / (gridsize - 1)) * (1.0 / (gridsize - 1));
  int gridsize_pow_2 = gridsize * gridsize;
  int C1 = (gridsize - 1);
  int C2 = gridsize * (gridsize - 1);
  int C3 = gridsize_pow_2 * (gridsize - 1);
  int Index;

  for (int k = 0; k < iter_max; k++)
  {
    for (int lz = gridsize_pow_2; lz < C3; lz += gridsize_pow_2)
    {
      for (int ly = gridsize; ly < C2; ly += gridsize)
      {
        for (int lx = 1; lx < C1; lx++)
        {
          Index = lz + ly + lx;
          aux_pot = _Onesixth_ * (C * rho[Index] + phi_old[Index + 1] + phi_old[Index - 1] + phi_old[Index + gridsize] + phi_old[Index - gridsize] + phi_old[Index + gridsize_pow_2] + phi_old[Index - gridsize_pow_2]);
          phi_new[Index] = _w_SOR_HEAD_ * aux_pot + (1.0 - _w_SOR_HEAD_) * phi_new[Index];
        }
      }
    }
    pointer_auxiliary = phi_old;
    phi_old = phi_new;
    phi_new = pointer_auxiliary;
  }
}

void jacobi(vtype *phi, const vtype *rho, int gridsize, int iter_max)
{

  vtype *phi_copy;
  int gridsize_pow_3 = gridsize * gridsize * gridsize;

  // Copy of phi
  phi_copy = (vtype *)malloc(gridsize_pow_3 * sizeof(vtype));
  memcpy(phi_copy, phi, gridsize_pow_3 * sizeof(vtype));

  // It depends on what buff update first
  if (iter_max % 2 == 0)
  {
    compute_jacobi(phi, phi_copy, rho, gridsize, iter_max);
  }
  else
  {
    compute_jacobi(phi_copy, phi, rho, gridsize, iter_max + 1);
  }
  free(phi_copy);
}
