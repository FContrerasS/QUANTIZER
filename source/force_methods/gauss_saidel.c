/*
 * gauss_saidel.c
 *
 * Solve the Discrete Poisson equation using the Gauss-Saidel method
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

#include "gauss_saidel.h"

void gauss_saidel(vtype *phi, const vtype *rho, int gridsize, int iter_max)
{
  vtype H_pow2 = (1.0L / (gridsize - 1)) * (1.0L / (gridsize - 1));
  vtype aux_pot;
  int grid_idx;

  int gridsize_pow_2 = gridsize * gridsize;
  int C1 = (gridsize - 1);
  int C2 = gridsize * (gridsize - 1);
  int C3 = gridsize_pow_2 * (gridsize - 1);

  for (int iter = 0; iter < iter_max; iter++)
  {
    for (int k = gridsize_pow_2; k < C3; k += gridsize_pow_2)
    {
      for (int j = gridsize; j < C2; j += gridsize)
      {
        for (int i = 1; i < C1; i++)
        {
          grid_idx = k + j + i;
          aux_pot = _Onesixth_ * (-H_pow2 * rho[grid_idx] + phi[grid_idx + 1] + phi[grid_idx - 1] + phi[grid_idx + gridsize] + phi[grid_idx - gridsize] + phi[grid_idx + gridsize_pow_2] + phi[grid_idx - gridsize_pow_2]);
          phi[grid_idx] = _w_SOR_HEAD_ * aux_pot + (1.0 - _w_SOR_HEAD_) * phi[grid_idx];
        }
      }
    }
  }
}
