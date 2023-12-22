/*
 * conjugate_gradient_multigrid.c
 *
 * Solve the Discrete Poisson equation using the conjugate gradient method
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

#include "conjugate_gradient_multigrid.h"

void conjugate_gradient_multigrid(vtype *phi, const vtype *rho, int gridsize, int conj_grad_iter)
{
  vtype H = 1.0 / (gridsize - 1);
  vtype one_over_H_pow_2 = 1.0 / (H * H);
  int gridsize_pow2 = gridsize * gridsize;
  int gridsize_pow3 = gridsize * gridsize * gridsize;

  vtype *r, *p, *Ap;
  r = (vtype *)calloc(gridsize_pow3, sizeof(vtype));
  p = (vtype *)malloc(gridsize_pow3 * sizeof(vtype));
  Ap = (vtype *)malloc(gridsize_pow3 * sizeof(vtype));

  vtype alpha, rsold, rsnew;

  int Index;
  int C1 = (gridsize - 1);
  int C2 = gridsize * (gridsize - 1);
  int C3 = gridsize_pow2 * (gridsize - 1);

  // Computing r
  for (int lz = gridsize_pow2; lz < C3; lz += gridsize_pow2)
  {
    for (int ly = gridsize; ly < C2; ly += gridsize)
    {
      for (int lx = 1; lx < C1; lx++)
      {
        Index = lz + ly + lx;
        r[Index] = rho[Index] + one_over_H_pow_2 * (6.0 * phi[Index] - phi[Index + 1] - phi[Index - 1] - phi[Index + gridsize] - phi[Index - gridsize] - phi[Index + gridsize_pow2] - phi[Index - gridsize_pow2]);
      }
    }
  }

  // Copy r to p
  memcpy(p, r, gridsize_pow3 * sizeof(vtype));

  // Computing rsold
  rsold = 0.0;
  for (int lz = gridsize_pow2; lz < C3; lz += gridsize_pow2)
  {
    for (int ly = gridsize; ly < C2; ly += gridsize)
    {
      for (int lx = 1; lx < C1; lx++)
      {
        Index = lz + ly + lx;
        rsold += r[Index] * r[Index];
      }
    }
  }

  // For cycle over mutually conjugate vectors p
  int cycle_size = conj_grad_iter < (gridsize - 2) * (gridsize - 2) * (gridsize - 2) ? conj_grad_iter : (gridsize - 2) * (gridsize - 2) * (gridsize - 2);
  for (int i = 0; i < cycle_size; i++)
  {

    // Computing Ap
    for (int lz = gridsize_pow2; lz < C3; lz += gridsize_pow2)
    {
      for (int ly = gridsize; ly < C2; ly += gridsize)
      {
        for (int lx = 1; lx < C1; lx++)
        {
          Index = lz + ly + lx;
          Ap[Index] = -one_over_H_pow_2 * (6.0 * p[Index] - p[Index + 1] - p[Index - 1] - p[Index + gridsize] - p[Index - gridsize] - p[Index + gridsize_pow2] - p[Index - gridsize_pow2]);
        }
      }
    }

    // Computing alpha
    alpha = 0.0;
    for (int lz = gridsize_pow2; lz < C3; lz += gridsize_pow2)
    {
      for (int ly = gridsize; ly < C2; ly += gridsize)
      {
        for (int lx = 1; lx < C1; lx++)
        {
          Index = lz + ly + lx;
          alpha += p[Index] * Ap[Index];
        }
      }
    }
    alpha = rsold / alpha;

    // New solution
    for (int lz = gridsize_pow2; lz < C3; lz += gridsize_pow2)
    {
      for (int ly = gridsize; ly < C2; ly += gridsize)
      {
        for (int lx = 1; lx < C1; lx++)
        {
          Index = lz + ly + lx;
          phi[Index] += alpha * p[Index];
        }
      }
    }

    // New rest
    for (int lz = gridsize_pow2; lz < C3; lz += gridsize_pow2)
    {
      for (int ly = gridsize; ly < C2; ly += gridsize)
      {
        for (int lx = 1; lx < C1; lx++)
        {
          Index = lz + ly + lx;
          r[Index] -= alpha * Ap[Index];
        }
      }
    }

    // Computing rsnew
    rsnew = 0.0;
    for (int lz = gridsize_pow2; lz < C3; lz += gridsize_pow2)
    {
      for (int ly = gridsize; ly < C2; ly += gridsize)
      {
        for (int lx = 1; lx < C1; lx++)
        {
          Index = lz + ly + lx;
          rsnew += r[Index] * r[Index];
        }
      }
    }

    if (rsnew / rsold < 1.0e-16)
    {
      break;
    }

    // Updating p
    for (int lz = gridsize_pow2; lz < C3; lz += gridsize_pow2)
    {
      for (int ly = gridsize; ly < C2; ly += gridsize)
      {
        for (int lx = 1; lx < C1; lx++)
        {
          Index = lz + ly + lx;
          p[Index] = r[Index] + rsnew / rsold * p[Index];
        }
      }
    }

    // Updating rsold
    rsold = rsnew;
  }

  free(r);
  free(p);
  free(Ap);
}
