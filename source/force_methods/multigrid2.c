/*
 * multigrid2.c
 *
 * Solve the Discrete Poisson equation using the multigrid2 method
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

#include "multigrid2.h"

static void multigrid_restriction2(const vtype *phi, vtype *cphi, int cgridsize)
{

    int gridsize = 2 * cgridsize - 1;
    int gridsize_pow_2 = gridsize * gridsize;
    int cgridsize_pow_2 = cgridsize * cgridsize;
    int Index_A, Index_B;

    for (int lz = 1; lz < cgridsize - 1; lz++)
    {
        for (int ly = 1; ly < cgridsize - 1; ly++)
        {
            for (int lx = 1; lx < cgridsize - 1; lx++)
            {
                Index_A = lz * cgridsize_pow_2 + ly * cgridsize + lx;
                Index_B = 2 * lz * gridsize_pow_2 + 2 * ly * gridsize + 2 * lx;
                cphi[Index_A] = phi[Index_B];
            }
        }
    }
}

static void multigrid_projection2(const vtype *dphic, vtype *dphi, int cgridsize)
{

    int gridsize = 2 * cgridsize - 1;
    int gridsize_pow_2 = gridsize * gridsize;
    int cgridsize_pow_2 = cgridsize * cgridsize;

    int fine_idx;

    int coarse_idx_x0;
    int coarse_idx_x1;
    int coarse_idx_y0;
    int coarse_idx_y1;
    int coarse_idx_z0;
    int coarse_idx_z1;

    int coarse_idx_i0_j0_k0;
    int coarse_idx_i1_j0_k0;
    int coarse_idx_i0_j1_k0;
    int coarse_idx_i1_j1_k0;
    int coarse_idx_i0_j0_k1;
    int coarse_idx_i1_j0_k1;
    int coarse_idx_i0_j1_k1;
    int coarse_idx_i1_j1_k1;

    for (int lz = 1; lz < gridsize - 1; lz+=2)
    {
        coarse_idx_z0 = lz >> 1;
        coarse_idx_z1 = (lz + 1) >> 1;
        for (int ly = 1; ly < gridsize - 1; ly+=2)
        {
            coarse_idx_y0 = ly >> 1;
            coarse_idx_y1 = coarse_idx_y0 + 1;
            for (int lx = 1; lx < gridsize - 1; lx+=2)
            {

                coarse_idx_x0 = lx >> 1;
                coarse_idx_x1 = coarse_idx_x0 + 1;
                fine_idx = lx + ly * gridsize + lz * gridsize_pow_2;

                coarse_idx_i0_j0_k0 = coarse_idx_x0 + coarse_idx_y0 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;
                coarse_idx_i1_j0_k0 = coarse_idx_x1 + coarse_idx_y0 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;
                coarse_idx_i0_j1_k0 = coarse_idx_x0 + coarse_idx_y1 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;
                coarse_idx_i1_j1_k0 = coarse_idx_x1 + coarse_idx_y1 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;
                coarse_idx_i0_j0_k1 = coarse_idx_x0 + coarse_idx_y0 * cgridsize + coarse_idx_z1 * cgridsize_pow_2;
                coarse_idx_i1_j0_k1 = coarse_idx_x1 + coarse_idx_y0 * cgridsize + coarse_idx_z1 * cgridsize_pow_2;
                coarse_idx_i0_j1_k1 = coarse_idx_x0 + coarse_idx_y1 * cgridsize + coarse_idx_z1 * cgridsize_pow_2;
                coarse_idx_i1_j1_k1 = coarse_idx_x1 + coarse_idx_y1 * cgridsize + coarse_idx_z1 * cgridsize_pow_2;

                dphi[fine_idx] = 0.125 * (dphic[coarse_idx_i0_j0_k0] + dphic[coarse_idx_i1_j0_k0] +
                                          dphic[coarse_idx_i0_j1_k0] + dphic[coarse_idx_i1_j1_k0] +
                                          dphic[coarse_idx_i0_j0_k1] + dphic[coarse_idx_i1_j0_k1] +
                                          dphic[coarse_idx_i0_j1_k1] + dphic[coarse_idx_i1_j1_k1]);

                // There are 4 cycles. The size are 1,6,12 and 8.
                // They are organized according the distance to the coarse grid point.

            }

            for (int lx = 2; lx < gridsize - 1; lx += 2)
            {

                coarse_idx_x0 = lx >> 1;
                fine_idx = lx + ly * gridsize + lz * gridsize_pow_2;

                coarse_idx_i0_j0_k0 = coarse_idx_x0 + coarse_idx_y0 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;
                coarse_idx_i0_j1_k0 = coarse_idx_x0 + coarse_idx_y1 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;
                coarse_idx_i0_j0_k1 = coarse_idx_x0 + coarse_idx_y0 * cgridsize + coarse_idx_z1 * cgridsize_pow_2;
                coarse_idx_i0_j1_k1 = coarse_idx_x0 + coarse_idx_y1 * cgridsize + coarse_idx_z1 * cgridsize_pow_2;

                dphi[fine_idx] = 0.25 * (dphic[coarse_idx_i0_j0_k0] +
                                        dphic[coarse_idx_i0_j1_k0]  +
                                        dphic[coarse_idx_i0_j0_k1] +
                                        dphic[coarse_idx_i0_j1_k1] );

                // There are 4 cycles. The size are 1,6,12 and 8.
                // They are organized according the distance to the coarse grid point.
            }
        }
        for (int ly = 2; ly < gridsize - 1; ly+=2)
        {
            coarse_idx_y0 = ly >> 1;
            for (int lx = 1; lx < gridsize - 1; lx += 2)
            {

                coarse_idx_x0 = lx >> 1;
                coarse_idx_x1 = coarse_idx_x0 + 1;
                fine_idx = lx + ly * gridsize + lz * gridsize_pow_2;

                coarse_idx_i0_j0_k0 = coarse_idx_x0 + coarse_idx_y0 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;
                coarse_idx_i1_j0_k0 = coarse_idx_x1 + coarse_idx_y0 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;
                coarse_idx_i0_j0_k1 = coarse_idx_x0 + coarse_idx_y0 * cgridsize + coarse_idx_z1 * cgridsize_pow_2;
                coarse_idx_i1_j0_k1 = coarse_idx_x1 + coarse_idx_y0 * cgridsize + coarse_idx_z1 * cgridsize_pow_2;

                dphi[fine_idx] = 0.25 * (dphic[coarse_idx_i0_j0_k0] + dphic[coarse_idx_i1_j0_k0] +
                                          dphic[coarse_idx_i0_j0_k1] + dphic[coarse_idx_i1_j0_k1]);

                // There are 4 cycles. The size are 1,6,12 and 8.
                // They are organized according the distance to the coarse grid point.
            }

            for (int lx = 2; lx < gridsize - 1; lx += 2)
            {

                coarse_idx_x0 = lx >> 1;
                fine_idx = lx + ly * gridsize + lz * gridsize_pow_2;

                coarse_idx_i0_j0_k0 = coarse_idx_x0 + coarse_idx_y0 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;
                coarse_idx_i0_j0_k1 = coarse_idx_x0 + coarse_idx_y0 * cgridsize + coarse_idx_z1 * cgridsize_pow_2;

                dphi[fine_idx] = 0.5 * (dphic[coarse_idx_i0_j0_k0] +
                                         dphic[coarse_idx_i0_j0_k1] );

                // There are 4 cycles. The size are 1,6,12 and 8.
                // They are organized according the distance to the coarse grid point.
            }
        }
    }

    for (int lz = 2; lz < gridsize - 1; lz+=2)
    {
        coarse_idx_z0 = lz >> 1;
        for (int ly = 1; ly < gridsize - 1; ly += 2)
        {
            coarse_idx_y0 = ly >> 1;
            coarse_idx_y1 = coarse_idx_y0 + 1;
            for (int lx = 1; lx < gridsize - 1; lx += 2)
            {

                coarse_idx_x0 = lx >> 1;
                coarse_idx_x1 = coarse_idx_x0 + 1;
                fine_idx = lx + ly * gridsize + lz * gridsize_pow_2;

                coarse_idx_i0_j0_k0 = coarse_idx_x0 + coarse_idx_y0 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;
                coarse_idx_i1_j0_k0 = coarse_idx_x1 + coarse_idx_y0 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;
                coarse_idx_i0_j1_k0 = coarse_idx_x0 + coarse_idx_y1 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;
                coarse_idx_i1_j1_k0 = coarse_idx_x1 + coarse_idx_y1 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;

                dphi[fine_idx] = 0.25 * (dphic[coarse_idx_i0_j0_k0] + dphic[coarse_idx_i1_j0_k0] +
                                          dphic[coarse_idx_i0_j1_k0] + dphic[coarse_idx_i1_j1_k0] );

                // There are 4 cycles. The size are 1,6,12 and 8.
                // They are organized according the distance to the coarse grid point.
            }

            for (int lx = 2; lx < gridsize - 1; lx += 2)
            {

                coarse_idx_x0 = lx >> 1;
                fine_idx = lx + ly * gridsize + lz * gridsize_pow_2;

                coarse_idx_i0_j0_k0 = coarse_idx_x0 + coarse_idx_y0 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;
                coarse_idx_i0_j1_k0 = coarse_idx_x0 + coarse_idx_y1 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;

                dphi[fine_idx] = 0.5 * (dphic[coarse_idx_i0_j0_k0] +
                                         dphic[coarse_idx_i0_j1_k0] );

                // There are 4 cycles. The size are 1,6,12 and 8.
                // They are organized according the distance to the coarse grid point.
            }
        }
        for (int ly = 2; ly < gridsize - 1; ly += 2)
        {
            coarse_idx_y0 = ly >> 1;
            for (int lx = 1; lx < gridsize - 1; lx += 2)
            {

                coarse_idx_x0 = lx >> 1;
                coarse_idx_x1 = coarse_idx_x0 + 1;
                fine_idx = lx + ly * gridsize + lz * gridsize_pow_2;

                coarse_idx_i0_j0_k0 = coarse_idx_x0 + coarse_idx_y0 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;
                coarse_idx_i1_j0_k0 = coarse_idx_x1 + coarse_idx_y0 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;

                dphi[fine_idx] = 0.5 * (dphic[coarse_idx_i0_j0_k0] + dphic[coarse_idx_i1_j0_k0] );

                // There are 4 cycles. The size are 1,6,12 and 8.
                // They are organized according the distance to the coarse grid point.
            }

            for (int lx = 2; lx < gridsize - 1; lx += 2)
            {
                coarse_idx_x0 = lx >> 1;
                fine_idx = lx + ly * gridsize + lz * gridsize_pow_2;

                coarse_idx_i0_j0_k0 = coarse_idx_x0 + coarse_idx_y0 * cgridsize + coarse_idx_z0 * cgridsize_pow_2;
                
                dphi[fine_idx] =  dphic[coarse_idx_i0_j0_k0];
                // There are 4 cycles. The size are 1,6,12 and 8.
                // They are organized according the distance to the coarse grid point.
            }
        }
    }
}

void multigrid2(vtype *phi, vtype *rho, int type_multigrid)
{

    // Testing V cycle

    int grid_idx;

    // Filling arrays dphic and rest
    int size, size_pow2, size_pow3, size_coarse;

    vtype one_over_H2;

    int C1, C2, C3;

    size = ((1 << lmin) + 1);
    size_pow3 = size * size * size; // size^3

    //** >> Pre - Smoothing **/
    if (solverPreS == 0)
    {
        gauss_saidel(phi, rho, size, _NiterPreS_);
    }
    else if (solverPreS == 1)
    {
        jacobi(phi, rho, size, _NiterPreS_);
    }
    else if (solverPreS == 2)
    {
        conjugate_gradient_multigrid(phi, rho, size,_NiterPreS_);
    }

    pp_phixx[0] = phi;
    pp_rhoxx[0] = rho;

    for (int lv = 1; lv < lmin - 1; lv++)
    {
        size = ((1 << (lmin - lv)) + 1);
        size_pow3 = size * size * size; // size^3
        memcpy(pp_phixx[lv], zeros_xx[lv], size_pow3 * sizeof(vtype));
    }

    for (int lv = 0; lv < lmin - 2; lv++)
    {
        size = ((1 << (lmin - lv)) + 1);
        size_pow2 = size * size;
        size_pow3 = size * size * size;
        size_coarse = (size + 1) / 2;
        one_over_H2 = (size - 1) * (size - 1);

        C1 = (size - 1);
        C2 = size * (size - 1);
        C3 = size_pow2 * (size - 1);

        //** >> Finding rest **/
        // rest_i = rho_i - Lap[phi_i]
        for (int k = size_pow2; k < C3; k += size_pow2)
        {
            for (int j = size; j < C2; j += size)
            {
                for (int i = 1; i < C1; i++)
                {
                    grid_idx = k + j + i;
                    pp_restxx[lv][grid_idx] = pp_rhoxx[lv][grid_idx] + one_over_H2 * (6 * pp_phixx[lv][grid_idx] - pp_phixx[lv][grid_idx + 1] - pp_phixx[lv][grid_idx - 1] - pp_phixx[lv][grid_idx + size] - pp_phixx[lv][grid_idx - size] - pp_phixx[lv][grid_idx + size_pow2] - pp_phixx[lv][grid_idx - size_pow2]);
                }
            }
        }

        //** >> Restriction: Interpoling fine to coarse **/
        multigrid_restriction2(pp_restxx[lv], pp_rhoxx[lv + 1], size_coarse);

        //** >> Smoothing at next level **/
        if (solverPreS == 0)
        {
            gauss_saidel(pp_phixx[lv + 1], pp_rhoxx[lv + 1], size_coarse, _NiterPreS_);
        }
        else if (solverPreS == 1)
        {
            jacobi(pp_phixx[lv + 1], pp_rhoxx[lv + 1], size_coarse, _NiterPreS_);
        }
        else if (solverPreS == 2)
        {
            conjugate_gradient_multigrid(pp_phixx[lv + 1], pp_rhoxx[lv + 1], size_coarse, _NiterPreS_);
        }
    }

    size = ((1 << 3) + 1);
    size_coarse = (size + 1) / 2;
    if (solverfinddphic == 0)
    {
        gauss_saidel(pp_phixx[lmin - 2], pp_rhoxx[lmin - 2], size_coarse, 2);
    }
    else if (solverfinddphic == 1)
    {
        jacobi(pp_phixx[lmin - 2], pp_rhoxx[lmin - 2], size_coarse, _Niterfinddphic_);
    }
    else if (solverfinddphic == 2)
    {
        conjugate_gradient_multigrid(pp_phixx[lmin - 2], pp_rhoxx[lmin - 2], size_coarse, _Niterfinddphic_);
    }

    for (int lv = lmin - 3; lv > -1; lv--)
    {
        size = ((1 << (lmin - lv)) + 1);
        size_pow2 = size * size;
        size_coarse = (size + 1) / 2;

        C1 = (size - 1);
        C2 = size * (size - 1);
        C3 = size_pow2 * (size - 1);

        //** >> Projection: Interpolation coarse to fine **/
        multigrid_projection2(pp_phixx[lv + 1], pp_dphixx[lv], size_coarse);

        //** >> Adding dphi in solution phi = phi + dphi **/
        for (int k = size_pow2; k < C3; k += size_pow2)
        {
            for (int j = size; j < C2; j += size)
            {
                for (int i = 1; i < C1; i++)
                {
                    grid_idx = k + j + i;
                    pp_phixx[lv][grid_idx] += pp_dphixx[lv][grid_idx];
                }
            }
        }

        //** >> Post - Smoothing **/
        if (solverPostS == 0)
        {
            gauss_saidel(pp_phixx[lv], pp_rhoxx[lv], size, _NiterPostS_);
        }
        else if (solverPostS == 1)
        {
            jacobi(pp_phixx[lv], pp_rhoxx[lv], size, _NiterPostS_);
        }
        else if (solverPostS == 2)
        {
            conjugate_gradient_multigrid(pp_phixx[lv], pp_rhoxx[lv], size, _NiterPostS_);
        }
    }
}
