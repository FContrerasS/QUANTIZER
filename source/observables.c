/*
 * observables.c
 *
 * Compute the energies of the simulation in a time-step given
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

#include "observables.h"

void observables(vtype *energies)
{

    energies[0] = 0; // Kinetic
    energies[1] = 0; // Potential
    energies[2] = 0; // Total

    vtype distance_x, distance_y, distance_z; // Axial distance between 2 particles
    vtype distance;                           // Distance between 2 particles
    vtype particle_v_pow_2;                   // Particle velocity to the power of 2

    vtype _one_over_User_sqrt_BoxSize_ = 1.0L / sqrt(_User_BoxSize_);
    // Returning to user units
    for (int i = 0; i < GL_no_ptcl; i++)
    {
        GL_ptcl_x[i] = (GL_ptcl_x[i] - 0.5) * _User_BoxSize_;
        GL_ptcl_y[i] = (GL_ptcl_y[i] - 0.5) * _User_BoxSize_;
        GL_ptcl_z[i] = (GL_ptcl_z[i] - 0.5) * _User_BoxSize_;

        GL_ptcl_vx[i] *= _one_over_User_sqrt_BoxSize_;
        GL_ptcl_vy[i] *= _one_over_User_sqrt_BoxSize_;
        GL_ptcl_vz[i] *= _one_over_User_sqrt_BoxSize_;
    }

    for (int i = 0; i < GL_no_ptcl; i++)
    {
        // Kinetic Energy
        particle_v_pow_2 = GL_ptcl_vx[i] * GL_ptcl_vx[i] + GL_ptcl_vy[i] * GL_ptcl_vy[i] + GL_ptcl_vz[i] * GL_ptcl_vz[i];
        // Adding Kinetic energy of the particle
        energies[0] += GL_ptcl_mass[i] * particle_v_pow_2;
        for (int j = 0; j < i; j++)
        {
            // Potential energy
            distance_x = GL_ptcl_x[i] - GL_ptcl_x[j];
            distance_y = GL_ptcl_y[i] - GL_ptcl_y[j];
            distance_z = GL_ptcl_z[i] - GL_ptcl_z[j];
            distance = distance_x * distance_x + distance_y * distance_y + distance_z * distance_z;
            distance = sqrt(distance);
            // Checking minimal distance accepted
            if (distance < _Min_Particle_Distance_For_Energy_Computation_)
            {
                distance = _Min_Particle_Distance_For_Energy_Computation_;
            }
            // Adding Potential energy of the particle
            energies[1] += GL_ptcl_mass[i] * GL_ptcl_mass[j] / (distance);
        }
    }
    energies[0] = energies[0] * 0.5;
    energies[1] = -_G_ * energies[1];
    energies[2] = energies[0] + energies[1];

    // Returning to code units
    for (int i = 0; i < GL_no_ptcl; i++)
    {
        GL_ptcl_x[i] = GL_ptcl_x[i] / _User_BoxSize_ + 0.5;
        GL_ptcl_y[i] = GL_ptcl_y[i] / _User_BoxSize_ + 0.5;
        GL_ptcl_z[i] = GL_ptcl_z[i] / _User_BoxSize_ + 0.5;

        GL_ptcl_vx[i] /= _one_over_User_sqrt_BoxSize_;
        GL_ptcl_vy[i] /= _one_over_User_sqrt_BoxSize_;
        GL_ptcl_vz[i] /= _one_over_User_sqrt_BoxSize_;
    }
}
