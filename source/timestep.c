/*
 * timestep.c
 *
 * Compute of the time-step considering velocity and/or acceleration
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

#include "timestep.h"

//* >> Local Functions
static vtype timestep_computation_vel_only_HEAD_ONLY(const struct node *ptr_node);
static vtype timestep_computation_vel_plus_accel_HEAD_ONLY(const struct node *ptr_node);
static vtype timestep_computation_vel_only(const struct node *ptr_node);
static vtype timestep_computation_vel_plus_accel(const struct node *ptr_node);

static vtype timestep_computation_vel_only_HEAD_ONLY(const struct node *ptr_node)
{

  vtype aux_v_pow2; // Auxiliary square velocity
  vtype mydt;       // The time-step of the node
  vtype myvmax;     // Maximum velocity of all particles

  vtype H = 1.0 / (1 << ptr_node->lv);

  myvmax = 0.0; // Minium velocity designated by the user

  for (int i = 0; i < GL_no_ptcl_final; i++)
  {
    //* >> Velocity at x-axis
    aux_v_pow2 = myabs(GL_ptcl_vx[i]);
    if (myvmax < aux_v_pow2)
    {
      myvmax = aux_v_pow2;
    }
    //* >> Velocity at y-axis
    aux_v_pow2 = myabs(GL_ptcl_vy[i]);
    if (myvmax < aux_v_pow2)
    {
      myvmax = aux_v_pow2;
    }
    //* >> Velocity at z-axis
    aux_v_pow2 = myabs(GL_ptcl_vz[i]);
    if (myvmax < aux_v_pow2)
    {
      myvmax = aux_v_pow2;
    }
  }

  mydt = myvmax < (_CFL_ * H / _MAX_dt_) ? _MAX_dt_ : (_CFL_ * H / myvmax);

  return mydt;
}

static vtype timestep_computation_vel_plus_accel_HEAD_ONLY(const struct node *ptr_node)
{
  vtype aux_dt;
  vtype mydt; // The time-step of the node

  vtype H = 1.0 / (1 << ptr_node->lv);

  vtype aux_const;
  vtype aux_epsilon;

  mydt = _MAX_dt_;
  for (int i = 0; i < GL_no_ptcl_final; i++)
  {
    //* >> x-axis
    if (myabs(GL_ptcl_ax[i]) < 1.0e-16)
    {

      if (myabs(GL_ptcl_vx[i]) > 1.0e-16)
      {
        aux_dt = _CFL_ * H / myabs(GL_ptcl_vx[i]);
        mydt = mydt < aux_dt ? mydt : aux_dt;
      }
    }
    else
    {
      aux_const = GL_ptcl_vx[i] * GL_ptcl_vx[i] + 2.0 * myabs(GL_ptcl_ax[i]) * H * _CFL_;
      if( 1 - GL_ptcl_vx[i] * GL_ptcl_vx[i] / aux_const  > 1.0e-8 )  //Expansion in Taylor to avoid the compution -b + sqrt(b^2 + c), with c << b
      {
        aux_dt = (-myabs(GL_ptcl_vx[i]) + mysqrt(aux_const)) / myabs(GL_ptcl_ax[i]);
      }
      else
      {
        aux_epsilon = myabs (GL_ptcl_ax[i] * H * _CFL_ / (GL_ptcl_vx[i] * GL_ptcl_vx[i]));
        aux_dt = myabs(GL_ptcl_vx[i]/GL_ptcl_ax[i]) * (aux_epsilon - aux_epsilon*aux_epsilon*0.5);
      }
      mydt = mydt < aux_dt ? mydt : aux_dt;
    }


    //* >> y-axis
    if (myabs(GL_ptcl_ay[i]) < 1.0e-16)
    {
      if (myabs(GL_ptcl_vy[i]) > 1.0e-16)
      {
        aux_dt = _CFL_ * H / myabs(GL_ptcl_vy[i]);
        mydt = mydt < aux_dt ? mydt : aux_dt;
      }
    }
    else
    {
      aux_const = GL_ptcl_vy[i] * GL_ptcl_vy[i] + 2.0 * myabs(GL_ptcl_ay[i]) * H * _CFL_;
      if( 1 - GL_ptcl_vy[i] * GL_ptcl_vy[i] / aux_const  > 1.0e-8 )  //Expansion in Taylor to avoid the compution -b + sqrt(b^2 + c), with c << b
      {
        aux_dt = (-myabs(GL_ptcl_vy[i]) + mysqrt(aux_const)) / myabs(GL_ptcl_ay[i]);
      }
      else
      {
        aux_epsilon = myabs (GL_ptcl_ay[i] * H * _CFL_ / (GL_ptcl_vy[i] * GL_ptcl_vy[i]));
        aux_dt = myabs(GL_ptcl_vy[i]/GL_ptcl_ay[i]) * (aux_epsilon - aux_epsilon*aux_epsilon*0.5);
      }
      mydt = mydt < aux_dt ? mydt : aux_dt;
    }


    //* >> z-axis
    if (myabs(GL_ptcl_az[i]) < 1.0e-16)
    {
      if (myabs(GL_ptcl_vz[i]) > 1.0e-16)
      {
        aux_dt = _CFL_ * H / myabs(GL_ptcl_vz[i]);
        mydt = mydt < aux_dt ? mydt : aux_dt;
      }
    }
    else
    {
      aux_const = GL_ptcl_vz[i] * GL_ptcl_vz[i] + 2.0 * myabs(GL_ptcl_az[i]) * H * _CFL_;
      if( 1 - GL_ptcl_vz[i] * GL_ptcl_vz[i] / aux_const  > 1.0e-8 )  //Expansion in Taylor to avoid the compution -b + sqrt(b^2 + c), with c << b
      {
        aux_dt = (-myabs(GL_ptcl_vz[i]) + mysqrt(aux_const)) / myabs(GL_ptcl_az[i]);
      }
      else
      {
        aux_epsilon = myabs (GL_ptcl_az[i] * H * _CFL_ / (GL_ptcl_vz[i] * GL_ptcl_vz[i]));
        aux_dt = myabs(GL_ptcl_vz[i]/GL_ptcl_az[i]) * (aux_epsilon - aux_epsilon*aux_epsilon*0.5);
      }
      mydt = mydt < aux_dt ? mydt : aux_dt;
    }



    // //* >> y-axis
    // if (myabs(GL_ptcl_ay[i]) <= 1.0e-12)
    // {
    //   if (myabs(GL_ptcl_vy[i]) > 1.0e-12)
    //   {
    //     aux_dt = _CFL_ * H / myabs(GL_ptcl_vy[i]);
    //     mydt = mydt < aux_dt ? mydt : aux_dt;
    //   }
    // }
    // else
    // {
    //   aux_dt = (-myabs(GL_ptcl_vy[i]) + mysqrt(GL_ptcl_vy[i] * GL_ptcl_vy[i] + 2.0 * myabs(GL_ptcl_ay[i]) * H * _CFL_)) / myabs(GL_ptcl_ay[i]);
    //   mydt = mydt < aux_dt ? mydt : aux_dt;
    // }

    // //* >> z-axis
    // if (myabs(GL_ptcl_az[i]) <= 1.0e-12)
    // {
    //   if (myabs(GL_ptcl_vz[i]) > 1.0e-12)
    //   {
    //     aux_dt = _CFL_ * H / myabs(GL_ptcl_vz[i]);
    //     mydt = mydt < aux_dt ? mydt : aux_dt;
    //   }
    // }
    // else
    // {
    //   aux_dt = (-myabs(GL_ptcl_vz[i]) + mysqrt(GL_ptcl_vz[i] * GL_ptcl_vz[i] + 2.0 * myabs(GL_ptcl_az[i]) * H * _CFL_)) / myabs(GL_ptcl_az[i]);
    //   mydt = mydt < aux_dt ? mydt : aux_dt;
    // }
  }

  return mydt;
}

static vtype timestep_computation_vel_only(const struct node *ptr_node)
{
  int lv; // Level of refinement

  vtype aux_v_pow2; // Auxiliary square velocity

  vtype H;

  int box_idx; // Box index of the node cell

  int ptcl_idx; // Particle grid_idx in the node

  vtype mydt;   // The time-step of the node
  vtype myvmax; // Maximum velocity of all particles

  lv = ptr_node->lv;
  H = 1.0 / (1 << lv);

  myvmax = 0.0; // Minium velocity designated by the user

  int counter_ptcl = 0;
  int total_ptcl = ptr_node->no_ptcl_outs_ref_zones;
  int cell_ptcl;
  int cell_idx = -1;

  //* >> Case no more child, the node is a leaf *//
  if (ptr_node->chn_size == 0)
  {
    while (counter_ptcl < total_ptcl)
    {
      cell_idx++;
      box_idx = ptr_node->ptr_box_idx[cell_idx];
      cell_ptcl = ptr_node->ptr_cell_struct[box_idx].ptcl_size;

      for (int j = 0; j < cell_ptcl; j++)
      {
        ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

        //* >> Velocity at x-axis
        aux_v_pow2 = myabs(GL_ptcl_vx[ptcl_idx]);
        if (myvmax < aux_v_pow2)
        {
          myvmax = aux_v_pow2;
        }
        //* >> Velocity at y-axis
        aux_v_pow2 = myabs(GL_ptcl_vy[ptcl_idx]);
        if (myvmax < aux_v_pow2)
        {
          myvmax = aux_v_pow2;
        }
        //* >> Velocity at z-axis
        aux_v_pow2 = myabs(GL_ptcl_vz[ptcl_idx]);
        if (myvmax < aux_v_pow2)
        {
          myvmax = aux_v_pow2;
        }
      }

      counter_ptcl += cell_ptcl;
    }
  }
  //* >> Case there are more children, the node is a branch *//
  else
  {
    while (counter_ptcl < total_ptcl)
    {
      cell_idx++;
      box_idx = ptr_node->ptr_box_idx[cell_idx];
      cell_ptcl = ptr_node->ptr_cell_struct[box_idx].ptcl_size;

      if (ptr_node->ptr_box[box_idx] == -3)
      {
        for (int j = 0; j < cell_ptcl; j++)
        {
          ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

          //* >> Velocity at x-axis
          aux_v_pow2 = myabs(GL_ptcl_vx[ptcl_idx]);
          if (myvmax < aux_v_pow2)
          {
            myvmax = aux_v_pow2;
          }
          //* >> Velocity at y-axis
          aux_v_pow2 = myabs(GL_ptcl_vy[ptcl_idx]);
          if (myvmax < aux_v_pow2)
          {
            myvmax = aux_v_pow2;
          }
          //* >> Velocity at z-axis
          aux_v_pow2 = myabs(GL_ptcl_vz[ptcl_idx]);
          if (myvmax < aux_v_pow2)
          {
            myvmax = aux_v_pow2;
          }
        }

        counter_ptcl += cell_ptcl;
      }
    }
  }

  // if (ptr_node->chn_size == 0)
  // {
  //     for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
  //     {
  //         box_idx = ptr_node->ptr_box_idx[cell_idx];

  //         for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
  //         {
  //             ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

  //             //* >> Velocity at x-axis
  //             aux_v_pow2 = myabs(GL_ptcl_vx[ptcl_idx]);
  //             if (myvmax < aux_v_pow2)
  //             {
  //                 myvmax = aux_v_pow2;
  //             }
  //             //* >> Velocity at y-axis
  //             aux_v_pow2 = myabs(GL_ptcl_vy[ptcl_idx]);
  //             if (myvmax < aux_v_pow2)
  //             {
  //                 myvmax = aux_v_pow2;
  //             }
  //             //* >> Velocity at z-axis
  //             aux_v_pow2 = myabs(GL_ptcl_vz[ptcl_idx]);
  //             if (myvmax < aux_v_pow2)
  //             {
  //                 myvmax = aux_v_pow2;
  //             }
  //         }
  //     }
  // }
  // //* >> Case there are more children, the node is a branch *//
  // else
  // {
  //     for (int cell_idx = 0; cell_idx < ptr_node->cell_size; cell_idx++)
  //     {
  //         box_idx = ptr_node->ptr_box_idx[cell_idx];

  //         if (ptr_node->ptr_box[box_idx] == -3)
  //         {
  //             for (int j = 0; j < ptr_node->ptr_cell_struct[box_idx].ptcl_size; j++)
  //             {
  //                 ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

  //                 //* >> Velocity at x-axis
  //                 aux_v_pow2 = myabs(GL_ptcl_vx[ptcl_idx]);
  //                 if (myvmax < aux_v_pow2)
  //                 {
  //                     myvmax = aux_v_pow2;
  //                 }
  //                 //* >> Velocity at y-axis
  //                 aux_v_pow2 = myabs(GL_ptcl_vy[ptcl_idx]) ;
  //                 if (myvmax < aux_v_pow2)
  //                 {
  //                     myvmax = aux_v_pow2;
  //                 }
  //                 //* >> Velocity at z-axis
  //                 aux_v_pow2 = myabs(GL_ptcl_vz[ptcl_idx]);
  //                 if (myvmax < aux_v_pow2)
  //                 {
  //                     myvmax = aux_v_pow2;
  //                 }
  //             }
  //         }
  //     }
  // }

  mydt = myvmax < (_CFL_ * H / _MAX_dt_) ? _MAX_dt_ : (_CFL_ * H / myvmax);

  return mydt;
}

static vtype timestep_computation_vel_plus_accel(const struct node *ptr_node)
{
  int box_idx;  // Box index of the node cell
  int ptcl_idx; // Particle grid_idx in the node

  vtype aux_dt;

  vtype mydt; // The time-step of the node


  vtype aux_const;
  vtype aux_epsilon;

  vtype H = 1.0 / (1 << ptr_node->lv);

  mydt = _MAX_dt_;

  int counter_ptcl = 0;
  int total_ptcl = ptr_node->no_ptcl_outs_ref_zones;
  int cell_ptcl;
  int cell_idx = -1;

  //* >> Case no more child, the node is a leaf *//
  if (ptr_node->chn_size == 0)
  {
    while (counter_ptcl < total_ptcl)
    {
      cell_idx++;
      box_idx = ptr_node->ptr_box_idx[cell_idx];
      cell_ptcl = ptr_node->ptr_cell_struct[box_idx].ptcl_size;

      for (int j = 0; j < cell_ptcl; j++)
      {
        ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

        // //* >> x-axis
        // if (myabs(GL_ptcl_ax[ptcl_idx]) <= 1.0e-12)
        // {
        //   if (myabs(GL_ptcl_vx[ptcl_idx]) > 1.0e-12)
        //   {
        //     aux_dt = _CFL_ * H / myabs(GL_ptcl_vx[ptcl_idx]);
        //     mydt = mydt < aux_dt ? mydt : aux_dt;
        //   }
        // }
        // else
        // {
        //   aux_dt = (-myabs(GL_ptcl_vx[ptcl_idx]) + mysqrt(GL_ptcl_vx[ptcl_idx] * GL_ptcl_vx[ptcl_idx] + 2.0 * myabs(GL_ptcl_ax[ptcl_idx]) * H * _CFL_)) / myabs(GL_ptcl_ax[ptcl_idx]);
        //   mydt = mydt < aux_dt ? mydt : aux_dt;
        // }

        //* >> x-axis
        if (myabs(GL_ptcl_ax[ptcl_idx]) < 1.0e-16)
        {
          if (myabs(GL_ptcl_vx[ptcl_idx]) > 1.0e-16)
          {
            aux_dt = _CFL_ * H / myabs(GL_ptcl_vx[ptcl_idx]);
            mydt = mydt < aux_dt ? mydt : aux_dt;
          }
        }
        else
        {
          aux_const = GL_ptcl_vx[ptcl_idx] * GL_ptcl_vx[ptcl_idx] + 2.0 * myabs(GL_ptcl_ax[ptcl_idx]) * H * _CFL_;
          if( 1 - GL_ptcl_vx[ptcl_idx] * GL_ptcl_vx[ptcl_idx] / aux_const  > 1.0e-8 )  //Expansion in Taylor to avoid the compution -b + sqrt(b^2 + c), with c << b
          {
            aux_dt = (-myabs(GL_ptcl_vx[ptcl_idx]) + mysqrt(aux_const)) / myabs(GL_ptcl_ax[ptcl_idx]);
          }
          else
          {
            aux_epsilon = myabs (GL_ptcl_ax[ptcl_idx] * H * _CFL_ / (GL_ptcl_vx[ptcl_idx] * GL_ptcl_vx[ptcl_idx]));
            aux_dt = myabs(GL_ptcl_vx[ptcl_idx]/GL_ptcl_ax[ptcl_idx]) * (aux_epsilon - aux_epsilon*aux_epsilon*0.5);
          }
          mydt = mydt < aux_dt ? mydt : aux_dt;
        }

        //* >> y-axis
        if (myabs(GL_ptcl_ay[ptcl_idx]) < 1.0e-16)
        {
          if (myabs(GL_ptcl_vy[ptcl_idx]) > 1.0e-16)
          {
            aux_dt = _CFL_ * H / myabs(GL_ptcl_vy[ptcl_idx]);
            mydt = mydt < aux_dt ? mydt : aux_dt;
          }
        }
        else
        {
          aux_const = GL_ptcl_vy[ptcl_idx] * GL_ptcl_vy[ptcl_idx] + 2.0 * myabs(GL_ptcl_ay[ptcl_idx]) * H * _CFL_;
          if( 1 - GL_ptcl_vy[ptcl_idx] * GL_ptcl_vy[ptcl_idx] / aux_const  > 1.0e-8 )  //Expansion in Taylor to avoid the compution -b + sqrt(b^2 + c), with c << b
          {
            aux_dt = (-myabs(GL_ptcl_vy[ptcl_idx]) + mysqrt(aux_const)) / myabs(GL_ptcl_ay[ptcl_idx]);
          }
          else
          {
            aux_epsilon = myabs (GL_ptcl_ay[ptcl_idx] * H * _CFL_ / (GL_ptcl_vy[ptcl_idx] * GL_ptcl_vy[ptcl_idx]));
            aux_dt = myabs(GL_ptcl_vy[ptcl_idx]/GL_ptcl_ay[ptcl_idx]) * (aux_epsilon - aux_epsilon*aux_epsilon*0.5);
          }
          mydt = mydt < aux_dt ? mydt : aux_dt;
        }

        //* >> z-axis
        if (myabs(GL_ptcl_az[ptcl_idx]) < 1.0e-16)
        {
          if (myabs(GL_ptcl_vz[ptcl_idx]) > 1.0e-16)
          {
            aux_dt = _CFL_ * H / myabs(GL_ptcl_vz[ptcl_idx]);
            mydt = mydt < aux_dt ? mydt : aux_dt;
          }
        }
        else
        {
          aux_const = GL_ptcl_vz[ptcl_idx] * GL_ptcl_vz[ptcl_idx] + 2.0 * myabs(GL_ptcl_az[ptcl_idx]) * H * _CFL_;
          if( 1 - GL_ptcl_vz[ptcl_idx] * GL_ptcl_vz[ptcl_idx] / aux_const  > 1.0e-8 )  //Expansion in Taylor to avoid the compution -b + sqrt(b^2 + c), with c << b
          {
            aux_dt = (-myabs(GL_ptcl_vz[ptcl_idx]) + mysqrt(aux_const)) / myabs(GL_ptcl_az[ptcl_idx]);
          }
          else
          {
            aux_epsilon = myabs (GL_ptcl_az[ptcl_idx] * H * _CFL_ / (GL_ptcl_vz[ptcl_idx] * GL_ptcl_vz[ptcl_idx]));
            aux_dt = myabs(GL_ptcl_vz[ptcl_idx]/GL_ptcl_az[ptcl_idx]) * (aux_epsilon - aux_epsilon*aux_epsilon*0.5);
          }
          mydt = mydt < aux_dt ? mydt : aux_dt;
        }
      }
      counter_ptcl += cell_ptcl;
    }
  }
  //* >> Case there are more children, the node is a branch *//
  else
  {
    while (counter_ptcl < total_ptcl)
    {
      cell_idx++;
      box_idx = ptr_node->ptr_box_idx[cell_idx];
      cell_ptcl = ptr_node->ptr_cell_struct[box_idx].ptcl_size;

      if (ptr_node->ptr_box[box_idx] == -3)
      {
        for (int j = 0; j < cell_ptcl; j++)
        {
          ptcl_idx = ptr_node->ptr_cell_struct[box_idx].ptr_ptcl[j];

          //* >> x-axis
          if (myabs(GL_ptcl_ax[ptcl_idx]) < 1.0e-16)
          {
            if (myabs(GL_ptcl_vx[ptcl_idx]) > 1.0e-16)
            {
              aux_dt = _CFL_ * H / myabs(GL_ptcl_vx[ptcl_idx]);
              mydt = mydt < aux_dt ? mydt : aux_dt;
            }
          }
          else
          {
            aux_const = GL_ptcl_vx[ptcl_idx] * GL_ptcl_vx[ptcl_idx] + 2.0 * myabs(GL_ptcl_ax[ptcl_idx]) * H * _CFL_;
            if( 1 - GL_ptcl_vx[ptcl_idx] * GL_ptcl_vx[ptcl_idx] / aux_const  > 1.0e-8 )  //Expansion in Taylor to avoid the compution -b + sqrt(b^2 + c), with c << b
            {
              aux_dt = (-myabs(GL_ptcl_vx[ptcl_idx]) + mysqrt(aux_const)) / myabs(GL_ptcl_ax[ptcl_idx]);
            }
            else
            {
              aux_epsilon = myabs (GL_ptcl_ax[ptcl_idx] * H * _CFL_ / (GL_ptcl_vx[ptcl_idx] * GL_ptcl_vx[ptcl_idx]));
              aux_dt = myabs(GL_ptcl_vx[ptcl_idx]/GL_ptcl_ax[ptcl_idx]) * (aux_epsilon - aux_epsilon*aux_epsilon*0.5);
            }
            mydt = mydt < aux_dt ? mydt : aux_dt;
          }

          //* >> y-axis
          if (myabs(GL_ptcl_ay[ptcl_idx]) < 1.0e-16)
         {
            if (myabs(GL_ptcl_vy[ptcl_idx]) > 1.0e-16)
            {
              aux_dt = _CFL_ * H / myabs(GL_ptcl_vy[ptcl_idx]);
              mydt = mydt < aux_dt ? mydt : aux_dt;
            }
         }
          else
         {
            aux_const = GL_ptcl_vy[ptcl_idx] * GL_ptcl_vy[ptcl_idx] + 2.0 * myabs(GL_ptcl_ay[ptcl_idx]) * H * _CFL_;
            if( 1 - GL_ptcl_vy[ptcl_idx] * GL_ptcl_vy[ptcl_idx] / aux_const  > 1.0e-8 )  //Expansion in Taylor to avoid the compution -b + sqrt(b^2 + c), with c << b
            {
              aux_dt = (-myabs(GL_ptcl_vy[ptcl_idx]) + mysqrt(aux_const)) / myabs(GL_ptcl_ay[ptcl_idx]);
            }
            else
            {
              aux_epsilon = myabs (GL_ptcl_ay[ptcl_idx] * H * _CFL_ / (GL_ptcl_vy[ptcl_idx] * GL_ptcl_vy[ptcl_idx]));
              aux_dt = myabs(GL_ptcl_vy[ptcl_idx]/GL_ptcl_ay[ptcl_idx]) * (aux_epsilon - aux_epsilon*aux_epsilon*0.5);
            }
            mydt = mydt < aux_dt ? mydt : aux_dt;
         }

          //* >> z-axis
          if (myabs(GL_ptcl_az[ptcl_idx]) < 1.0e-16)
          {
            if (myabs(GL_ptcl_vz[ptcl_idx]) > 1.0e-16)
            {
              aux_dt = _CFL_ * H / myabs(GL_ptcl_vz[ptcl_idx]);
              mydt = mydt < aux_dt ? mydt : aux_dt;
            }
          }
          else
          {
            aux_const = GL_ptcl_vz[ptcl_idx] * GL_ptcl_vz[ptcl_idx] + 2.0 * myabs(GL_ptcl_az[ptcl_idx]) * H * _CFL_;
            if( 1 - GL_ptcl_vz[ptcl_idx] * GL_ptcl_vz[ptcl_idx] / aux_const  > 1.0e-8 )  //Expansion in Taylor to avoid the compution -b + sqrt(b^2 + c), with c << b
            {
              aux_dt = (-myabs(GL_ptcl_vz[ptcl_idx]) + mysqrt(aux_const)) / myabs(GL_ptcl_az[ptcl_idx]);
            }
            else
            {
              aux_epsilon = myabs (GL_ptcl_az[ptcl_idx] * H * _CFL_ / (GL_ptcl_vz[ptcl_idx] * GL_ptcl_vz[ptcl_idx]));
              aux_dt = myabs(GL_ptcl_vz[ptcl_idx]/GL_ptcl_az[ptcl_idx]) * (aux_epsilon - aux_epsilon*aux_epsilon*0.5);
            }
            mydt = mydt < aux_dt ? mydt : aux_dt;
          }
        }
        counter_ptcl += cell_ptcl;
      }
    }
  }
  return mydt;
}

int timestep(vtype *ptr_dt)
{

  //* >> Time-step computing *//

  vtype dt_min; // minimum dt
  vtype aux_dt; // Auxiliary time

  dt_min = _MAX_dt_;

  if (lmin < lmax)
  {
    int total_iter = 0;
    for (int lv = GL_tentacles_level_max; lv > -1; lv--)
    {
      total_iter += GL_tentacles_size[lv];
    }

    int lv_idx[total_iter], node_idx[total_iter];  
    int idx = 0;
    for (int lv = GL_tentacles_level_max; lv > -1; lv--)
    {
      for (int i = 0; i < GL_tentacles_size[lv]; i++)
      {
        lv_idx[idx] = lv;
        node_idx[idx] = i;
        idx += 1;
      }
    }

    

    if (time_step_method == 1)
    {

      #pragma omp parallel for reduction(min:dt_min)
      for(int hh = 0; hh < total_iter; hh++ )
      {
        aux_dt = timestep_computation_vel_plus_accel(GL_tentacles[lv_idx[hh]][node_idx[hh]]);
        dt_min = dt_min > aux_dt ? aux_dt : dt_min;
      }

      // for (int lv = GL_tentacles_level_max; lv > -1; lv--)
      // {
      //   // Number of parent of the level = GL_tentacles_size[lv];

      //   //* >> For cycle over parent nodes *//
      //   for (int i = 0; i < GL_tentacles_size[lv]; i++)
      //   {
      //     // ptr_node = GL_tentacles[lv][i];

      //     //* >> Computing the time step of the node
      //     aux_dt = timestep_computation_vel_plus_accel(GL_tentacles[lv][i]);

      //     if (dt_min > aux_dt)
      //     {
      //       dt_min = aux_dt;
      //     }
      //   }
      // }
    }
    else // case time_step_method == 0
    {

      #pragma omp parallel for reduction(min:dt_min)
      for(int hh = 0; hh < total_iter; hh++ )
      {
        aux_dt = timestep_computation_vel_only(GL_tentacles[lv_idx[hh]][node_idx[hh]]);
        dt_min = dt_min > aux_dt ? aux_dt : dt_min;
      }

      // for (int lv = GL_tentacles_level_max; lv > -1; lv--)
      // {
      //   // NUmber of parents of the level = GL_tentacles_size[lv];

      //   //* >> For cycle over parent nodes *//
      //   for (int i = 0; i < GL_tentacles_size[lv]; i++)
      //   {
      //     // ptr_node = GL_tentacles[lv][i];

      //     //* >> Computing the time step of the node
      //     aux_dt = timestep_computation_vel_only(GL_tentacles[lv][i]);

      //     if (dt_min > aux_dt)
      //     {
      //       dt_min = aux_dt;
      //     }
      //   }
      // }
    }
  }
  else
  {
    if (time_step_method == 1)
    {
      aux_dt = timestep_computation_vel_plus_accel_HEAD_ONLY(GL_tentacles[0][0]);
    }
    else // case time_step_method == 0
    {
      aux_dt = timestep_computation_vel_only_HEAD_ONLY(GL_tentacles[0][0]);
    }

    if (dt_min > aux_dt)
    {
      dt_min = aux_dt;
    }
  }

  *ptr_dt = dt_min;

  return _SUCCESS_;
}