/**
 * Copyright (C) (2010-2025) Vadim Biktashev, Irina Biktasheva et al. 
 * (see ../AUTHORS for the full list of contributors)
 *
 * This file is part of Beatbox.
 *
 * Beatbox is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Beatbox is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Beatbox.  If not, see <http://www.gnu.org/licenses/>.
 */

/* List of all existing devices */

/*
LEGEND
------------------------------------------------------
D(device_name) will run sequentially and with MPI.
S(device_name) will run sequentially only.
*/

/**************** UNIVERSAL *************************/

D(activation)
D(byteout)
D(clock)
D(ctlpoint)
D(d_dt)
D(d_dx)
D(d_dy)
D(diff)
D(dump)
D(ecg)
D(elliptic)
D(euler)
D(flip)
D(globalfunc)
D(grad2d)
D(k_func)
D(k_print)
D(k_poincare)
D(l_abs)  /* u[v1]=abs(u[v0])    */
D(l_add)  /* u[v1]=u[v0]+a       */
D(l_axpy) /* u[v2]=a*u[v0]+u[v1] */
D(l_copy) /* u[v1]=u[v0]         */
D(l_max)  /* u[v2]=max(u[v0],u[v1]) */
D(l_min)  /* u[v2]=max(u[v0],u[v1]) */
D(l_move) /* move along the grid */
D(l_mult) /* u[v1]=u[v0]*a       */
D(l_prod) /* u[v2]=u[v0]*u[v1]   */
D(l_scal) /* u[v1]=u[v0]*a       */
D(l_sum)  /* u[v2]=u[v0]+u[v1]   */
D(l_swap) /* u[v1] <=> u[v0]     */
D(load)
D(localfunc)
D(mpipaint)
D(neum2d)
D(neum3d)
D(pause)
D(poincare)
D(pgmout)
D(ppmout)
D(pw_mult)
D(reduce)
D(regular)
D(rk4)
D(rushlarsen)
D(sample)
D(shell)
D(singz)
D(stop)
D(vtkout2)

/**************** SEQUENTIAL ONLY ******************/

S(adi1d)
S(adi2d)
S(adi3d)
S(bytein)
S(ezpaint)
S(ezstep)
S(ezview)
S(imgout)
S(inner_bc)
S(k_clock)
S(k_imgout)
S(k_draw)
S(k_paint)
S(k_paintgl)
S(k_plot)
S(matout)
S(pgmin)
S(ppmin)
S(record)
S(skrecord)
S(torx)
S(tory)
S(torz)
S(update)
S(screen_dump)

/******************** obsolete **********************/

D(diff2dv)
D(diff3dv)
D(diffstep)
D(diffold)
