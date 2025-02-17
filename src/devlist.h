/**
 * Copyright (C) (2010-2021) Vadim Biktashev, Irina Biktasheva et al. 
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
D(diff)
D(dump)
D(ecg)
D(elliptic)
D(euler)
D(grad2d)
D(k_func)
D(k_print)
D(k_poincare)
D(load)
D(neum2d)
D(neum3d)
D(pause)
D(poincare)
D(ppmout)
D(pw_mult)
D(reduce)
D(rk4)
D(rushlarsen)
D(sample)
D(shell)
D(singz)
D(stop)
D(vtkout2)

/**************** SEQUENTIAL ONLY ******************/

S(adi3d)
S(bytein)
S(ezpaint)
S(ezstep)
S(ezview)
S(imgout)
S(k_clock)
S(k_imgout)
S(k_draw)
S(k_paint)
S(k_paintgl)
S(k_plot)
S(matout)
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
