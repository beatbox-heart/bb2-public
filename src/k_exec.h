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

/* Execute a k-program previously compiled by k_comp.h snippet */
{
  int icode;
  p_vb *result;
  pp_fn the_code;
  for (icode=0;icode<ncode;icode++) {
    the_code = code[icode];
    if (!the_code) break;
    result = execute(the_code);
    CHK(NULL);
    memcpy(data[icode],result,sizetable[res_type(the_code)]);
  }
}


