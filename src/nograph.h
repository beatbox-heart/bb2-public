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

#define  NOGRAPH_DUMMY(device_name) \
  typedef void STR; \
  RUN_HEAD(device_name) {} RUN_TAIL(device_name) \
  DESTROY_HEAD(device_name) {} DESTROY_TAIL(device_name) \
  CREATE_HEAD(device_name) { MESSAGE("/* WARNING: this version of Beatbox is compiled without graphics, device %s is unavailable and is replaced with a dummy */",#device_name);} CREATE_TAIL(device_name,0)
