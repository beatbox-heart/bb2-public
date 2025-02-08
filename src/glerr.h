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

static char *glerr(GLenum err)
{
  switch (err) {
  case GL_NO_ERROR: return "no error";
  case GL_INVALID_ENUM: return "invalid enum";
  case GL_INVALID_VALUE: return "invalid value";
  case GL_INVALID_OPERATION: return "invalid operation";
  case GL_STACK_OVERFLOW: return "stack overflow";
  case GL_STACK_UNDERFLOW: return "stack underflow";
  case GL_OUT_OF_MEMORY: return "out of memory";
  case GL_TABLE_TOO_LARGE: return "table too large";
  default: return "unknown error";
  }      
}
  
