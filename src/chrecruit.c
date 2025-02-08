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

/*
 * Weight correction (>=0) 
 * and recruitment as in Cushing-Horwood 1994 
 * Finite operation (NOT evolution in time)
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"
#include "beatbox.h"

#include "state.h"
#include "device.h"
#include "bikt.h"
#include "rhs.h"
#include "qpp.h"

				/* indices of dynamic variables */
enum {
  #define _(n,v,d,unit) var_##n, 
  #include "chrecruitvar.h"
  #undef _
  Nvar /* total number of variables */
};

typedef struct {
				/* parameters */
  #define _(n,v,d,unit) real n;
  #include "chrecruitpar.h"
  #undef _
				/* external currents */
  #define _(n,v,d,unit) real I##n;
  #include "chrecruitvar.h"
  #undef _
} rhspar;

RHS_HEAD(chrecruit,Nvar) {
  #define _(n,v,d,unit) DEVICE_CONST(real,n)
  #include "chrecruitpar.h"
  #undef _
  #define _(n,v,d,unit) real n = u[var_##n];
  #include "chrecruitvar.h"
  #undef _

  W=max(W,Wmin);
  if (W<=0) {
    W=0; 
    N=0;
  } else if (W>metwgt || a>metage) {
    FN+=N;
    FB+=N*W;
    W=0;
    N=0;    
  } 

  #define _(n,v,d,unit) u[var_##n]=n;
  #include "chrecruitvar.h"
  #undef _
  #define _(n,v,d,unit) du[var_##n]=0;
  #include "chrecruitvar.h"
  #undef _
} RHS_TAIL(chrecruit)

RHS_CREATE_HEAD(chrecruit) {
  int i;
  #define _(n,v,d,unit) ACCEPTP(n,v,RNONE,RNONE);
  #include "chrecruitpar.h"
  #undef _

  #define _(n,v,d,unit) ACCEPTP(I##n,0,RNONE,RNONE);
  #include "chrecruitvar.h"
  #undef _
  
  MALLOC(*u,Nvar*sizeof(real));
  for(i=0;i<Nvar;i++) {
    #define _(n,v,d,unit) (*u)[var_##n]=v;
    #include "chrecruitvar.h"
    #undef _
  }
} RHS_CREATE_TAIL(chrecruit,Nvar)
