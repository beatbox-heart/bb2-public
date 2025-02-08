/**
 * Copyright (C) (2010-2016) Vadim Biktashev, Irina Biktasheva et al. 
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

/* TS 2016: Definitions related to Markov chain description of ionic channels */

/* maximal number of allowed subchains in a Markov chain models */
#define MAX_SUBCHAINS 4

#define SUBCHAIN(fun_tr, index, min, max, incr, sc)	\
  sbch->trans_rates_mat = fun_tr;		\
  sbch->i_control = index;			\
  sbch->tmin = min;				\
  sbch->tmax = max;				\
  if (min > max) {					\
  URGENT_MESSAGE("min > max for subchain tabulation of %s", #fun_tr);	\
  ABORT("");								\
  }						\
  sbch->tincr = incr;				\
  sbch->scale = sc;				\
  ch->num_sub += 1;				\
  sbch+=1;

#define TR_MAT(channel,from, to, direct, reverse)			\
  tr_mat[channel##_##to*NM_##channel+channel##_##from]=direct;		\
  tr_mat[channel##_##from*NM_##channel+channel##_##from]-=direct;	\
  tr_mat[channel##_##from*NM_##channel+channel##_##to]=reverse;		\
  tr_mat[channel##_##to*NM_##channel+channel##_##to]-=reverse;

#define CHANNEL_TR_MATRIX(name) int name(real * u, real *tr_mat)
typedef CHANNEL_TR_MATRIX(TransRatesMat);

/* subchannel of Markov chain model */
typedef struct {
  TransRatesMat *trans_rates_mat; /* function to get transition rates matrix */
  int i_control;	  /* index of controling dynamical variable */
  real tmin;		  /* minimal value for tabulation */
  real tmax;		  /* maximal value for tabulation */
  real tincr;		  /* increment in the tabulation */
  int scale;		  /* 0 for linar, 1 for logarithmic (must be reflected in TransRatesMat) */
} subchain_str;

/* general structure of ionic channel to embarace Markov chain and gate model */
typedef struct {
  int dimension;		/* dimensionality of the model */
  int num_sub;			/* number of subchains */
  subchain_str * subchain;	/* subchains of the model */
  /* TODO: add conservation variable */
} channel_str;


