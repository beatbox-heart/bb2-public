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


#ifndef INITIAL_H
#define INITIAL_H

#include "p2c.h"
#include "extern.h"

#define depth           10   /* depth is no. of diffusion layers */
#define Cmdlen          80

#define Infinite        100000L   /* Infinite is beyond likely Tend */

#define Faraday         965.0e2   /* Faraday constant */
#define delta3          0.0001   /* small value */
#define delta4          0.00001   /* smaller value */
#define delta6          0.0000001   /* very small value */

#define Version         4.0   /* Program version */

#define Update          "1st March 1993"



typedef char Vsub;

typedef double Ardsub[12];
typedef double Ar10[10];
typedef char CmdIndex;
typedef char Cmdbuffer[Cmdlen];
typedef char pstring[13];
typedef char nstring[21];
typedef char lstring[81];
typedef char mstring[57];
typedef char llstring[65];
typedef char sstring[6];

EXTERN nstring File0, File1, File2, File3, File4, File5; /* names of files */
EXTERN FILE *Filin2,*Filin3,*Filin4,*Filout1,*Filout2,*Filout4;

/* DEVICE_CONSTANT AND COMMON PARAMETERS */
#define FILEP FILE *
#define __AA(type,name) EXTERN type name;
#include "noblecns.h"
#include "nobleglb.h"
#undef __AA

/* screen displayed information on computation being run */
EXTERN char model[256];

/* interval between restart files */
EXTERN double restrep;


/* first restart file number */
EXTERN short restfileno;



EXTERN unsigned short hour, minute, second, sec100;


EXTERN char licence[256], site[256], contact[256];
EXTERN char Filin1_NAME[_FNSIZE];
EXTERN char Filin2_NAME[_FNSIZE];
EXTERN char Filin3_NAME[_FNSIZE];
EXTERN char Filin4_NAME[_FNSIZE];
EXTERN char Filout1_NAME[_FNSIZE];
EXTERN char Filout2_NAME[_FNSIZE];
EXTERN char Filout3_NAME[_FNSIZE];
EXTERN char Filout4_NAME[_FNSIZE];

/* Define power function PoN */

/*extern double PoN(double ARG, double POW);*/
#define PoN(ARG,POW) pow((ARG),(POW))

/* standard Pascal does not have a power function. This defines one */



/* Initiate all global variables described above
(? - but for some obviously safe) */
extern void InitiateAll(void);

#undef EXTERN

#endif /*INITIAL_H*/
