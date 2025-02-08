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


#ifndef SCR_H
#define SCR_H

/*extern void Logo(FILE **fil, Char *name);*/
/*extern void Textline2(Char *txt1, Char *txt2, short n, short l);*/
/*extern void SetHead(Char *F, Char *S);*/
/*extern void FileOpen(FILE **fil, Char *name, boolean b);*/
#define freset(f,name) if(NULL!=(f=fopen(name,"rt")))assign(f,name);else exit(printf("cannot reset %s\n",name))
#define frewrite(f,name) if(NULL!=(f=fopen(name,"wt")))assign(f,name);else exit(printf("cannot rewrite %s\n",name))
/*extern void ReadFileText(FILE **filin, FILE **filout, Char *name);*/
/*extern Char SReadKey(void);*/

#endif /*SCR_H*/
