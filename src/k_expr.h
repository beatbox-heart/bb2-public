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

/* version of k_comp with expressions not assignements */

#undef SEPARATORS
#define SEPARATORS ";"
#define BLANK " \n\t\r"
  BEGINBLOCK("list=",buf); {
  int icode;
  char *pcode;
  char *s1=strdup(buf);
  for(S->ncode=0;strtok(S->ncode?NULL:s1,SEPARATORS);S->ncode++);
  if (!S->ncode) EXPECTED_ERROR("no expressions in \"%s\"",buf);
  FREE(s1);

  if NOT(S->code=calloc(S->ncode,sizeof(pp_fn)))
    ABORT("not enough memory for code array of %d",S->ncode);

  for(icode=0;icode<S->ncode;icode++) {
    if NOT(pcode=strtok(icode?NULL:buf,SEPARATORS)) EXPECTED_ERROR("internal error");
    S->code[icode] = compile(pcode,loctb,t_real); CHK(pcode);
    MESSAGE3("\x01""\n %s=%s%c",pdata,pcode,SEPARATORS[0]);
  }
} ENDBLOCK;
#undef SEPARATORS
#undef BLANK


