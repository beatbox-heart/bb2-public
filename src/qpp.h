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

/*====================================================================*/
/* Exported functions prototypes                                      */
/*====================================================================*/

/*! 
 * Initiation of constants used in checking limits etc.
 * 	\return 1 for success, 0 for failure.
 */
int init_const (void);

/*! Clear the symbol table and deallocate variables.
 * 	\return 1 for success, 0 for failure.
 */
int term_const (void);

/*! Opens for reading an included <file> or a string [macro].
 * 	\param name The name of the [macro] or the <file>
 * 	\return 1 for success, 0 for failure.
 */
int qopen (char *name);

/*! Close <file> or [macro] previously opened by qopen.
 *	\return 1 for success, 0 for failure.
*/
int qclose (void);

/*! Points w to first word in s, and return pointer to next after
 *	\param s String to search.
 *	\param w Pointer to the first word.
 *      \param delim String containing possible delimiter characters.
 *	\return Pointer to the rest of the search string.
 */
char *first_word (char *s, char **w,char *delim);

/*! Read to next terminator or end of file.
 *      \param s string where the result will be placed.
 *      \param strlen length of \a s.
 *	\return Number of bytes read, or 0 if end of file (error) .
*/
int read_command (char *s, int strlen);

/*! Evaluates a k-expression.
 *	\param resaddr The address into which the result will be written.
 *	\param restype The datatype of the result.
 *	\param expr The expression to be evaluated.
 * 	\return 1 for success, 0 for failure.
 */
int calc (void *resaddr, int restype, char *expr);

/*! Defines (and possibly initiates) a global k-variable in the default user table (deftb).
 *	\param s The string to be parsed.
 * 	\return 1 for success, 0 for failure.
 */
int def (char *s);

/*! Defines (and possibly initiates) a local k-variable in the given symbol table.
 *	\param s The string to be parsed.
 *	\param table The symbol table into which the variable is to be inserted.
 * 	\return 1 for success, 0 for failure.
 */
int def_local (char *s, p_tb table);

/*! Define a device name as a global strig macro. 
 *	\param d The device whose name is to be defined.
 * 	\return 1 for success, 0 for failure.
 */
int def_dev (Device *d);

/*! Obtain a reference to a device from global string macro containing its name.
 *	\param n The name of the device to be found.
 *	\param d Address of the pointer to the named device, if found.
 * 	\return 1 for success, 0 for failure.
 */
int get_dev (Name n, Device **d);

/*! Find a key in buffer, return to the next char after key or NULL if not 
 *  \param key the string to be found, typically like "name="
 *  \param s string to be searched for thekey
 *  \return Pointer to the next character after the key or NULL if none.
 */
char *find_key (const char *key, char *s);

/*! Accepts an integer value, identified by name, from a given 
 *      parameter string.
 *	\param name Name of the parameter whose value is to be read.
 *	\param var Address of variable where value is to be assigned.
 *	\param deflt Default value, in case no parameter called "name" 
 *                   is found. If INONE, a script parameter MUST be found.
 *	\param minv Minimum allowed value. If INONE, no minimum is prescribed.
 *	\param maxv Maximum allowed value. If INONE, no maximum is prescribed.
 *	\param w Parameter string to search.
 * 	\return 1 for success, 0 for failure.
 */
int accepti (const char *name,int *var,int deflt, int minv, int maxv, char *w);

/*! Accepts an INT value, identified by name, from a given 
 *      parameter string. The value may be linked to k-variable or k-expression. 
 *	\param name Name of the parameter whose value is to be read.
 *	\param varptr pointer to the address of the accepted value, if fixed
 *	\param varval address of k-variable containing the current value
 *	\param varcode address of k-code producing the current value
 *	\param deflt Default value, in case no parameter called "name" 
 *                   is found. If INONE, a script parameter MUST be found.
 *	\param minv Minimum allowed value. If INONE, no minimum is prescribed.
 *	\param maxv Maximum allowed value. If INONE, no maximum is prescribed.
 *	\param w Parameter string to search.
 * 	\return 1 for success, 0 for failure.
 */
int acceptik (const char *name, 
	      INT **varptr, INT *varval, char **varsrc, pp_fn *varexe, p_tb vartab,
	      INT deflt, INT minv, INT maxv, char *w);

/*! Accepts an entry of an integer parameter array.
 *	\param mask for the parameter name e.g. a%d for a0, a1, etc.
 *	\param arr the integer array whose value is to be accepted. 
 *	\param i the enumerator of the array so i=1 gives a1 in beatbox script and arr[1] in C code.
 *	\param deflt Default value for the identified array element; mandatory if INONE.
 *	\param minv Minimum allowed value. If INONE, no minimum is prescribed.
 *	\param maxv Maximum allowed value. If INONE, no maximum is prescribed.
 *	\param w Parameter string to search.
 * 	\return 1 for success, 0 for failure.
 */
int acceptie (const char *mask,int *arr,int i,int deflt,int minv,int maxv,char *w);

/*! Accepts a long value, identified by name, from a given parameter string.
 *	\param name Name of the parameter whose value is to be read.
 *	\param var Address of variable where value is to be assigned.
 *	\param deflt Default value, in case no parameter called "name" 
 *                   is found. If LNONE, a script parameter MUST be found.
 *	\param minv Minimum allowed value. If LNONE, no minimum is prescribed.
 *	\param maxv Maximum allowed value. If LNONE, no maximum is prescribed.
 *	\param w Parameter string to search.
 * 	\return 1 for success, 0 for failure.
 */
int acceptl (const char *name,long *var,long deflt,long minv,long maxv,char *w);

/*!
 *	Accepts a real (double) value, identified by name, from a given 
 *      parameter string.
 *	\param name Name of the parameter whose value is to be read.
 *	\param var Address of variable where value is to be assigned.
 *	\param deflt Default value, in case no parameter called "name" is 
 *                   found. If RNONE, a script parameter MUST be found.
 *	\param minv Minimum allowed value. If RNONE, no minimum is prescribed.
 *	\param maxv Maximum allowed value. If RNONE, no maximum is prescribed.
 *	\param w Parameter string to search.
 *	\return Flag indicating success. 1 if successful, 0 otherwise.
 */
int acceptr (const char *name,real *var,real deflt,real minv,real maxv,char *w);

/*! Accepts a REAL value, identified by name, from a given parameter string. 
 *	The value may be linked to k-variable or k-expression. 
 *	\param name Name of the parameter whose value is to be read (in)
 *	\param varptr pointer to the address of the accepted value, if fixed (out)
 *	\param varval address (in) of k-variable containing the current value (out)
 *	\param varcode address of k-code producing the current value (out)
 *      \param vartab address of k-symtable for the k-code (in)
 *	\param deflt Default value, in case no parameter "name" is found.
 *	       If INONE, the parameter is mandatory (in)
 *	\param minv Minimum allowed value. If INONE, no minimum is prescribed (in)
 *	\param maxv Maximum allowed value. If INONE, no maximum is prescribed (in)
 *	\param w Parameter string to search (in)
 * 	\return 1 for success, 0 for failure
 */
int acceptrk (const char *name, 
	      REAL **varptr, REAL *varval, char **varsrc, pp_fn *varexe, p_tb vartab,
	      real deflt, real minv, real maxv, char *w);

/*! Accepts an entry of a real parameter array.
 *	\param mask for the parameter name e.g. a%d for a0, a1, etc.
 *	\param arr the integer array whose value is to be accepted. 
 *	\param i the enumerator of the array so i=1 gives a1 in beatbox script and arr[1] in C code.
 *	\param deflt Default value for the identified array element; mandatory if INONE.
 *	\param minv Minimum allowed value. If INONE, no minimum is prescribed.
 *	\param maxv Maximum allowed value. If INONE, no maximum is prescribed.
 *	\param w Parameter string to search.
 * 	\return 1 for success, 0 for failure.
 */
int acceptre (const char *mask,real *arr,int i,real deflt,real minv,real maxv,char *w);

#if defined(_rhs) || defined(_ionic)
/*! RHS ONLY: accepts a parameter or dependent parameter from the script.
 * 	NB functionally similar to acceptkr, but faster. 
 *	\param name  Name of the parameter to be found in the script.
 *	\param var   Address of variable where value is to be assigned.
 *	\param deflt Default value, in case no parameter called "name" 
 *                  is found. If RNONE, a script parameter MUST be found.
 *	\param minv Minimum allowed value. If RNONE, no minimum is prescribed.
 *	\param maxv Maximum allowed value. If RNONE, no maximum is prescribed.
 *	\param w Parameter string to search.
 *	\param v Address of the dependent-parameter-describing Var structure.
 *	\param iv Address of the variable indicating the index in v of the 
 *                current dependent parameter.
 *	\param v0 Lower variable bound for the calling device's Space 
 *                structure. Used to ensure that the v's src fields point to 
 *                variable 0.
 * 	\return 1 for success, 0 for failure.
 */
int acceptp (const char *name,real *var,real deflt,real minv,real maxv,char *w,Var *v,int *iv,int v0);
#endif

/*! Accepts a string, identified by name, from a given parameter string.
 *	\param name Name of the parameter whose value is to be read.
 *	\param var Address of variable where value is to be assigned.
 *	\param deflt Default value, in case no parameter called "name" 
 *             is found. If NULL, a script parameter MUST be found.
 *	\param w Parameter string to search.
 * 	\return 1 for success, 0 for failure.
 */
int accepts (const char *name,char *var,const char *deflt, char *w);

/*! Accepts a string of known max length, identified by name, from a given parameter string.
 *	\param name Name of the parameter whose value is to be read.
 *	\param var Address of variable where value is to be assigned.
 *      \param len maximal allowed length of the string.
 *	\param deflt Default value, in case no parameter called "name" 
 *             is found. If NULL, a script parameter MUST be found.
 *	\param w Parameter string to search.
 * 	\return 1 for success, 0 for failure.
 */
int acceptsn (const char *name,char *var,int len,const char *deflt, char *w);

/*! Accepts a file handle, identified by parameter name, from a given 
 *      parameter string.
 *	\param name Name of the parameter whose value is to be read.
 *	\param mode Access mode in which the file is to be opened. Syntax is 
 *                  identical to a fopen() call:
 *		- r or rb Open existing ﬁle for reading. 
 *		- w or wb Create ﬁle or wipe existing ﬁle before writing. 
 *		- a or ab Append to end of existing ﬁle, creating if necessary. 
 *		- rt or rbt or rtb Open existing ﬁle for update—reading and 
 *                writing. 
 *		- wt or wbt or wtb Create ﬁle or wipe existing ﬁle before 
 *                updating. 
 *		- at or abt or atb Append—Open or create ﬁle for update, 
 *                writing at end of ﬁle. 
 *	\param deflt Default filename, in case no parameter called "name" is 
 *             found in the input script. If NULL, a filename MUST be found in the script.
 *	\param fname Address at which to store the filename.
 *	\param f Pointer to the open file.
 *	\param w Parameter string to search.
 * 	\return 1 for success, 0 for failure.
 */
int acceptf (const char *name,const char *mode,const char *deflt,char *fname,FILE **f, char *w);

/*! Accepts a code-string; a k-expression in input script that will be evaluated at run time. 
 *	\param name Name of the parameter whose value is to be read.
 *	\param type Data type returned by the code-string when executed. 
 *	\param type Address of variable where value is to be assigned.
 *	\param deflt Default value, in case no parameter called "name" is 
 *             found. If NULL, a script parameter MUST be found.
 *	\param codestring Address at which the executable code should be stored.
 *	\param code Pointer to the compiled code block.
 *	\param w Parameter string to search.
 * 	\return 1 for success, 0 for failure.
 */
int acceptc (const char *name,const k_type type,const char *deflt,char *codestring,void **code, char *w);

/*! Accepts a REAL k-variable from a given parameter string.
 *	\param name Name of the parameter whose value will be the k-variable. 
 *	\param writable if 1, the k-variable must be writable (not read-only). 
 *	\param deflt Default value  (NULL: mandatory, "": omit the variable).
 *	\param varaddr address of the k-variable.
 *	\param varname name of the k-variable.
 *	\param w Parameter string to search.
 * 	\return 1 for success, 0 for failure.
 */
int acceptkr (const char *key,const int writable,const char *deflt,REAL **varaddr,char *varname,char *w);

/*! Accepts an INT k-variable from a given parameter string.
 *	\param name Name of the parameter whose value will be the k-variable. 
 *	\param writable if 1, the variable must be writable (not read-only). 
 *	\param deflt Default value  (NULL: mandatory, "": omit the variable).
 *	\param varaddr address of the k-variable.
 *	\param varname name of the k-variable.
 *	\param w Parameter string to search.
 * 	\return 1 for success, 0 for failure.
 */
int acceptki (const char *key,const int writable,const char *deflt,INT **varaddr,char *varname,char *w);

/*! "Opens" a {block} in Beatbox script.
 *	\param name Name of the block in the script, e.g. name={BLOCK}.
 *	\param var Address of the string where the block content is copied. 
 *	\param w Parameter string to search.
 * 	\return 1 for success, 0 for failure.
 */
int acceptb (const char *name,char *var, char *w);

/*! "Closes" a block opened by acceptb. 
 *	This is purely syntactic, i.e. does not do anything. 
 */
void close_block (void);

/*! Check for state dimensions (xmax,ymax,zmax) in the given string.
 *	\param parameterString The string to search.
 * 	\return 1 for yes, 0 for no.
 */
int stateDimensionsExist (char *parameterString);

/*! Check for any of the space parameters in the given string (check out Space structure). 
 *	\param parameterString The string to search.
 * 	\return 1 for yes, 0 for no.
 */
int spaceParametersExist (char *parameterString);

/*! Check for any of the space parameters that are beyond a straight box shape.
 *	\param parameterString The string to search.
 * 	\return 1 for yes, 0 for no.
 */
int spaceIsBox (char *parameterString);

/*! Accepts a device's Space constraints from a given string.
 *	\param s Address of Space structure to be initialised.
 *	\param w String to search.
 * 	\return 1 for success, 0 for failure.
 */
int accept_space (Space *s, char *w);

/*! (Re)make the list of points according to the given Space description.
 *	\param s Address of Space structure to be modified.
 * 	\return 1 for success, 0 for failure.
 */
int list_space (Space *s);

#ifdef X11
/*! Accepts a device's BGIWindow description from a given string.
 *	\param s Address of BGIWindow structure to be initialised.
 *	\param w String to search.
 * 	\return 1 for success, 0 for failure.
 */
int accept_window (BGIWindow *s, char *w);
#endif

/*! 
 * The smallest real bigger than val. 
 */
real RSUCC (real val);
/*! 
 * The biggest real smaller than val. 
 */
real RPRED (real val);

/*! Interpolated value of the field in given layer.
 * TODO: does it really belong here? 
 *	\param _x The x coordinate of the point.  
 *	\param _y The y coordinate of the point.  
 *	\param _z The z coordinate of the point.  
 *	\param _v Rounded down to integer, is the layer.
 *	\return The interpolated value: trilinear within grid, const outside.
 */
double _u (double _x, double _y, double _z, double _v);

/*====================================================================*/
/* Macros-shortcuts for the exported functions                        */
/*====================================================================*/

/*! Accept an integer parameter via accepti(), and define an eponymous local variable. 
 *	\param b Parameter name/key to be found.
 *	\param c Default value.
 *	\param d Minimum value.
 *	\param e Maximum value.
 */
#define ACCEPTI(b,c,d,e) if (!accepti(#b"=",&(S->b),c,d,e,w)) return(0); int b=S->b;

/*! Accepts an INT value possibly linked to k-variable or k-expression via acceptik(), with deftb (global vars only).
 * 	No local variable defined: unfeasible. 
 *	\param b Parameter name/key to be found.
 *	\param c Default value.
 *	\param d Minimum value.
 *	\param e Maximum value.
 */
#define ACCEPTIK(b,c,d,e) if (!acceptik(#b"=",&(S->b##ptr),&(S->b),&(S->b##src),&(S->b##exe),deftb,c,d,e,w)) return(0); 

/*! Accept an entry of an integer parameter array via acceptie(). 
 * 	No local variable defined: unfeasible. 
 *	\param m Mask for the parameter name e.g. arr%d for arr0, arr1, etc.
 *	\param a The integer array whose element's value is to be accepted. 
 *	\param i The enumerator of the array so i=1 gives arr1 in beatbox script and a[1] in C code.
 *	\param c Default value.
 *	\param d Minimum value.
 *	\param e Maximum value.
 */
#define ACCEPTIE(m,a,i,c,d,e) if (!acceptie(m"=",a,i,c,d,e,w)) return(0)

/*! Accept a long integer parameter via acceptl(), and define an eponymous local variable. 
 *	\param b Parameter name/key to be found.
 *	\param c Default value.
 *	\param d Minimum value.
 *	\param e Maximum value.
 */
#define ACCEPTL(b,c,d,e) if (!acceptl(#b"=",&(S->b),c,d,e,w)) return(0); long b=S->b;

/*! Accept a real parameter via acceptr(), and define an eponymous local variable. 
 *	\param b Parameter name/key to be found.
 *	\param c Default value.
 *	\param d Minimum value.
 *	\param e Maximum value.
 *	\sa acceptr()
 */
#define ACCEPTR(b,c,d,e) if (!acceptr(#b"=",&(S->b),c,d,e,w)) return(0); real b=S->b;

/*! Accepts a REAL value possibly linked to k-variable or k-expression via acceptrk(), with deftb (global vars only). 
 * 	No local variable defined: unfeasible. 
 *	\param b Parameter name/key to be found.
 *	\param c Default value.
 *	\param d Minimum value.
 *	\param e Maximum value.
 */
#define ACCEPTRK(b,c,d,e)  if (!acceptrk(#b"=",&(S->b##ptr),&(S->b),&(S->b##src),&(S->b##exe),deftb,c,d,e,w)) return(0);
/* Same, with loctb (including local vars x,y,z,u[] and geom[] */
#define ACCEPTRL(b,c,d,e)  if (!acceptrk(#b"=",&(S->b##ptr),&(S->b),&(S->b##src),&(S->b##exe),S->loctb,c,d,e,w)) return(0);
/* Same, using alias 'a' */
#define ACCEPTRA(a,b,c,d,e) if(!acceptrk(#a"=",&(S->b##ptr),&(S->b),&(S->b##src),&(S->b##exe),S->loctb,c,d,e,w)) return(0);
/* NB macros to work with these recalculable parameters are defined in device.h */
 
/*! Accept an entry of a real parameter array via acceptre(). 
 * 	No local variable defined: unfeasible. 
 *	\param m Mask for the parameter name e.g. arr%d for arr0, arr1, etc.
 *	\param a The integer array whose element's value is to be accepted. 
 *	\param i The enumerator of the array so i=1 gives arr1 in beatbox script and a[1] in C code.
 *	\param c Default value.
 *	\param d Minimum value.
 *	\param e Maximum value.
 */
#define ACCEPTRE(m,a,i,c,d,e) if (!acceptre(m"=",a,i,c,d,e,w)) return(0);

#if defined _rhs || defined _ionic
/*! RHS ONLY: accept a parameter or dependent parameter via acceptp(), and define an eponymous local variable. 
 *	\param b Parameter name/key to be found.
 *	\param c Default value.
 *	\param d Minimum value.
 *	\param e Maximum value.
 */
#define ACCEPTP(b,c,d,e) if (!acceptp(#b"=",&(S->b),c,d,e,w,var,&ivar,v0)) return(0); real b=S->b;
#endif

/*! Accept a string via accepts(), and define an eponymous local variable.
 *	\param b Parameter name/key to be found.
 *	\param c Default value.
 */
#define ACCEPTS(b,c)     if (!accepts(#b"=",&(S->b[0]),c,w)) return(0); char *b=&(S->b[0]);

/*! Accept a string of limited length via acceptsn(), and define an eponymous local variable.
 *	\param b Parameter name/key to be found.
 *	\param c Default value.
 */
#define ACCEPTSN(b,n,c)     if (!acceptsn(#b"=",&(S->b[0]),n,c,w)) return(0); char *b=&(S->b[0]);

/*! Accept a file via acceptf(), and define eponymous local variables.
 *	\param b Parameter name/key to be found.
 *	\param c Access mode.
 *	\param d Default value.
 *	\sa acceptf()
 */
/* #define ACCEPTF(b,c,d)   if (!acceptf(#b"=",c,d,&(S->b##name[0]),&(S->b),w)) {\ */
/*   DEBUG("not acceptf\n"); return(0); \ */
/*   }; DEBUG("yes acceptf, name='%s'\n",S->b##name); FILE *b=S->b; char *b##name=&(S->b##name[0]); */
#define ACCEPTF(b,c,d)   if (!acceptf(#b"=",c,d,&(S->b##name[0]),&(S->b),w)) return(0); \
  FILE *b=S->b; char *b##name=&(S->b##name[0]);


/*! Shortcut for entry in a device structure, suitable for accepting by ACCEPTF */
#define NFILE(a) FILE *a; char a##name[MAXPATH]

/*! Accept a k-expression code string via acceptc().
 * 	No local variable defined.
 *	\param b Parameter name/key to be found.
 *	\param c k_type to be returned by the code string.
 *	\param d Default code string.
 */
#define ACCEPTC(b,c,d)   if (!acceptc(#b"=",c,d,&(S->b[0]),&(S->b##code),w)) return(0);

/*! Accept a REAL k_variable via acceptkr(). 
 * 	No local variable defined.
 *	\param b Variable name/key to be found.
 * 	\param c The variable has to be writable.
 *	\param d The name of the default variable. 
 */
#define ACCEPTKR(b,c,d)   if (!acceptkr(#b"=",c,d,&(S->b),&(S->b##name[0]),w)) return(0);
/*! Shortcut for k-variable entry in a device structure, suitable for accepting by ACCEPTKR */
/*! Shortcut for k-variable entry in a device structure, suitable for accepting by ACCEPTKR */
#define KREAL(a) REAL *a; char a##name[maxname]

/*! Accept a INT k_variable via acceptki(). 
 * 	No local variable defined.
 *	\param b Variable name/key to be found.
 * 	\param c The variable has to be writable.
 *	\param d The name of the default variable. 
 */
#define ACCEPTKI(b,c,d)   if (!acceptki(#b"=",c,d,&(S->b[0]),&(S->b##name),w)) return(0);

/*! Accept a code block via acceptb().
 *	\param a Parameter name/key to be found.
 *	\param b Address to which code block should be assigned.
 */
#define BEGINBLOCK(a,b)    if (!acceptb (a,b,w)) return(0)

/*! Synonymous with close_block(). 
 *  Syntactic purposes only, no action. 
 */
#define ENDBLOCK	   close_block()


#ifdef X11
/*! Accept a device's BGIWindow description via accept_window(), and define eponymous local variable. 
 *	\param b Address of BGIWindow structure to be initialised.
 */
#define ACCEPT_WINDOW(b) if (!accept_window(&(S->b),w)) return (0); BGIWindow b=S->b;
#endif

/*! Similar to assert() but cannot be screened */
/* #define ASSERT(p) { if(0==(p)) EXPECTED_ERROR("Assertion failed:\n %s",#p); } - defined in error.h */

/*====================================================================*/
/* Other exported objects					      */
/*====================================================================*/
#include "extern.h"
			
/*! Maximum depth of nested includes in script files.*/
#define MAXDEPTH 10

/*! Current depth of include. */
EXTERN int depth;		

/*! Names of includes. */
EXTERN char inname[MAXDEPTH][MAXPATH];

/*! Coordinates in file of includes. */
EXTERN long inpline[MAXDEPTH], inppos[MAXDEPTH]; 

/*! Integer null value (-MAXINT+1). */
EXTERN int  INONE;

/*! Long null value (-MAXLONG+1). */
EXTERN long LNONE;

/*! Real (double) null value (-MAXREAL). */
EXTERN real RNONE;

/*! Real (double) "infinite" value (+MAXREAL). */
EXTERN REAL real_inf;

/*! Machine epsilon */
EXTERN real macheps;

#undef EXTERN

/*! "null" filename */
#define null "null"

/*! Command terminators */
#define TERMINATOR ";$"

/*! Escape character */
#define ESCAPE "#"

/*! Token separators.
 * NB: decimal point ".", and colon ":" are used in numbers and file ids, "," in functions 
 */
#define SEPARATORS " ;\t\r\n!$"

/*! Marks the beginning of a reference to a file. */
#define INCLUDEBEGIN '<'
/*! Marks the end of a reference to a file. */
#define INCLUDEEND '>'

/*! Marks the beginning of a reference to a string macro. */
#define PASTEBEGIN '['
/*! Marks the end of a reference to a string macro. */
#define PASTEEND ']'

/* reference to a system command output */
/*! marks the beginning of a reference to a system command output. */
#define CATCHBEGIN '`'
/*! marks the end of a reference to a system command output. */
#define CATCHEND '`'

/* Block of parameters */
/*! marks the beginning of a block of parameters. */
#define BLOCKBEGIN '{'
/*! marks the end of a block of parameters. */
#define BLOCKEND '}'

/* String with blanks */
/*! marks the end of a string with blanks. */
#define STRBEGIN '\"'
/*! marks the end of a string with blanks. */
#define STREND '\"'

/*! Linking parameter to a level of dynamic variables */
#define AT '@'

/* Links device parameter to a k-expression */
#define LINK '~'

/*! How default values are shown in res file */
#define DFLT "!"

