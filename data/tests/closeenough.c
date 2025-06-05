#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void error(char *fmt, ...) {
  va_list argptr;
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);
  exit(1);
}

/* #define error(fmt,...) {fprintf(stderr,fmt,__VA_ARGS__);exit(1);} */


static double verysmall=2.e-4;
#define BUFLEN 4096
#define SEP " \t\r\b\n,:;"

int main(int argc, char **argv) {
  FILE *in1, *in2;
  char buf1[BUFLEN], buf2[BUFLEN];
  char *p1, *p2, *last1, *last2;
  int key;
  double val1, val2, diff;
  int line, field, fieldok;
  
  if (argc!=3) error("%s: need 2 arguments not %d\n",argv[0],argc-1);
  /* printf("comparing %s against %s\n",argv[1],argv[2]); */
  if (NULL==(in1=fopen(argv[1],"r"))) error("could not open %s for reading\n",argv[1]);
  if (NULL==(in2=fopen(argv[2],"r"))) error("could not open %s for reading\n",argv[2]);

  for (line=1; !feof(in1) && !feof(in2); line++) {
    /* fprintf(stderr,"line %d\n",line); */
    if (NULL==fgets(buf1,BUFLEN,in1) && ferror(in1)) error("error reading line %d from %s\n",line,argv[1]);
    if (NULL==fgets(buf2,BUFLEN,in2) && ferror(in2)) error("error reading line %d from %s\n",line,argv[2]);
    fieldok=1;
    for(field=1;fieldok;field++) {
      /* fprintf(stderr,"field %d\n",field); */
      p1=strtok_r(field==1?buf1:NULL,SEP,&last1);
      p2=strtok_r(field==1?buf2:NULL,SEP,&last2);
      key=(2*(p1!=NULL) + (p2!=NULL));
      /* fprintf(stderr,"p1=%p p2=%p key=%d line=%d field=%d\n",p1,p2,key,line,field); */
      switch (key) {
      case 0:
	fieldok=0;
	break;
      case 1:
	error("line %d in %s has field %d but %s does not have it\n",line,argv[2],field,argv[1]);
      case 2:
	error("line %d in %s has field %d but %s does not have it\n",line,argv[1],field,argv[2]);
      case 3:
	if (1!=sscanf(p1,"%lg",&val1)) 
	  error("cannot read field %d in line %d of %s as a real\n",field,line,argv[1]);
	if (1!=sscanf(p2,"%lg",&val2)) 
	  error("cannot read field %d in line %d of %s as a real\n",field,line,argv[2]);
	diff=val1-val2;
	diff=(diff>0)?diff:-diff;
	if (diff>verysmall)
	  error("values in line %d field %d differ by %lg > %lg\n",line,field,diff,verysmall);
	fieldok=1;
      }	
    }
  }
  fclose(in2);
  fclose(in1);
  return 0;
}
