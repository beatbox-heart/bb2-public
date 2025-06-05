#include <stdio.h>
#include <stdlib.h>
/* Create the stepwise colormap for fig 3.           */
/* 1 millisecond corresponds to input pixel value 2. */
/* The colours are chosen arbitrarily thus           */
/* do not coincide with the paper.                   */

int main(void) {
  /*           0   1   2   3   4   5   6   7   8   9  10, 11, 12 */
  int r[13]={  0,  0,  0,  0,  0,  0, 10, 50, 90,130,170,210,250};
  int g[13]={ 10, 50, 90,130,170,210,250,210,170,130, 90, 50, 10};
  int b[13]={250,210,170,130, 90, 50, 10,  0,  0,  0,  0,  0,  0};
  int i, j;
  printf("P3\n");
  printf("255 1\n");
  printf("255\n");
  printf("  0  0  0\n"); /* activation time 0 == never activated so black */
  for (i=1;i<256;i++) {
    if (i%20 >=18) 
      printf("  0  0  0\n"); /* black contour */
    else {
      j=i/20;
      printf("%d %d %d\n",r[j],g[j],b[j]); /* colour band */
    }
  }
  /* the remaining 2 ms would have been cropped to within 255 anyway */
  return 0;
}
