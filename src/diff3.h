switch(dim) {
 case 1: {
#define _(x) (x),(0),(0)
#define STC STC1D(ANISO,MANY,VARIA)
   DEBUG("1D version with GEOM=%d ANISO=%d MANY=%d VARIA=%d\n",GEOM,ANISO,MANY,VARIA);   
#include "diff1d.h"
#undef STC
#undef _
 } break; /* case 1 */
 case 2: {
#define _(x,y) (x),(y),(0)
#define STC STC2D(ANISO,MANY,VARIA)
   DEBUG("2D version with GEOM=%d ANISO=%d MANY=%d VARIA=%d\n",GEOM,ANISO,MANY,VARIA);   
#include "diff2d.h"
#undef STC
#undef _
 } break; /* case 2 */
 case 3: {
#define _(x,y,z) (x),(y),(z)
#define STC STC3D(ANISO,MANY,VARIA)
   DEBUG("3D version with GEOM=%d ANISO=%d MANY=%d VARIA=%d\n",GEOM,ANISO,MANY,VARIA);   
#include "diff3d.h"
#undef STC
#undef _
 } break; /* case 3 */
 default:
   EXPECTED_ERROR("unexpected dimension %d",dim);
   break;
 } /* switch dim */
