#ifndef DEFS_H
#define DEFS_H

/* $Id: defs.h,v 1.2 2008/06/03 17:57:56 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

typedef short int traj_type;
typedef short int count_type;

const double tiny= 1e-10;

#define IS_ZERO(x) ((x) < tiny && (x) > -tiny)

// from /usr/local/matlab-6.5/extern/include/matrix.h:
struct mxArray_tag;
typedef struct mxArray_tag mxArray;


#endif   // DEFS_H
