#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TENSION0        /* for a cardinal spline tension of 0 */
/*#define TENSION1 */        /* for a cardinal spline tension of 1 */


void cardinal(float *cpx, float *cpy, int ncp, float *x, float *y, int npoints)
    /* for each of the npoints values in the array x, returns in y the value of the
     * cardinal spline at that point.  */
{
    int     i,j;
    float    x1, x2, x3;


    for (i = 0; i < npoints ; i++) {
        /* find the segment of the spline containing the current point */
        j = 1;
        while (((cpx[j] > x[i]) || (cpx[j+1] < x[i])) && (j <= ncp - 2)){
            j++;
        }
        if (j == ncp - 2) {
            fprintf(stderr, "Error in cardinal: point outside spline bounds\n");
            exit(1);
        }
        /* compute the value of the spline at this point */
        x1 = (float) (x[i] - cpx[j]) / (float) (cpx[j+1] - cpx[j]);
        x2 = x1 * x1;
        x3 = x2 * x1;
#ifdef TENSION0
        y[i] = (cpy[j-1]*(-.5*x3 + x2 -.5*x1) + cpy[j]*(1.5*x3 - 2.5*x2 + 1) + 
                   cpy[j+1]*(-1.5*x3 + 2.0*x2 + .5*x1) + cpy[j+2]*(.5*x3 - .5*x2));  
#endif
#ifdef TENSION1
        y[i] = (cpy[j]*(2*x3 - 3*x2 + 1) + cpy[j+1]*(-2*x3 + 3.0*x2));
#endif
        if (y[i] < 0) y[i] = 0;  
    }
}
