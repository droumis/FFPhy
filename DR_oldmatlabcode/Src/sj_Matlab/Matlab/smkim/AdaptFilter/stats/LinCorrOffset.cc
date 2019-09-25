/* $Id: LinCorrOffset.cc,v 1.1 2007/01/12 04:57:17 chengs Exp $
   
   implementation of class AFilter::LinCorrOffset
*/

#include <math.h>

#include "LinCorrOffset.h"
#include "../aux/mexAux.h"

/* ************************************************************
                          class LinCorrOffset
   ************************************************************ */

using namespace AFilter;

const int fac[]= {-1,1,-1,1};
void LinCorrOffset::setOffset()
{
    double C= 0, maxC= 0;
    int bestshift= 0;
    for(int s=0; s<nbiny; ++s) {
        for (int j= 0; j < nbiny; j++) 
            y[j]= miny+(mod(j+s,nbiny)+0.5)*sy;
        C= fabs(calcCorrCoef());
        if(C > maxC) {
            maxC= C;
            bestshift= s;
        }
    }
    offset= y[bestshift];
    for (int j= 0; j < nbiny; j++) 
        y[j]= miny+(mod(j+bestshift,nbiny)+0.5)*sy;
    printf("traj= %d, shift= %d, opt C= %.2f\n", traj, 
            bestshift, calcCorrCoef());
}


//void LinCorrOffset::setOffset()
//{
//    double offset= lim->getVarMax("offset", ana);
//    for (int j= 0; j < nbiny; j++) {
//        y[j]= fmod((yreal[j]-offset), 2*M_PI);
//        if(y[j]<0) y[j]+= 2*M_PI;
//        cerr << j << ":\t " << yreal[j] << "\t" << y[j] << endl;
//    }
//}

