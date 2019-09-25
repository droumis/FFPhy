/* $Id: LinCorrShift.cc,v 1.3 2008/08/24 19:47:06 chengs Exp $
   
   implementation of class AFilter::LinCorrShift
*/

#include <math.h>

#include "LinCorrShift.h"
#include "../aux/mexAux.h"

/* ************************************************************
                          class LinCorrShift
   ************************************************************ */

using namespace AFilter;

LinCorrShift::LinCorrShift(const TData *d, const AdaptModel *m, const StatLimits *l)
        : LinCorr(d,m,l) 
{
    yreal=0;
};

LinCorrShift::~LinCorrShift() 
{
    if(yreal) delete[] yreal;

};

void LinCorrShift::zeroOut()
{
    LinCorr::zeroOut();
    yreal=0;
}

void LinCorrShift::allocate()
{
    LinCorr::allocate();
    yreal= new double[nbiny];
}

void LinCorrShift::dealloc()
{
    if(yreal) delete[] yreal;
    LinCorr::dealloc();
}


const int fac[]= {-1,1,-1,1};
void LinCorrShift::setOffset()
{
    double C= 0, maxC= 0;
    int bestshift= 0;
    for(int s=0; s<nbiny; ++s) {
        for (int j= 0; j < nbiny; j++) 
            y[j]= miny+(mod(j+s,nbiny)+0.5)*sy;
//        C= calcCorrCoef();
//        cerr << "s= " << s << "\t C= " << C << endl;
        C= fac[traj]*(calcCorrCoef());
//        C= fabs(calcCorrCoef());
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


/*
// set off set to min activity
void LinCorrShift::setOffset()
{
    int bestshift= 0;
    double min= py[0];
    for (int j= 1; j < nbiny; j++) 
        if(py[j] < min) {
            min= py[j];
            bestshift= j;
        }
    for (int j= 0; j < nbiny; j++) 
        y[j]= miny+(mod(j+bestshift,nbiny)+0.5)*sy;
        
    printf("traj= %d, shift= %d, py= %.2f, opt C= %.2f\n", traj, 
            bestshift, min, calcCorrCoef());
}
*/

void LinCorrShift::run()
{
    // calc mutual info.

    // init analyses 
    nana= lim->getNAnalyses();
    init(nana, 2, lim->getNSteps());

    double time;

    for (ana=0; ana<nana; ana++) { // do all analyses
        traj= (traj_type) round(lim->getVarMax("traj", ana));
        dealloc();
        minx= lim->getVarMin("linpos", ana);
        maxx= lim->getVarMax("linpos", ana);
        nbinx= (int) floor((maxx-minx)/sx);
        miny= lim->getVarMin("phase", ana);
        maxy= lim->getVarMax("phase", ana);
        nbiny= (int) floor((maxy-miny)/sy);
        allocate();


        setXY();
        for (int j= 0; j < nbiny; j++) yreal[j]= y[j];

        // get average
        for (int i= 0; i < nbinx; i++)
            for (int j= 0; j < nbiny; j++)
                p[i][j]= 0;
        for (int is= 0; is < nsteps[ana]; is++) { // loop over timesteps
            time= lim->getXMax(ana, is);
            if(time < dyn2->getMinTime()) time= dyn2->getMinTime();
            if(time > dyn2->getMaxTime()) time= dyn2->getMaxTime();

            for (int i= 0; i < nbinx; i++) {
                for (int j= 0; j < nbiny; j++) {
                    p[i][j]+= dyn2->evalAtTime(data->timeToIndex(time), x[i], y[j], traj);
                }
            }
        }  // get average
        normNmarg();
        setOffset();

        for (int is= 0; is < nsteps[ana]; is++) { // loop over timesteps

            time= lim->getXMax(ana, is);

            if(time < dyn2->getMinTime()) time= dyn2->getMinTime();
            if(time > dyn2->getMaxTime()) time= dyn2->getMaxTime();

            // get firing rate distribution
            for (int i= 0; i < nbinx; i++) {
                for (int j= 0; j < nbiny; j++) {
                    p[i][j]= dyn2->evalAtTime(data->timeToIndex(time), x[i], yreal[j], traj);
                }
            }
            normNmarg();
            // calc and save corr. coef.
            set(ana, 0, is, calcCorrCoef());
            set(ana, 1, is, offset);
        }
    }
}

