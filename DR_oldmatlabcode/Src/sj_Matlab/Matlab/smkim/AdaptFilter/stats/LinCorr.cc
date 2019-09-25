/* $Id: LinCorr.cc,v 1.6 2008/08/24 19:47:06 chengs Exp $
   
   implementation of class AFilter::LinCorr
*/

#include <math.h>

#include "LinCorr.h"
#include "../aux/mexAux.h"

/* ************************************************************
                          class LinCorr
   ************************************************************ */

using namespace AFilter;

LinCorr::LinCorr(const TData *d, const AdaptModel *m, const StatLimits *l)
    : VStat2d(d,m,l)
{
//    zeroOut();
    x=y=0;
    nbinx= nbiny= 0;
    p= 0;
    px= py= 0;
    sx= sy= 0;
    minx= maxx= 0;
    miny= maxy= 0;
}

LinCorr::~LinCorr()
{
//    dealloc();
    if(p) {
        for (int i= 0; i < nbinx; i++) delete[] p[i];
        delete[] p;
    }
    if(px) delete[] px;
    if(py) delete[] py;
    if(x) delete[] x;
    if(y) delete[] y;
}

void LinCorr::zeroOut()
{
    x=y=0;
    nbinx= nbiny= 0;
    p= 0;
    px= py= 0;
    sx= sy= 0;
    minx= maxx= 0;
    miny= maxy= 0;
}

void LinCorr::allocate()
{
    p= new double*[nbinx];
    for (int i= 0; i < nbinx; i++) p[i]= new double[nbiny];
    px= new double[nbinx];
    py= new double[nbiny];
    x= new double[nbinx];
    y= new double[nbiny];
}

void LinCorr::dealloc()
{
    if(p) {
        for (int i= 0; i < nbinx; i++) delete[] p[i];
        delete[] p;
    }
    if(px) delete[] px;
    if(py) delete[] py;
    if(x) delete[] x;
    if(y) delete[] y;
//    zeroOut();
}

void LinCorr::setXY()
{
    for(int i= 0; i < nbinx; i++) x[i]= minx+(i+0.5)*sx;
    for(int j= 0; j < nbiny; j++) y[j]= miny+(j+0.5)*sy;
}

void LinCorr::normNmarg()
{
    // initialize for this calculation
    double norm= 0;
    for (int i= 0; i < nbinx; i++) px[i]= 0;
    for (int j= 0; j < nbiny; j++) py[j]= 0;

    // calc norm
    for (int i= 0; i < nbinx; i++)
        for (int j= 0; j < nbiny; j++)
            norm+= p[i][j];
    // normalize and get marginal distributions
    for (int i= 0; i < nbinx; i++) {
        for (int j= 0; j < nbiny; j++) {
            p[i][j]/= norm;
            px[i]+= p[i][j];
            py[j]+= p[i][j];
        }
    }
}

double LinCorr::calcCorrCoef()
{
    // calc means 
    double mx=0, my=0;
    for (int i= 0; i < nbinx; i++) mx+= x[i]*px[i]; 
    for (int j= 0; j < nbiny; j++) my+= y[j]*py[j]; 

    // calc variances
    double Vx=0, Vy=0;
    for (int i= 0; i < nbinx; i++) Vx+= sq(x[i]-mx) *px[i]; 
    for (int j= 0; j < nbiny; j++) Vy+= sq(y[j]-my) *py[j]; 

    // calculate linear correlation coeff
    double C=0;
    for (int i= 0; i < nbinx; i++) {
        for (int j= 0; j < nbiny; j++) {
            C+= (x[i]-mx)*(y[j]-my)*p[i][j];
        }
    }
    C= C/sqrt(Vx*Vy);
    return C;
}


void LinCorr::run()
{
// calc mutual info.

    // init analyses 
    nana= lim->getNAnalyses();
    init(nana, 1, lim->getNSteps());

    double time;
    traj_type traj;

    for (int ana=0; ana<nana; ana++) { // do all analyses
        dealloc();
        minx= lim->getVarMin("linpos", ana);
        maxx= lim->getVarMax("linpos", ana);
        nbinx= (int) floor((maxx-minx)/sx);
        miny= lim->getVarMin("phase", ana);
        maxy= lim->getVarMax("phase", ana);
        nbiny= (int) floor((maxy-miny)/sy);
        allocate();

        setXY();

        traj= (traj_type) round(lim->getVarMax("traj", ana));
        for (int is= 0; is < nsteps[ana]; is++) { // loop over timesteps

            time= lim->getXMax(ana, is);
            
            if(time < dyn2->getMinTime()) time= dyn2->getMinTime();
            if(time > dyn2->getMaxTime()) time= dyn2->getMaxTime();
//            if(time < dyn2->getMinTime() || time > dyn2->getMaxTime()) {
//                set(ana, 0, is, mxGetNaN());
//                continue;
//            }
	    
            // get firing rate distribution
            for (int i= 0; i < nbinx; i++) {
                for (int j= 0; j < nbiny; j++) {
                    p[i][j]= dyn2->evalAtTime(data->timeToIndex(time), x[i], y[j], traj);
                }
            }
            normNmarg();

            // calc and save corr. coef.
            set(ana, 0, is, calcCorrCoef());
        }
    }
}

void LinCorr::setMexInput(const mxArray *in)
{
    VStat::setMexInput(in);
    MX_FieldScalar(in, "sizex", sx, double);
    MX_FieldScalar(in, "sizey", sy, double);
}
