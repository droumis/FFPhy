/* $Id: Moments2d.cc,v 1.6 2008/08/24 19:47:06 chengs Exp $
   
   implementation of class AFilter::Moments2d
*/

#include <math.h>

#include "Moments2d.h"
#include "../aux/mexAux.h"

/* ************************************************************
                          class Moments2d
   ************************************************************ */

using namespace AFilter;

Moments2d::Moments2d(const TData *d, const AdaptModel *m, const StatLimits *l)
    : VStat2d(d,m,l)
{
    zeroOut();
}

Moments2d::~Moments2d()
{
    dealloc();
}

void Moments2d::zeroOut()
{
    nbinx= nbiny= 0;
    p= 0;
    px= py= 0;
    x= y= 0;
}

void Moments2d::allocate()
{
    p= new double*[nbinx];
    for (int i= 0; i < nbinx; i++) p[i]= new double[nbiny];
    px= new double[nbinx];
    py= new double[nbiny];
    x= new double[nbinx];
    y= new double[nbiny];
}

void Moments2d::dealloc()
{
    for (int i= 0; i < nbinx; i++) delete[] p[i];
    delete[] p;
    delete[] px;
    delete[] py;
    delete[] x;
    delete[] y;
    zeroOut();
}

void Moments2d::run()
{
// calc mutual info.
//    hippo_Mark;
    const double tiny= 1e-15;

    // init analyses 
    nana= lim->getNAnalyses();
    init(nana, 7, lim->getNSteps());

    double time;
    traj_type traj;
    double norm;
    double meanx, meany, varx, vary, medianx, mediany, I;

    for (int ana=0; ana<nana; ana++) { // do all analyses
        dealloc();
        double minx= lim->getVarMin("linpos", ana);
        double maxx= lim->getVarMax("linpos", ana);
        nbinx= (int) floor((maxx-minx)/sx);
        double miny= lim->getVarMin("phase", ana);
        double maxy= lim->getVarMax("phase", ana);
        nbiny= (int) floor((maxy-miny)/sy);
        allocate();
        traj= (traj_type) round(lim->getVarMax("traj", ana));
        for (int is= 0; is < nsteps[ana]; is++) { // loop over timesteps
//    hippo_Mark;

            time= lim->getXMax(ana, is);
//                hippo_Print(time);
            
            if(time < dyn2->getMinTime()) time= dyn2->getMinTime();
            if(time > dyn2->getMaxTime()) time= dyn2->getMaxTime();

            // initialize for this calculation
            for (int i= 0; i < nbinx; i++) { px[i]= 0; x[i]= minx+(i+0.5)*sx; }
            for (int j= 0; j < nbiny; j++) { py[j]= 0; y[j]= miny+(j+0.5)*sy; }
            norm= 0;
	    
            // get firing rate distribution
            for (int i= 0; i < nbinx; i++) {
                for (int j= 0; j < nbiny; j++) {
                    p[i][j]= dyn2->evalAtTime(data->timeToIndex(time), x[i], y[j], traj) * sx * sy;
//                    hippo_Assert((p[i][j] >= 0), "firing rates cannot be negative");
                    if (p[i][j] < 0) p[i][j]= 0;
                    norm+= p[i][j];
//            hippo_Print(p[i][j]);
                }
            }

//            hippo_Print(norm);
            // normalize and get marginal distributions
            for (int i= 0; i < nbinx; i++) {
                for (int j= 0; j < nbiny; j++) {
                    if(norm > tiny) p[i][j]/= norm;
                    px[i]+= p[i][j];
                    py[j]+= p[i][j];
                }
            }

            // calculate mutual information
            meanx= meany= varx= vary= medianx= mediany= I= 0;
            for (int i= 0; i < nbinx; i++) meanx+= x[i]*px[i];
            for (int i= 0; i < nbinx; i++) varx+= sq(x[i]-meanx)*px[i];

            if(!circy) {
                for (int j= 0; j < nbiny; j++) meany+= y[j]*py[j];
                for (int j= 0; j < nbiny; j++) vary+= sq(y[j]-meany)*py[j];
            } else {
//                hippo_Mark;
                double tmpx=0, tmpy=0;
                for (int j= 0; j < nbiny; j++) {
//                    hippo_Print(py[j]);
//                    hippo_Print(y[j]);
                    tmpx+= cos(y[j])*py[j];
                    tmpy+= sin(y[j])*py[j];
                }
                meany= atan2(tmpy, tmpx);
                if(meany<0) meany+= 2*M_PI;
                vary= 1-sqrt(tmpx*tmpx+tmpy*tmpy);
            }

//                hippo_Mark;
            for (int i= 0; i < nbinx; i++) {
                for (int j= 0; j < nbiny; j++) {
                    if(p[i][j] > tiny) {
                       I+= p[i][j]*log2(p[i][j]/(px[i]*py[j]));
                    }

                }
            }

            // save result
            set(ana, 0, is, meanx);
            set(ana, 1, is, meany);
            set(ana, 2, is, varx);
            set(ana, 3, is, vary);
            set(ana, 4, is, medianx);
            set(ana, 5, is, mediany);
            set(ana, 6, is, I);
        } // for is
    }  // for ana
}

void Moments2d::setMexInput(const mxArray *in)
{
    VStat::setMexInput(in);
    dealloc();
    MX_FieldScalar(in, "sizex", sx, double);
    MX_FieldScalar(in, "sizey", sy, double);
    int tmp;
    MX_FieldScalarDefault(in, "circy", tmp, 0, int);
    circy= tmp ? true: false;
    allocate();
}
