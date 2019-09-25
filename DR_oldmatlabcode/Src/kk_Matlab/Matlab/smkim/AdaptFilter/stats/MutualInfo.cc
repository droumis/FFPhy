/* $Id: MutualInfo.cc,v 1.5 2008/08/24 19:47:06 chengs Exp $
   
   implementation of class AFilter::MutualInfo
*/

#include <math.h>

#include "MutualInfo.h"
#include "../aux/mexAux.h"

/* ************************************************************
                          class MutualInfo
   ************************************************************ */

using namespace AFilter;

MutualInfo::MutualInfo(const TData *d, const AdaptModel *m, const StatLimits *l)
    : VStat2d(d,m,l)
{
    zeroOut();
}

MutualInfo::~MutualInfo()
{
    dealloc();
}

void MutualInfo::zeroOut()
{
    nbinx= nbiny= 0;
    p= 0;
    px= py= 0;
}

void MutualInfo::allocate()
{
    p= new double*[nbinx];
    for (int i= 0; i < nbinx; i++) p[i]= new double[nbiny];
    px= new double[nbinx];
    py= new double[nbiny];
}

void MutualInfo::dealloc()
{
    for (int i= 0; i < nbinx; i++) delete[] p[i];
    delete[] p;
    delete[] px;
    delete[] py;
    zeroOut();
}

void MutualInfo::run()
{
// calc mutual info.
//    hippo_Mark;
    const double tiny= 1e-15;

    // init analyses 
    nana= lim->getNAnalyses();
    init(nana, 1, lim->getNSteps());

    double time;
    traj_type traj;
    double norm;
    double I;

    for (int ana=0; ana<nana; ana++) { // do all analyses
        double minx= lim->getVarMin("linpos", ana);
        double maxx= lim->getVarMax("linpos", ana);
        double dx= (maxx-minx)/nbinx;
        double miny= lim->getVarMin("phase", ana);
        double maxy= lim->getVarMax("phase", ana);
        double dy= (maxy-miny)/nbiny;
        traj= (traj_type) round(lim->getVarMax("traj", ana));
        for (int is= 0; is < nsteps[ana]; is++) { // loop over timesteps

            time= lim->getXMax(ana, is);
            
            if(time < dyn2->getMinTime()) time= dyn2->getMinTime();
            if(time > dyn2->getMaxTime()) time= dyn2->getMaxTime();
//            if(time < dyn2->getMinTime() || time > dyn2->getMaxTime()) {
//                set(ana, 0, is, mxGetNaN());
//                continue;
//            }
            // initialize for this calculation
            for (int i= 0; i < nbinx; i++) px[i]= 0;
            for (int j= 0; j < nbiny; j++) py[j]= 0;
            norm= 0;
	    
            // get firing rate distribution
            for (int i= 0; i < nbinx; i++) {
                for (int j= 0; j < nbiny; j++) {
                    p[i][j]= dyn2->evalAtTime(data->timeToIndex(time), minx+(i+0.5)*dx, miny+(j+0.5)*dy, traj) * dx * dy;
                    norm+= p[i][j];
                    hippo_Assert((p[i][j] >= 0), "firing rates cannot be negative");
                }
            }

            // normalize and get marginal distributions
            for (int i= 0; i < nbinx; i++) {
                for (int j= 0; j < nbiny; j++) {
                    if(norm > tiny) p[i][j]/= norm;
                    px[i]+= p[i][j];
                    py[j]+= p[i][j];
                }
            }

            // calculate mutual information
            I=0;
            for (int i= 0; i < nbinx; i++) {
                for (int j= 0; j < nbiny; j++) {
                    if(p[i][j] > tiny) {
                       I+= p[i][j]*log2(p[i][j]/(px[i]*py[j]));
//                       if(px[i] < tiny || py[j] < tiny)
//                cerr << px[i] << ", " << py[j] << " | \t"; 
                    }
                }
            }
//            cerr << I << "\n";
            // save result
            set(ana, 0, is, I);
//            hippo_Print(I);
        }
    }
}

void MutualInfo::setMexInput(const mxArray *in)
{
    VStat::setMexInput(in);
    dealloc();
    MX_FieldScalar(in, "nbinx", nbinx, int);
    MX_FieldScalar(in, "nbiny", nbiny, int);
    allocate();
}
