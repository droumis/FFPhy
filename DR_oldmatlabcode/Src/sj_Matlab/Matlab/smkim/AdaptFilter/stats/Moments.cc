/* $Id: Moments.cc,v 1.3 2008/08/24 19:47:06 chengs Exp $
   
   implementation of class AFilter::Moments
*/

#include <math.h>

#include "Moments.h"
#include "../aux/mexAux.h"

/* ************************************************************
                          class Moments
   ************************************************************ */

using namespace AFilter;

Moments::Moments(const TData *d, const AdaptModel *m, const StatLimits *l)
    : VStat1d(d,m,l)
{
//    hippo_Mark;
    zeroOut();
}

Moments::~Moments()
{
    dealloc();
}

void Moments::zeroOut()
{
    nbinx= 0;
    p= 0;
    x= 0;
}

void Moments::allocate()
{
    p= new double[nbinx];
    x= new double[nbinx];
}

void Moments::dealloc()
{
    delete[] p;
    delete[] x;
    zeroOut();
}

void Moments::run()
{
// calc mutual info.
//    hippo_Mark;
    const double tiny= 1e-15;

    // init analyses 
    nana= lim->getNAnalyses();
    init(nana, 4, lim->getNSteps());

    double time;
    traj_type traj;
    double area, mean, std, skew;

    for (int ana=0; ana<nana; ana++) { // do all analyses
        dealloc();
        double minx= lim->getVarMin("linpos", ana);
        double maxx= lim->getVarMax("linpos", ana);
        nbinx= (int) floor((maxx-minx)/sx);
        allocate();
        traj= (traj_type) round(lim->getVarMax("traj", ana));
        for (int is= 0; is < nsteps[ana]; is++) { // loop over timesteps
//    hippo_Mark;

            time= lim->getXMax(ana, is);
//                hippo_Print(time);
            
            if(time < dyn1->getMinTime()) time= dyn1->getMinTime();
            if(time > dyn1->getMaxTime()) time= dyn1->getMaxTime();

            // initialize for this calculation
            for (int i= 0; i < nbinx; i++) { x[i]= minx+(i+0.5)*sx; }

            area= 0;
	    
            // get firing rate distribution
            for (int i= 0; i < nbinx; i++) {
                p[i]= dyn1->evalAtTime(data->timeToIndex(time), x[i], traj) * sx;
//                hippo_Assert((p[i]>= 0), "firing rates should not be negative");
                if (p[i] < 0) p[i]= 0;
                area+= p[i]; // calculate area
            }

//            hippo_Print(area);
            if (area<tiny) area= tiny;

            // normalize firing rate function to get PDF
            for (int i= 0; i < nbinx; i++) p[i]/= area;

            // calculate mean, std, skewness
            mean= std= skew= 0;
            for (int i= 0; i < nbinx; i++) mean+= x[i]*p[i];
//            hippo_Print(mean);
            for (int i= 0; i < nbinx; i++) {
                std+= sq(x[i]-mean)*p[i];
                skew+= pow((x[i]-mean),3)*p[i];
            }
            std= sqrt(std);
            skew= skew/pow(std,3);

//            hippo_Print(std);
//            hippo_Print(skew);

            // save result
            set(ana, 0, is, area);
            set(ana, 1, is, mean);
            set(ana, 2, is, std);
            set(ana, 3, is, skew);
        } // for is
    }  // for ana
}

void Moments::setMexInput(const mxArray *in)
{
    VStat1d::setMexInput(in);
    dealloc();
    MX_FieldScalar(in, "sizex", sx, double);
    allocate();
}
