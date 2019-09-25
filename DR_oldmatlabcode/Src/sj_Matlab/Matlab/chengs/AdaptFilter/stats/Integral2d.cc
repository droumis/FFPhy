/* $Id: Integral2d.cc,v 1.4 2008/08/24 19:47:06 chengs Exp $
   
   implementation of class AFilter::Integral2d
*/

#include <math.h>

#include "Integral2d.h"
#include "../aux/mexAux.h"

/* ************************************************************
                          class Integral2d
   ************************************************************ */

using namespace AFilter;

Integral2d::Integral2d(const TData *d, const AdaptModel *m, const StatLimits *l)
    : VStat2d(d,m,l)
{
    zeroOut();
}

Integral2d::~Integral2d()
{
}

void Integral2d::zeroOut()
{
    dx= dy= 0;
}

void Integral2d::run()
{
// calc integral

    // init analyses 
    nana= lim->getNAnalyses();
    init(nana, 1, lim->getNSteps());

    double time;
    traj_type traj;

    for (int ana=0; ana<nana; ana++) { // do all analyses
        double minx= lim->getVarMin("linpos", ana);
        double maxx= lim->getVarMax("linpos", ana);
        double miny= lim->getVarMin("phase", ana);
        double maxy= lim->getVarMax("phase", ana);
        traj= (traj_type) round(lim->getVarMax("traj", ana));
        for (int is= 0; is < nsteps[ana]; is++) { // loop over timesteps

            time= lim->getXMax(ana, is);
            
            if(time < dyn2->getMinTime()) time= dyn2->getMinTime();
            if(time > dyn2->getMaxTime()) time= dyn2->getMaxTime();
//            if(time < dyn2->getMinTime() || time > dyn2->getMaxTime()) {
//                set(ana, 0, is, mxGetNaN());
//                continue;
//            }
	    
            double I= 0;
            // get firing rate distribution @@ use better algorithm
            for (double x= minx; x < maxx; x+= dx)
                for (double y= miny; y < maxy; y+= dy)
                    I+= dyn2->evalAtTime(data->timeToIndex(time), x, y, traj);

            I*= (dx*dy);
            set(ana, 0, is, I);
        }
    }
}

void Integral2d::setMexInput(const mxArray *in)
{
    VStat::setMexInput(in);
    MX_FieldScalar(in, "dx", dx, double);
    MX_FieldScalar(in, "dy", dy, double);
//    hippo_Print(dx);
//    hippo_Print(dy);
}
