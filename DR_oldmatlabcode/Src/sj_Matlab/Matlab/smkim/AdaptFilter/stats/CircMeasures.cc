/* $Id: CircMeasures.cc,v 1.3 2008/08/24 19:47:06 chengs Exp $
   
   implementation of class AFilter::CircMeasures
*/

#include <math.h>

#include "CircMeasures.h"
#include "../aux/mexAux.h"
#include "../aux/numerics.h"

/* ************************************************************
                          class CircMeasures
   ************************************************************ */

using namespace AFilter;

CircMeasures::CircMeasures(const TData *d, const AdaptModel *m, const StatLimits *l)
    : LinCorr(d,m,l)
{
}

CircMeasures::~CircMeasures()
{
}


void CircMeasures::run()
{
    // init analyses 
    nana= lim->getNAnalyses();
    init(nana, 5, lim->getNSteps());

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
	    
            // get firing rate distribution
            for (int i= 0; i < nbinx; i++) {
                for (int j= 0; j < nbiny; j++) {
                    p[i][j]= dyn2->evalAtTime(data->timeToIndex(time), x[i], y[j], traj);
                }
            }
            normNmarg();

            // calc variables
            double change= 0, D= 0, mTheta= 0;
            double prevm=0, m, d, prevd= 0, tmp;
            double change2= 0;
            double change3= 0;
            for (int i= 0; i < nbinx; i++) {
                circstats(nbiny, y, p[i], m, d);
//                if(i) change+= angdist(m, prevm);
                if(i) {
                    tmp= angdist(m, prevm);
                    change+= tmp;
                    change2+= min((1-d),1-prevd)*tmp;
                    change3+= (2-d-prevd)/2*tmp;
                }
                prevm= m;
                prevd= d;
                D+= px[i]*d;
//                printf("i= %2d, m= %.4f\n", i, m);
            }
//            change/= (nbinx-1);
            circstats(nbiny, y, py, mTheta, d);


            // set outputs
            set(ana, 0, is, change); // mean changes
            set(ana, 1, is, D); // mean dispersion
            set(ana, 2, is, mTheta); // overall mean phase
            set(ana, 3, is, change2); // mean changes
            set(ana, 4, is, change3); // mean changes
        }
    }
}

