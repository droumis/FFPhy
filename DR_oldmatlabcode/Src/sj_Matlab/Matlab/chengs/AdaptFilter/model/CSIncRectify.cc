/* $Id: CSIncRectify.cc,v 1.2 2008/06/27 00:21:28 chengs Exp $
   
   Sen Cheng, Wed Oct  6 08:39:54 PDT 2004
   implementation of class CSIncRectify
*/

#include "CSIncRectify.h"

/* ************************************************************
                          class CSIncRectify
   ************************************************************ */

using namespace AFilter;

CSIncRectify::CSIncRectify()
{
}

CSIncRectify::~CSIncRectify()
{
}

// move to time t
void CSIncRectify::moveTo(int t) const
{
    hippo_Assert(t >= startindex && t <= endindex, "t out of bounds");

    if(t== tCurr) return;

    double tmp;
    if(t > tCurr) { // move forward
        for( ; tCurr < t; tCurr+= 1) {
            super->getActiveIndices(tCurr, nt, trajs, inds);
            if(nt <= 0) continue;
            double *u= updates+nUpdate*tCurr;
            for(int it= 0; it < nt; it++) 
                for(int i= 0; i < nUpdate; i++)  {
                    tmp= param[trajs[it]][inds[i]]+= u[i];
                    if(tmp<0) {
                        u[i]-= tmp;
                        param[trajs[it]][inds[i]]= 0;
                    }
                }
        }
    } else {       // move backwards
        for(tCurr-= 1; tCurr >= t; tCurr-= 1) {
            super->getActiveIndices(tCurr, nt, trajs, inds);
            if(nt <= 0) continue;
            double *u= updates+nUpdate*tCurr;
            for(int it= 0; it < nt; it++) 
                for(int i= 0; i < nUpdate; i++) {
                    tmp= param[trajs[it]][inds[i]]-= u[i];
                    if(tmp<0) {
                        u[i]+= tmp;
                        param[trajs[it]][inds[i]]= 0;
                    }
            }
        }
        tCurr+= 1;
    }

    hippo_Assert(t==tCurr, "current time messed up");
}

