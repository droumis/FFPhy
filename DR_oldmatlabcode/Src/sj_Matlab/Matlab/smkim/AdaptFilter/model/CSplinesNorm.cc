/* $Id: CSplinesNorm.cc,v 1.6 2006/07/04 15:30:19 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include <stdlib.h>
#include <string.h>

#include "../aux/hippoIO.h" 
#include "../aux/mexAux.h" 
#include "CSplinesNorm.h"
//#include "CardSplines.h"

/* ************************************************************
                          class CSplinesNorm
   ************************************************************ */

using namespace AFilter;

CSplinesNorm::CSplinesNorm()
{
    zero_out();
}

CSplinesNorm::CSplinesNorm(CSplinesNorm &cs)
    : CSplines(cs)
{
    zero_out();
}

CSplinesNorm::~CSplinesNorm()
{
    freeDNdTheta();
}

void CSplinesNorm::zero_out() 
{
    dNdTheta= 0;
    norm_last.factor= 0;
    norm_last.t= -1;
    norm_last.traj= -1;
}

void CSplinesNorm::allocDNdTheta() 
{
    hippo_Assert(ncpx, "control point are not initialized");
    freeDNdTheta();
    dNdTheta= new double*[nFct];
    for(int t= 0; t < nFct; t++) dNdTheta[t]= new double[ncpx[t]];
}

void CSplinesNorm::freeDNdTheta() 
{
    if(dNdTheta) {
        for(int t= 0; t < nFct; t++) delete[] dNdTheta[t];
        delete[] dNdTheta;
    }
    dNdTheta= 0;
}

void CSplinesNorm::calc_norm(int t, traj_type traj) const
{
    if(norm_last.t!=t || norm_last.traj!=traj)  {
        double n= fctCurr[traj]->integral();
        hippo_Assert(n > 0, "Integral became negative");
        norm_last.factor= 1.0/n;
        norm_last.t= t;
        norm_last.traj= traj;
//        hippo_Print(n);
    }
}

void CSplinesNorm::initObj()
{
    CSplines::initObj();
    nUpdate= 0;             // will be set to largest number of control points
    for(int traj= 0; traj < nFct; traj++)
        if(ncpx[traj]>nUpdate) nUpdate= ncpx[traj];

    // initialize the gradient of the norm w.r.t. the parameters,
    // the gradient can be pre-computed since it is independent of the 
    // parameters
    allocDNdTheta();
    const double s= 0.5/12.0;  // tension parameter t=0.5
    for(int traj= 0; traj < nFct; traj++) {
        int n= ncpx[traj]-3;
        double *x= cpx[traj];
        dNdTheta[traj][0]= -s*(x[2]-x[1]);	
        dNdTheta[traj][1]= (s+0.5)*(x[2]-x[1]) -s*(x[3]-x[2]);
        dNdTheta[traj][2]= (s+0.5)*(x[3]-x[1]) -s*(x[4]-x[3]);
        for(int j= 3; j < n; j++) 
            dNdTheta[traj][j]= 0.5*(x[j+1]-x[j-1]) 
                +s*(-x[j+2] +2*x[j+1] -2*x[j-1] +x[j-2]);
        dNdTheta[traj][n]= (s+0.5)*(x[n+1]-x[n-1]) -s*(x[n-1]-x[n-2]);
        dNdTheta[traj][n+1]= (s+0.5)*(x[n+1]-x[n]) -s*(x[n]-x[n-1]);
        dNdTheta[traj][n+2]= -s*(x[n+1]-x[n]);
    }
}

bool  CSplinesNorm::getParam(int t, double z[], int a[]) const
{
    hippo_Empty;
    return 0;
}

void CSplinesNorm::getAllParam(int t, double z[]) const
{
    hippo_Empty;
//    CSFunction::getAllParam(t, z);
//    calc_norm(t, traj);
//    for(int i= 0; i < ncpx[traj]; i++) out[i]*= norm_last.factor;
}

double CSplinesNorm::eval(double x, traj_type traj) const
{
//    hippo_Mark;
    double result= CSplines::eval(x,traj);
    if(traj < 0) return 1;
    calc_norm(tCurr, traj);
    if(result <0) result= 0;
    return norm_last.factor*result;
}

double CSplinesNorm::evalGrad(double x, traj_type traj, 
        double diff[]) const
{
    double result;
    result= this->CSplines::evalGrad(x, traj, diff);
    if(traj < 0)  return 1;

    calc_norm(tCurr, traj);
    double tmp[4];
    for(int i=0; i<4; i++) tmp[i]= diff[i];
    for(int i=0; i<nUpdate; i++) {
//        cerr << i << ": \t" << diff+i << "\t" << dNdTheta[traj]+i << "\n";
        diff[i]= -norm_last.factor*norm_last.factor*result*dNdTheta[traj][i];
    }

    traj_type tsegx;
    tsegx= csegx[tCurr];
    for(int i=0; i<4; i++)
        diff[tsegx-1+i]+= norm_last.factor*tmp[i];

    if(result <0) result= 0;
    return  result* norm_last.factor;
}

double CSplinesNorm::getIntegral(int timestep, traj_type traj) const
{
    hippo_Mark;
    moveTo(timestep);
    calc_norm(timestep, traj);
    return norm_last.factor*fctCurr[traj]->integral();
}

