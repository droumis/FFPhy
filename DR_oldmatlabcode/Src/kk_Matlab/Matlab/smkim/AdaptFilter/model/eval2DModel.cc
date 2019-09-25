/* $Id: eval2DModel.cc,v 1.3 2008/08/24 19:47:05 chengs Exp $
   eval2DModel.c 
  
*/
#include <string.h>
#include <math.h>
#include <iostream>

#include "../aux/defs.h"
#include "../aux/mexAux.h"
#include "../aux/hippoIO.h"

#include "../model/modelFactory.h"
#include "../model/PosPhase_Isi.h"


using namespace std;
using namespace AFilter;

void Usage(void) ;

/******************************************************************************
  INTERFACE FUNCTION to be called from matlab
*/

/* 
   [output]= eval2DModel(data, model, opt)
   output: x, y, z, t
   opt:    traj, nt, nx, ny
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const double tiny= 1e-10;
    /* Check numbers of arguments */
    if ( (nlhs != 1) | ((nrhs != 3))) {
        mexErrMsgTxt("incorrect number of parameters, type \"help eval2DModel\" for usage\n");
    }

    // read inputs and initialize objects

    TData *data= new TData(prhs[0]);
    AdaptModel *am= allocModel(prhs[1], data);
    VDynamics *dyn= dynamic_cast<IsiModel*>(am)->getSpatialModel();
    hippo_Assert(dyn->getDim()==2 | dyn->getDim()==0, "This function requires 2-dim spatial model.");
    VDynamics2d * m= dynamic_cast<VDynamics2d*>(dyn);

    int *traj, ntraj, nt, nx, ny;
    double *t;
    const mxArray *opts= prhs[2];
    mxArray *mxtmp;
    double *dptr;

    // read in options
    mxtmp= mxGetField(opts, 0, "traj");
    hippo_Assert(mxtmp, "couldn't read field 'traj'");
    MX_Assign(mxtmp, dptr);
    ntraj= mxGetM(mxtmp)*mxGetN(mxtmp);
    traj= new int[ntraj];
    for(int i=0; i<ntraj; i++) {
        traj[i]= (int) round(dptr[i]); 
//        cerr << traj[i] << endl;
    }
//    MX_FieldScalar(opts, "nt", nt , int);
    MX_FieldAssign(opts, "t", t, nt);
    MX_FieldScalar(opts, "nx", nx , int);
    MX_FieldScalar(opts, "ny", ny , int);

//    hippo_Print(ntraj)
//    hippo_Print(nt)
//    for(int it=0; it<nt; it++) {
//        cerr << it << "\t" << t[it] << endl;
//    }

    // evaluate model
//    double *t= new double[nt], *x= new double[nx], *y= new double[ny], 
    double *x= new double[nx], *y= new double[ny], 
        *z= new double[nx*ny];
//    double mint= m->getMinTime()+tiny, maxt= m->getMaxTime()-tiny, 
    double minx= m->getMinX()+tiny, maxx= m->getMaxX()-tiny, 
           miny= m->getMinY()+tiny, maxy= m->getMaxY()-tiny;
//    hippo_Print(minx);
//    hippo_Print(maxx);

//    if(nt==1) t[0]= mint;
//    else for(int it=0; it<nt; it++) t[it]= it*(maxt-mint)/(nt-1) + mint;
    if(nx==1) x[0]= minx;
    else for(int ix=0; ix<nx; ix++) x[ix]= ix*(maxx-minx)/(nx-1) + minx;
    if(ny==1) y[0]= miny;
    else for(int iy=0; iy<ny; iy++) y[iy]= iy*(maxy-miny)/(ny-1) + miny;

    // temporary storage for z
    mxArray *zout= mxCreateCellMatrix(ntraj, nt);

//    double dx, dy, dt;
    for(int itraj=0; itraj<ntraj; itraj++) {
        for(int it=0; it<nt; it++) {
            for(int ix=0; ix<nx; ix++) {
                for(int iy=0; iy<ny; iy++) {
                      z[iy*nx+ix]= m->evalAtTime(data->timeToIndex(t[it]), x[ix], y[iy], traj[itraj]); 
                }
            }
            mxArray *mxtmp;
            double *dblPtr;
            MX_CreateDoubleAssign(mxtmp, nx, ny, dblPtr);
            memcpy(dblPtr, z, nx*ny * sizeof(double));
            mxSetCell(zout, it*ntraj+itraj, mxtmp);
        }
    }

    // generate output structures
    mxArray *out= mxCreateStructMatrix(1,1, 0,  0);
    mxAddField(out, "t");
    MX_AssignDouble(out, "t", 1, nt, t);
    mxAddField(out, "x");
    MX_AssignDouble(out, "x", 1, nx, x);
    mxAddField(out, "y");
    MX_AssignDouble(out, "y", 1, ny, y);
    mxAddField(out, "z");
    mxSetField(out, 0, "z", zout);

    // set output
    plhs[0]= out;

    // free allocated memory
//    delete[] t; 
    delete[] x; delete[] y; delete[] z;
    delete m;
    delete data;
    return;
}
