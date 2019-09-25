/* $Id: CSplines2d.cc,v 1.10 2008/08/24 19:47:03 chengs Exp $
   File name   : model/CSplines2d.cc
   authors     : Sen Cheng
   created     : Fri 08 Oct 2004 04:16:33 PM PDT
  
 */

#include <stdlib.h>
#include <string.h>

#include "../aux/hippoIO.h" 
#include "../aux/mexAux.h" 
#include "CSplines2d.h"
//#include "CardSplines2d.h"

/* ************************************************************
                          class CSplines2d
   ************************************************************ */

using namespace AFilter;

CSplines2d::CSplines2d()
{
    // init all data members
    zero_out();
}

CSplines2d::CSplines2d(CSplines2d &)
{
    // stil needs to be implemented
    zero_out();
    hippo_Empty;
}

CSplines2d::~CSplines2d()
{
    dealloc();
}

void CSplines2d::zero_out()
{
    nUpdate= 16; T= 0; nFct=0; id=0;

    ncpx= ncpy= 0; cpx= cpy= 0; csegx= csegy= 0;

    allocated.controls= false;
    periodicY= false;
}

void CSplines2d::dealloc() 
{
    if(allocated.controls) {
        if(nFct>1) {
            if (ncpx) delete[] ncpx; 
            if (ncpy) delete[] ncpy; 
            if (cpx) delete[] cpx; 
            if (cpy) delete[] cpy; 
        } else if(nFct==1) {
            if (ncpx) delete ncpx; 
            if (ncpy) delete ncpy; 
            if (cpx) delete cpx; 
            if (cpy) delete cpy; 
        }

        if (csegx) delete[] csegx;
        if (csegy) delete[] csegy;
        if(nPerFct) delete[] nPerFct; // is allocated in child class

    }
    for (int i = 0; i < nFct; i++) delete fctCurr[i];
//    if(fctCurr) delete[] fctCurr;
    fctCurr.erase(fctCurr.begin(), fctCurr.end());

    zero_out();
//    hippo_Mark;
}

void CSplines2d::init(int _startindex, int _endindex, int _T, 
        double *_x, double *_y, int _nFct, traj_type *_id) 
{
//    hippo_Print(_id);
    VDynamics2d::init(_startindex, _endindex, _T, _x, _y, _nFct, _id); 
    initFctCurr(); // fct obj's have to be allocated before calling CSFunction::init

    if(!nPerFct) nPerFct= new int[nFct];
//    hippo_Print(nPerFct);
    for(int i=0; i<nFct; i++) nPerFct[i]= ncpx[i]*ncpy[i];

    CSFunction::init(_startindex, _endindex, _T, _nFct, _id);

    if(!csegx) csegx= new int[T];
    if(!csegy) csegy= new int[T];
    allocated.controls= true;


//    for(int i=0; i<nFct; i++) {
//        hippo_Print(cpx[i]);
//        hippo_Print(cpy[i]);
//    }

//    hippo_Mark;
//    for(int i=0; i<nFct; i++) cerr << ncpy[i] << "\t";
//    cerr << endl;
//    for(int i=0; i<ncpy[0]; i++)
//        cerr << i << '\t' << cpy[0][i] 
//        << '\t' << cpy[1][i] 
//        << '\t' << cpy[2][i] 
//        << '\t' << cpy[3][i] << '\n';


    int tid;
    for(int i=startindex; i <= endindex; i++) {
        tid= id? id[i] : 0;
        if (tid== -1) {
            csegx[i]= csegy[i]= -1;
        } else {
            // find segments for all values in advance to speed up filtering
            // algo
            fctCurr[tid]->findSegments(x[i], y[i], csegx[i], csegy[i]);
            hippo_Assert(csegx[i]>0 && csegx[i] < ncpx[tid]-1, "should not get here");
            hippo_Assert(csegy[i]>=0 && csegy[i] < ncpy[tid], "should not get here");
        }
    }

    // check wether mapping is correct
    if(mapping) {
        vector<mapItem>::iterator iter;
        for(iter= map.begin(); iter!= map.end(); ++iter) {
            tid= iter->t1;
            int i= 0;
            while( cpx[tid][i]< iter->min & i<ncpx[tid]) i++;
            while( cpx[tid][i]< iter->max & i<ncpx[tid]) {
                hippo_Assert(cpx[tid][i]== cpx[iter->t2][i], 
                        "mapped control points must be located at same values");
                i++;
            }
        }
    }

}   


double CSplines2d::getMinX() const 
{ 
//    hippo_Mark;
    if(!xmin.valid) {
        xmin= DBL_MAX;
        for(int traj= 0; traj < nFct; traj++) {
            if(fctCurr[traj]->getMinX() < xmin) xmin= fctCurr[traj]->getMinX();
//    hippo_Print(cpx[traj][1]);
        }
//    hippo_Print(xmin);
    }
    return xmin;
};

double CSplines2d::getMaxX() const 
{ 
//    hippo_Mark;
    if(!xmax.valid) {
        xmax= -DBL_MAX;
        for(int traj= 0; traj < nFct; traj++) {
            if(fctCurr[traj]->getMaxX() > xmax) xmax= fctCurr[traj]->getMaxX();
        }
    }
    return xmax;
};

double CSplines2d::getMinY() const 
{ 
//        hippo_Mark;
    if(!ymin.valid) {
        ymin= DBL_MAX;
        for(int traj= 0; traj < nFct; traj++) {
            if(fctCurr[traj]->getMinY() < ymin) ymin= fctCurr[traj]->getMinY();
        }
    }
    return ymin;
};

double CSplines2d::getMaxY() const 
{ 
//        hippo_Mark;
    if(!ymax.valid) {
        ymax= -DBL_MAX;
        for(int traj= 0; traj < nFct; traj++) {
            if(fctCurr[traj]->getMaxY() > ymax) ymax= fctCurr[traj]->getMaxY();
        }
    }
    return ymax;
};


void CSplines2d::initFctCurr()
{
    fctCurr.erase(fctCurr.begin(), fctCurr.end());
    fct.erase(fct.begin(), fct.end());
    for (int i = 0; i < nFct; i++) {
        CSplinesFct2d *tmp= periodicY ? new CSplinesFctLinCirc : new CSplinesFct2d;
        fctCurr.push_back(tmp);
        fct.push_back(tmp);
        // init. the current values 
        tmp->init(ncpx[i], ncpy[i], cpx[i], cpy[i]);
    }
}

void CSplines2d::getActiveIndices(int t, int &ntraj, traj_type trajs[], int ind[]) const
{
    ntraj= 0;
    trajs[0]= id? id[t] : 0;
    int tsegx= csegx[t];
    int tsegy= csegy[t];
    if (trajs[0] != -1 && tsegx!=-1 && tsegy!= -1) {
        ntraj++;
        fctCurr[trajs[0]]->getPartialIndices(ind, tsegx, tsegy);
        if(!mapping) return;
        vector<mapItem>::const_iterator i;
        for(i= map.begin(); i!= map.end(); ++i) {
            if(i->t1== trajs[0] & x[t]>i->min & x[t]<i->max) {
                // add mapping   
            ntraj++; 
            trajs[ntraj-1]= i->t2; 
            }
        }
    }
}

//double CSplines2d::eval(int t, double _x, double _y, traj_type traj) const {
//    if(t < startindex || t > endindex)
//        cerr << startindex << "\t" << endindex << "\t" << t << endl;

//    hippo_Assert(t >= startindex && t <= endindex, "t out of bounds");
//    hippo_Assert(T != -1, "Parameters was not properly initialized");

//    moveTo(t);
//    return eval(x,y,traj);
//}

double CSplines2d::eval(double _x, double _y, traj_type traj) const 
{
    if(traj < 0) return 1;

    hippo_Assert(traj >= 0 && traj < nFct, "traj out of bound");

    double result= op->op(fctCurr[traj]->eval(_x,_y));
    hippo_Assert(isfinite(result), "result is infinite");
    return result;
}

void  CSplines2d::eval(double _x, double _y, traj_type _id, Grad &g) const 
{
    if(_id < 0) { g.markInvalid(nUpdate); return; }

    double xtmp[4], ytmp[4];
    g.f= fctCurr[_id]->evalPartial(_x, _y, xtmp, ytmp);
//    g.f= fctCurr[id]->evalPartial(x, y, xtmp, ytmp, 
//            csegx[tCurr], csegy[tCurr]);
    int k=0;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            g.df[k++]= xtmp[i]*ytmp[j];
    op->op(nUpdate, g);
}

double CSplines2d::evalGrad(double _x, double _y,
        traj_type traj, double diff[]) const
{
//    cerr << ", cpy[1][0]= " << cpy[1][0] << ", cpy[1][1]= " << cpy[1][1] << endl;
    if(traj < 0) {
        for(int i=0; i<nUpdate; i++) diff[i]= 0;
        return 1;
    }
//    cerr << "t= " << t << "\tx= " << x << "\tp= " << y << "\n";
//    cerr << "\ttraj= " << traj << "\t segx= " << csegx[t] << "\t segy= " << csegy[t] << endl;
//    cerr << ", cpy[1][0]= " << cpy[1][0] << ", cpy[1][1]= " << cpy[1][1] << endl;
//    hippo_Print(t);
    double xtmp[4], ytmp[4];
    double result= fctCurr[traj]->evalPartial(_x, _y, xtmp, ytmp);
//    double result= fctCurr[traj]->evalPartial(x, y, xtmp, ytmp, 
//            csegx[tCurr], csegy[tCurr]);
    int k=0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
                diff[k]= xtmp[i]*ytmp[j];
            k++;
        }
    }
    op->op(result, nUpdate, diff);
    return result;
}

CSplines2d* CSplines2d::link_copy()
{
    CSplines2d *v= new CSplines2d;
    CSFunction::copyInit(v);
    CSplines2d::copyInit(v);
    v->store= store->alloc_new();
    v->store->init(startindex, endindex, v);

    return v;
}

void CSplines2d::copyInit(CSplines2d *v) const
{
    v->ncpx= ncpx;
    v->ncpy= ncpy;
    v->cpx= cpx;
    v->cpy= cpy;
    v->csegx= csegx;
    v->csegy= csegy;
    v->allocated.controls= false;

    v->periodicY= periodicY;

    v->initFctCurr();
}

void CSplines2d::setMexInput(const mxArray *in)
{
//    hippo_Mark;
    CSFunction::setMexInput(in);
    int tmp;
    // cp is matlab cell array with column vectors
    mxArray *cellPtr, *mxtmp;
//    double *dblPtr;

    cellPtr= mxGetField(in,0,"cpx");
    hippo_Assert(cellPtr, "no control points defined, errors from now on");
    if(mxIsCell(cellPtr)) {
        nFct = mxGetNumberOfElements(cellPtr);
        ncpx= new int[nFct];
        cpx= new double*[nFct];
        for (int i = 0; i < nFct; i++) {
            mxtmp = mxGetCell(cellPtr, i);
            ncpx[i] = mxGetN(mxtmp) * mxGetM(mxtmp);
            cpx[i]  = mxGetPr(mxtmp);
        }
    } else {
        nFct= 1;
        ncpx= new int;
        cpx= new double*;
        ncpx[0] = mxGetN(cellPtr) * mxGetM(cellPtr);
        cpx[0]  = mxGetPr(cellPtr);
    }
    allocated.controls= true;

    cellPtr= mxGetField(in,0,"cpy");
    hippo_Assert(cellPtr, "no control points defined");
    if(mxIsCell(cellPtr)) {
        tmp = mxGetNumberOfElements(cellPtr);
        hippo_Assert(nFct== tmp, "cpx and cpy have different number of traj");
        ncpy= new int[nFct];
        cpy= new double*[nFct];
        for (int i = 0; i < nFct; i++) {
            mxtmp = mxGetCell(cellPtr, i);
            ncpy[i] = mxGetN(mxtmp) * mxGetM(mxtmp);
            cpy[i]  = mxGetPr(mxtmp);
        }
    } else {
        hippo_Assert(nFct== 1, "cpx and cpy have different number of traj");
        ncpy= new int;
        cpy= new double*;
        ncpy[0] = mxGetN(cellPtr) * mxGetM(cellPtr);
        cpy[0]  = mxGetPr(cellPtr);
    }

    nTotal= 0;
    for (int i = 0; i < nFct; i++) {
        nTotal+= ncpx[i]*ncpy[i];
    }

    mxtmp= mxGetField(in, 0, "periodic");
    if(mxtmp) {
        char *tmpstring;
        MX_StringAllocAssign(mxtmp, tmpstring);
        if (strcmp(tmpstring,"y")==0) {
            periodicY= true;
        } else  MX_Error("do not understand 'periodic= %s'", tmpstring);
        mxFree(tmpstring);
    } else periodicY= false;

//    if(csegx) delete[] csegx;
//    mxtmp= mxGetField(in, 0, "csegx");
//    if(mxtmp) {
//        MX_FieldAssign(in, "csegx", dblPtr, tmp);
//        hippo_Assert(tmp== T, "csegx must have T elements");
//        csegx= new int[T];
//        for(int t=0; t < T; t++) csegx[t]= (int) round(dblPtr[t]);
//    } 

//    if(csegy) delete[] csegy;
//    mxtmp= mxGetField(in, 0, "csegy");
//    if(mxtmp) {
//        MX_FieldAssign(in, "csegy", dblPtr, tmp);
//        hippo_Assert(tmp== T, "csegy must have T elements");
//        csegy= new int[T];
//        for(int t=0; t < T; t++) csegy[t]= (int) round(dblPtr[t]);
//    }

    id= data? data->fieldID : 0;

}    

mxArray* CSplines2d::getMexOutput() 
{
    mxArray *out= CSFunction::getMexOutput();
    mxAddField(out, "name");
    mxAddField(out, "cpx");
    mxAddField(out, "cpy");
//    mxAddField(out, "csegx");
//    mxAddField(out, "csegy");

    mxArray *mxtmp= mxCreateString("CSplines2d");
    hippo_Assert(mxtmp, "couldn't allocate space for string");
    mxSetField(out, 0, "name", mxtmp);

    double *dblPtr;
    mxArray *cellPtr;
    cellPtr= mxCreateCellMatrix(1,nFct);
    for (int i = 0; i < nFct; i++) {
        MX_CreateDoubleAssign(mxtmp, 1,ncpx[i], dblPtr);
        memcpy(dblPtr, cpx[i], ncpx[i] * sizeof(double));
        mxSetCell(cellPtr, i, mxtmp);
    }
    mxSetField(out, 0, "cpx", cellPtr);

    cellPtr= mxCreateCellMatrix(1,nFct);
    for (int i = 0; i < nFct; i++) {
        MX_CreateDoubleAssign(mxtmp, 1,ncpy[i], dblPtr);
        memcpy(dblPtr, cpy[i], ncpy[i] * sizeof(double));
        mxSetCell(cellPtr, i, mxtmp);
    }
    mxSetField(out, 0, "cpy", cellPtr);

    if (periodicY) {
        mxArray *mxtmp= mxCreateString("y");
        hippo_Assert(mxtmp, "couldn't allocate space for string");
        mxAddField(out, "periodic");
        mxSetField(out, 0, "periodic", mxtmp);
    }
//    MX_CreateDoubleAssign(mxtmp, 1, T, dblPtr);
//    for(int t=0; t < T; t++) dblPtr[t]=  csegx[t];
//    mxSetField(out, 0, "csegx", mxtmp);

//    MX_CreateDoubleAssign(mxtmp, 1, T, dblPtr);
//    for(int t=0; t < T; t++) dblPtr[t]=  csegy[t];
//    mxSetField(out, 0, "csegy", mxtmp);


    return out;
}
