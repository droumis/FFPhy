#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "../aux/hippoIO.h"
#include "../aux/numerics.h"
#include "CSplinesFct.h"

using namespace std;
using namespace AFilter;

/* ************************************************************
   class  CSFct
 ************************************************************ */

CSFct::CSFct() 
    : param(0), allocated(false), nParam(0)
{};

CSFct::~CSFct() {
    dealloc();
};

void CSFct::allocParam(int n) {
    if(allocated) {
        if(n == nParam) return;
        else dealloc();
    }
    nParam= n;
    hippo_CAlloc(param,nParam,double);
    allocated= true;
//    hippo_Print(nParam);
};

void CSFct::dealloc() {
    if(allocated && param) {
        free(param);
        param= 0;
        allocated= false;
        nParam= 0;
    }
};



/* ************************************************************
   class  CSplinesFct
 ************************************************************ */

CSplinesFct::CSplinesFct() 
: ncpx(0), cpx(0)
{
}

CSplinesFct::CSplinesFct(int _ncpx, double* _cpx)
{
    init(_ncpx,  _cpx);
}

CSplinesFct::~CSplinesFct()
{
    release();
}

void CSplinesFct::init(int _ncpx, double* _cpx)
{
    hippo_Assert(_ncpx >= 4, "Cardinal splines need at least 4 control points");
    ncpx= _ncpx;
    cpx= _cpx;
    allocParam(ncpx);
}

void CSplinesFct::getPartialIndices(int a[], int tsegx)
{
    for (int i = -1; i < 3; i++) (*a++)= tsegx+i;
}

void CSplinesFct::release()
{
//    ncpx= 0;
//    if(allocatedTheta && theta) {
//        free(theta);
//    }
//    cpx= theta= NULL;    
}

int CSplinesFct::findSegment(double x) const
{
    hippo_Assert(param, "CSplinesFct was not initialized properly");

    hippo_Assert(x >= cpx[1], "x coordinate smaller than control points");
    hippo_Assert(x < cpx[ncpx-2], "x coordinate larger than control points");
#ifdef DEBUG
    if(x < cpx[1]){
        hippo_Print(x);
        for(int i=0; i < ncpx; i++) {
            hippo_Print(cpx[i]);
        }
    }
#endif        
//    if(x<0) {
//         hippo_Print(x);
//    }
    //     hippo_Print(cpx[1]);
    int csegx= 1;
    while (csegx < ncpx-1 && cpx[csegx+1] <= x) {
        csegx++;
    } 
//    hippo_Print(csegx);
    hippo_Assert(x >= cpx[csegx] && x < cpx[csegx+1], "wrong segment");
#ifdef DEBUG
    if(csegx >= ncpx-2){
        hippo_Print(x);
        for(int i=0; i < ncpx; i++) {
            hippo_Print(cpx[i]);
        }
    }
#endif        

    //
//    printf("x= %.3f, csegx= %3.d, [%.3f, %.3f]\n", x, csegx, cpx[csegx], cpx[csegx+1]);
    return csegx;
}

double CSplinesFct::eval(double x, int csegx) const
{
    double xtmp[4];
    return evalPartial(x, xtmp, csegx);
}

double CSplinesFct::evalPartial(double x, double partialx[], int csegx) const
{
    hippo_Assert(param, "CSplinesFct was not initialized properly");
    double x1, x2, x3;

    if (csegx== -1) csegx= findSegment(x);


    x1 =  (x - cpx[csegx]) / (cpx[csegx+1] - cpx[csegx]);
    x2 = x1 * x1;
    x3 = x2 * x1;

    hippo_Assert(x1>=0 && x1 <= 1, "x1 out of bounds");

    /* calculate x * Mc for the cardinal spline */
    partialx[0] = (-.5*x3 + x2 -.5*x1);
    partialx[1] = (1.5*x3 - 2.5*x2 + 1);
    partialx[2] = (-1.5*x3 + 2.0*x2 + .5*x1);
    partialx[3] = (.5*x3 - .5*x2); 

    //@@
//    printf("%.2f, %.2f, %.2f, %.2f\n", param[csegx-1], param[csegx], param[csegx+1], param[csegx+2]) ;

    // x * Mc * c
    return  param[csegx-1]*partialx[0] + param[csegx]*partialx[1] +
        param[csegx+1]*partialx[2] + param[csegx+2]*partialx[3];
}

double CSplinesFct::integral() const
{
    double integral= 0, tmp;
    for(int i= 1; i <= ncpx-3; i++) {
        tmp= 6.5*(param[i]+param[i+1]) - 0.5*(param[i-1]+param[i+2]);
        integral+= (cpx[i+1]-cpx[i])*tmp/12.0;
    }
    return integral;
}

/* ************************************************************
   class  CSplinesFct2d
 ************************************************************ */

CSplinesFct2d::CSplinesFct2d() 
: cpx(0), cpy(0), ncpx(0), ncpy(0)
{
}
CSplinesFct2d::CSplinesFct2d(int _ncpx, int _ncpy, double* _cpx, double* _cpy)
{
    init( _ncpx,  _ncpy,  _cpx,  _cpy);
}

CSplinesFct2d::~CSplinesFct2d()
{
}

void CSplinesFct2d::init(int _ncpx, int _ncpy, double* _cpx, double* _cpy)
{
    hippo_Assert(_ncpx >= 4, "Cardinal splines need at least 4 control points");
    hippo_Assert(_ncpy >= 4, "Cardinal splines need at least 4 control points");
    ncpx= _ncpx;
    ncpy= _ncpy;
    cpx= _cpx;
    cpy= _cpy;
    allocParam(ncpx*ncpy);
//    hippo_Print(nParam);
//    hippo_Mark;
//    for(int i=0; i<ncpy; i++)
//        cerr << i << '\t' << cpy[i] << '\n';
}

void CSplinesFct2d::release()
{
    ncpx= ncpy= 0;
    cpx= cpy= NULL;
}

void CSplinesFct2d::findSegments(double x, double y, int &csegx, int &csegy) const
{
    csegx = 1;
    csegy = 1;

    hippo_Assert(param, "CSplinesFct was not initialized properly");

    // check whether x and y values are within control points 
//    hippo_Assert(x >= cpx[1], "x coordinate smaller than control points");
//    hippo_Assert(y >= cpy[1], "y coordinate smaller than control points");
//    hippo_Assert(x < cpx[ncpx-2], "x coordinate larger than control points");
//    hippo_Assert(y < cpy[ncpy-2], "y coordinate larger than control points");
    if (!(x >= cpx[1])) {
        cerr << "x (" << x << ") below lowest control point (" <<
            cpx[1] << ")" << ", line " << __LINE__ << ", file " __FILE__ << endl;
        throw 1;
    } else if (!(x < cpx[ncpx-2])) {
        cerr << "x (" << x << ") above highest control point (" <<
            cpx[ncpx-2]<< ")" << ", line " << __LINE__ << ", file " __FILE__ << endl;
        throw 1;
    } else if (!(y >= cpy[1])) {
        cerr << "y (" << y << ") below lowest control point (" <<
            cpy[1] << ")" << ", line " << __LINE__ << ", file " __FILE__ << endl;
        throw 1;
    } else if (!(y < cpy[ncpy-2])) {
        cerr << "y (" << y << ") above highest control point (" <<
            cpy[ncpy-2] << ")" << ", line " << __LINE__ << ", file " __FILE__ << endl;
        throw 1;
    }


    // find correct segment
    while (csegx < ncpx-1 && cpx[csegx+1] <= x) csegx++;
    while (csegy < ncpy-1 && cpy[csegy+1] <= y) csegy++;

    hippo_Assert(x >= cpx[csegx] && x <= cpx[csegx+1], "wrong x segment");
    hippo_Assert(y >= cpy[csegy] && y <= cpy[csegy+1], "wrong y segment");
}

double CSplinesFct2d::eval(double x, double y, int csegx, int csegy) const
{
    double xtmp[4], ytmp[4];
    return evalPartial(x, y, xtmp, ytmp, csegx, csegy);
}

double CSplinesFct2d::evalPartial(double x, double y, 
        double partialx[], double partialy[],
        int csegx, int csegy) const
{
    hippo_Assert(param, "CSplinesFct was not initialized properly");

    double x1, x2, x3;
    double y1, y2, y3;

    //     hippo_Print(csegx);
    //     hippo_Print(csegy);
    /* find the segment of the spline corresponding to the current x value. */
    if (csegx==-1 || csegy==-1) 
        findSegments(x,y,csegx,csegy);
    //    hippo_Mark;

    x1 = (x - cpx[csegx]) / (cpx[csegx+1] - cpx[csegx]);
    x2 = x1 * x1;
    x3 = x2 * x1;

    y1 = (y - cpy[csegy]) / (cpy[csegy+1] - cpy[csegy]);
    y2= y1*y1;
    y3= y1*y2;

//    if (!(x1>=0 && x1 <= 1))
//        cerr << "x1 out of bounds (" << x1 << "), should be between 0 and 1" <<
//            ", line " << __LINE__ << ", file " __FILE__ << endl;
//    if (!(y1>=0 && y1 <= 1))
//        cerr << "y1 out of bounds (" << y1 << "), should be between 0 and 1" <<
//            ", line " << __LINE__ << ", file " __FILE__ << endl;
    hippo_Assert(x1>=0 && x1 <= 1, "x1 out of bounds");
    hippo_Assert(y1>=0 && y1 <= 1, "y1 out of bounds");
    //       hippo_Mark;

    /* calculate x * Mc for the cardinal spline */
    partialx[0] = (-.5*x3 + x2 -.5*x1);
    partialx[1] = (1.5*x3 - 2.5*x2 + 1);
    partialx[2] = (-1.5*x3 + 2.0*x2 + .5*x1);
    partialx[3] = (.5*x3 - .5*x2); 

    /* calculate y * Mc for the cardinal spline */
    partialy[0] = (-.5*y3 + y2 -.5*y1);
    partialy[1] = (1.5*y3 - 2.5*y2 + 1);
    partialy[2] = (-1.5*y3 + 2.0*y2 + .5*y1);
    partialy[3] = (.5*y3 - .5*y2); 

#ifdef DEBUG
    //@@
    // large partial
    bool large= false;
    for (int q=0; q < 4; q++) 
        if (partialx[q] > 1 || partialy[q] > 1) { large= true; break; }
    if (large) {
        for (int q=0; q < 4; q++) 
            printf("dx[%d]= %f,\tdy[%d]= %f\n", q, partialx[q], q, partialy[q]);
        printf("\n");
    }
#endif

    // calculate lambda(x,phase)= x * Mc * C * Mc * y
    double val= 0;
    for (int n=0; n < 4; n++) {
        for (int m=0; m < 4; m++) {
            val += partialx[n]*param[getLinear(csegx-1+n,csegy-1+m)]*partialy[m];
        }
    }
    hippo_Assert(isfinite(val), "val is infinite");
    return val;
}

void CSplinesFct2d::getPartialIndices(int a[], int csegx, int csegy)
{
    for (int i = -1; i < 3; i++)
        for (int j = -1; j < 3; j++)
            *a++= getLinear(csegx+i,csegy+j);
}


/* ************************************************************
   class  CSplinesFctLinCirc
 ************************************************************ */

/*
CSplinesFctLinCirc::CSplinesFctLinCirc() 
{
    // nothing to do
}

CSplinesFctLinCirc::CSplinesFctLinCirc(int _ncpx, int _ncpy, double* _cpx, double* _cpy)
    : CSplinesFct2d( _ncpx,  _ncpy,  _cpx,  _cpy)
{
    // nothing to do
}

CSplinesFctLinCirc::~CSplinesFctLinCirc()
{
    // nothing to do
}
*/

void CSplinesFctLinCirc::findSegments(double x, double y, int &csegx, int &csegy) const
{
    csegx = 1;
    csegy = 0;

    hippo_Assert(param, "CSplinesFct was not initialized properly");

    // check whether x and y values are within control points 
    if (!(x >= cpx[0])) {
        cerr << "x (" << x << ") below lowest control point (" <<
            cpx[0] << ")" << ", line " << __LINE__ << ", file " __FILE__ << endl;
        throw 1;
    } else if (!(x < cpx[ncpx-2])) {
        cerr << "x (" << x << ") above highest control point (" <<
            cpx[ncpx-2]<< ")" << ", line " << __LINE__ << ", file " __FILE__ << endl;
        throw 1;
    } else if (!(y >= cpy[0])) {
        cerr << "y (" << y << ") below lowest control point (" <<
            cpy[0] << ")" << ", line " << __LINE__ << ", file " __FILE__ << endl;
        throw 1;
    } else if (!(y < cpy[ncpy-1])) {
        cerr << "y (" << y << ") above highest control point (" <<
            cpy[ncpy-1] << ")" << ", line " << __LINE__ << ", file " __FILE__ << endl;
        throw 1;
    }


    // find correct segment
    while (csegx < ncpx-1 && cpx[csegx+1] <= x) csegx++;
    while (csegy < ncpy-1 && cpy[csegy+1] <= y) csegy++;

//    if(x < cpx[csegx] || x > cpx[csegx+1]) {
//        hippo_Print(y);
//        hippo_Print(csegy);
//        hippo_Print(cpy[csegy]);
//        hippo_Print(cpy[csegy+1]);
//    }
    hippo_Assert(x >= cpx[csegx] && x <= cpx[csegx+1], "wrong x segment");
    hippo_Assert(y >= cpy[csegy] && y <= cpy[csegy+1], "wrong y segment");

    hippo_Assert(csegx >= 1 && csegx <= ncpx-3, "x index out of bounds");
    hippo_Assert(csegy >= 0 && csegy <= ncpy-2, "y index out of bounds");
}

/*
double CSplinesFctLinCirc::evalPartial(double x, double y, 
        double partialx[], double partialy[],
        int csegx, int csegy) const
{
//    cerr << "\t\tx= " << x << "\tp= " << y << "\n";
//    cerr << "\t\t segx= " << csegx << "\t segy= " << csegy << endl;
//    cerr << ", seg= " << 0 << ", cpy[0]= " << cpy[0] << ", cpy[0+1]= " << cpy[0+1] << endl;

    hippo_Assert(param, "CSplinesFct was not initialized properly");

    double x1, x2, x3;
    double y1, y2, y3;

//     find the segment of the spline corresponding to the current x value. 
    if (csegx==-1 || csegy==-1)  {
        findSegments(x,y,csegx,csegy);
//        hippo_Print(x);
//        hippo_Print(y);
    }

    x1 = (x - cpx[csegx]) / (cpx[csegx+1] - cpx[csegx]);
    x2 = x1 * x1;
    x3 = x2 * x1;

    y1 = (y - cpy[csegy]) / (cpy[csegy+1] - cpy[csegy]);
    y2= y1*y1;
    y3= y1*y2;

#ifdef DEBUG
    if (!(x1>=0 && x1 <= 1)) {
        cerr << "x1 out of bounds (" << x1 << "), should be between 0 and 1" <<
            ", line " << __LINE__ << ", file " __FILE__ << endl;
        cerr << "x= " << x << ", seg= " << csegx << ", cpx[seg]= " << cpx[csegx] << ", cpx[seg+1]= " << cpx[csegx+1] << endl;
        throw(1);
    }
    if (!(y1>=0 && y1 <= 1)) {
        cerr << "y1 out of bounds (" << y1 << "), should be between 0 and 1" <<
            ", line " << __LINE__ << ", file " __FILE__ << endl;
        cerr << "y= " << y << ", seg= " << csegy << ", cpy[seg]= " << cpy[csegy] << ", cpy[seg+1]= " << cpy[csegy+1] << endl;
    hippo_Print(cpy);
        throw(1);
    }
#endif

//     calculate x * Mc for the cardinal spline 
    partialx[0] = (-.5*x3 + x2 -.5*x1);
    partialx[1] = (1.5*x3 - 2.5*x2 + 1);
    partialx[2] = (-1.5*x3 + 2.0*x2 + .5*x1);
    partialx[3] = (.5*x3 - .5*x2); 

//     calculate y * Mc for the cardinal spline 
    partialy[0] = (-.5*y3 + y2 -.5*y1);
    partialy[1] = (1.5*y3 - 2.5*y2 + 1);
    partialy[2] = (-1.5*y3 + 2.0*y2 + .5*y1);
    partialy[3] = (.5*y3 - .5*y2); 

    //@@
    // large partial
    bool large= false;
    for (int q=0; q < 4; q++) 
        if (partialx[q] > 1 || partialy[q] > 1) { large= true; break; }
    if (large) {
        for (int q=0; q < 4; q++) 
            printf("dx[%d]= %f,\tdy[%d]= %f\n", q, partialx[q], q, partialy[q]);
        printf("\n");
    }

    // calculate lambda(x,phase)= x * Mc * C * Mc * p
    double val= 0;
    for (int n=0; n < 4; n++) {
        for (int m=0; m < 4; m++) {
            int tmp= mod((csegy-1+m), ncpy-1);
            hippo_Assert(tmp >= 0 && tmp <= ncpy-1, "y index out of bounds");

            val += partialx[n]*param[getLinear(csegx-1+n, csegy-1+m)]*partialy[m];
        }
    }
    hippo_Assert(isfinite(val), "val is infinite");
    return val;
}

*/

//ostream& operator << (ostream &os, CSplinesFctLinCirc &cs)
//{
//    os << (CSplinesFct2d) cs;
//    return os;
//}
