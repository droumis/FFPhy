/* $Id: VDynamics.h,v 1.10 2008/09/01 18:24:22 chengs Exp $
   File name   : model/VDynamics.h
   authors     : Sen Cheng
   created     : Sun 17 Oct 2004 10:48:58 AM PDT
  
 */
#ifndef VDYNAMICS_H
#define VDYNAMICS_H

#include "../aux/defs.h"
#include "../aux/numerics.h"
#include "../aux/mexAux.h"
#include "../aux/hippoIO.h"
#include "../aux/TData.h"
#include "../aux/Grad.h"

namespace AFilter {

/** Abstract class for dynamics, i.e. function changing in time.
 * 
 * @author Sen Cheng
 */
class VDynamics 
{
private:
protected:
    TData *data;
    
    int startindex, endindex;
    int T;

    mutable int tCurr;
//    mutable int tIni;
    
    int  nFct;	      // number of functions
    traj_type *id;            // vector of function id's

//    mutable lazyVar<double> zmin, zmax;
//    bool reverseUpdate;
protected:
    void copyInit(VDynamics *v) const {
        v->startindex= startindex; v->endindex= endindex; v->T= T; v->nFct= nFct; v->id= id; 
    }
public:
//    VDynamics() : startindex(0), endindex(0), T(0), tCurr(-1), tIni(-1), nFct(0), id(0), reverseUpdate(true) {};
    VDynamics() : startindex(0), endindex(0), T(0), nFct(0), id(0){};
    VDynamics(VDynamics &) {};
    virtual ~VDynamics() {};

    virtual const char* getName() const =0;

    virtual int getDim() const =0;
    
    inline int getStartindex() const { return startindex; }
    inline int getEndindex() const { return endindex; }
    inline int getT() { return T; }
    inline double getMinTime() const { return data->timearray[startindex]; };
    inline double getMaxTime() const { return data->timearray[endindex]; };

    virtual int getNFct() const { return nFct; }
    virtual int getNParam() const {hippo_Virtual; return 0; };
    virtual int getNAllParam() const {hippo_Virtual; return 0; };

    virtual bool getParam(int t, double [], int a[]) const {hippo_Virtual; return false; };
    virtual void getAllParam(int t, double []) const {hippo_Virtual; };
    virtual void setAllParam(int t, double []) {hippo_Virtual; };
    virtual void setAllParam(int t, double ) {hippo_Virtual; };
    virtual void setValues(int t, double ) {hippo_Virtual; };

    virtual double getMinZ() const { hippo_Virtual; } ;
    virtual double getMaxZ() const { hippo_Virtual; } ;

    virtual void init(int _startindex, int _endindex, int _T, int _nFct, 
                traj_type* _id) 
    {
        startindex= _startindex; endindex= _endindex; T= _T;
        hippo_Assert(_nFct== nFct, "model must be defined with same number of trajectories as are in data");
        nFct= _nFct; 
        id= _id; 
    }
//    virtual void setInitial(double) {hippo_Virtual;}
//    virtual void setInitialToCurrent() {hippo_Virtual;}
//    virtual void setInitialEst(int t= -1) const
// {hippo_Virtual;};
    
    virtual void moveTo(int timeindex) const =0;

    virtual double evalAtTime(int timeindex) const =0; 

    virtual void evalAtTime(int timeindex, Grad &results) const =0;

    virtual double evalGradAtTime(int timeindex, double partial[]) const =0;

    virtual double evalNoMove(int timeindex, double partial[]) const =0;

    virtual void evalNoMove(int t, Grad &results) const =0;

    virtual void updateDynamics(int tStart, int tEnd, double diff[]) =0;

//    virtual void setNoReverseUpdate() { reverseUpdate= false; };
//    virtual void setReverseUpdate() { reverseUpdate= true; };

    virtual double getIntegral(int timestep, traj_type traj) const 
    { hippo_Virtual; return 0;}
    
    virtual void allocForConverge() { hippo_Virtual; }
    virtual void copyForConverge() { hippo_Virtual; }
    virtual bool hasConverged() { hippo_Virtual; return false; }
    virtual void averageWithCopy() { hippo_Virtual; }
    virtual VDynamics* link_copy() { hippo_Virtual; return 0; };
    
    TData* getData() const { return data; };
    void setData(TData *d) { data= d; };
    void setId(traj_type *i) { id= i; };
    
    virtual void setMexInput(const mxArray *in) { 
        MX_FieldScalarDefault(in, "startindex", startindex, 0, int);
        MX_FieldScalarDefault(in, "endindex", endindex, 0, int);
    };
    virtual mxArray* getMexOutput() { 
        const char *fields[]= { "startindex", "endindex"};
        mxArray *out= mxCreateStructMatrix(1,1,2,fields);
        mxSetField(out,0,"startindex",mxCreateDoubleScalar(startindex));
        mxSetField(out,0,"endindex",mxCreateDoubleScalar(endindex));
        return out;
    };

}; // class VDynamics

}; // namespace AFilter 

#endif   // VDYNAMICS_H
