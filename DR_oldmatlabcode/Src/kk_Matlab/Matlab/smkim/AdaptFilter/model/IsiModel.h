#ifndef IsiModel_H
#define IsiModel_H

/* $Id: IsiModel.h,v 1.2 2008/09/01 18:24:22 chengs Exp $
   
   Sen Cheng,  Thu Aug 21 09:33:53 PDT 2008
*/

#include "../model/AdaptModel.h"
#include "../model/VDynamics1d.h"
#include "../model/VDynamics2d.h"
#include "../aux/defs.h"

//struct mxArray_tag;
//typedef struct mxArray_tag mxArray;

/* ************************************************************
                          class IsiModel
   ************************************************************ */

namespace AFilter {

/** Abstract base class for models with a "spatial" and an isi component
 */
class IsiModel : public AdaptModel {
protected:
    VDynamics *xp;
    VDynamics1d *isi;
    int nxp, nisi;
    int totalnxp;

    double max_isi;
    traj_type *isiID;

    struct {
        bool isiID;
    } allocated;
    
protected:
    IsiModel(TData *_data=0);

    virtual void init_spatial(int n) =0;
public:
    virtual ~IsiModel();

    virtual const char* getName() const =0;

    inline double getMaxIsi() const { return max_isi; }
    inline VDynamics1d* getIsiModel() const { return isi; }
    inline VDynamics* getSpatialModel() const { return xp; }

    VDynamics1d* get1dSpatialModel() const; 
    VDynamics2d* get2dSpatialModel() const; 

    virtual void init();
//    virtual void update(int t);

    virtual void moveTo(int t) const { xp->moveTo(t); isi->moveTo(t); }
    virtual double eval(int timeindex) const;

    /** calculate f(x), log(x), df/dz, d(log f)/dz, d^2f/dz^2, d^2(log f)/d^2z 
    */
    virtual void eval(int t, Grad &g) const;

    /** calculate f(x), log(x), df/dz, d(log f)/dz, d^2f/dz^2, d^2(log f)/d^2z 
    */
    virtual void evalNoMove(int timeindex, Grad &g) const;

    /** calculate f(x) and df/dz, where z are the parameters */
    virtual double evalGrad(int timeindex, double partial[]) const;

    /// For filtering algorithm. Evaluate without moving.
    virtual double evalNoMove(int timeindex, double partial[]) const;

//    /** calculate log f(x) and d(log f)/dz, where z are the parameters */
//    virtual double evalLogGrad(int timeindex, double partial[]) const;

    virtual int getNParam() const { return nxp+nisi; }
    virtual int getNAllParam() const { return xp->getNAllParam()+isi->getNAllParam();}; 

    virtual bool getParam(int t, double z[], int a[]) const;
    virtual void getAllParam(int t, double p[]) const;

    virtual void setParam(int t, double z[]);
    virtual void setAllParam(int t, double x) {
        xp->setAllParam(t,x); isi->setAllParam(t,x);
    }
    virtual void setAllParam(int t, double x[]) {
        xp->setAllParam(t,x); isi->setAllParam(t,x+totalnxp);
    }

//    virtual void setInitialTo(double x) {xp->setInitial(x); isi->setInitial(x); }
    virtual void setInitialEst(int t);
//    virtual void setInitialEst();
//    virtual void setInitialFrom(AdaptModel *);
//    virtual void setInitialToCurrent()
//    {xp->setInitialToCurrent(); isi->setInitialToCurrent(); }
//    virtual void resetToInitial();

    inline virtual void updateModel(int tStart, int tEnd, double diff[]) {
        updateModel(tStart, tEnd, diff, -1);
    }

    // fixing= {-1, 0, 1} for fixing no component, the first (spatial) component
    // and the second (isi) component
    virtual void updateModel(int tStart, int tEnd, double diff[], int fixing);

//    virtual void setNoReverseUpdate() { xp->setNoReverseUpdate();  isi->setNoReverseUpdate(); }
//    virtual void setReverseUpdate() { xp->setReverseUpdate();  isi->setReverseUpdate(); }


    virtual void allocForConvergence();
    virtual void copyForConverge();
    virtual bool hasConverged();
    virtual void averageWithCopy() {
        hippo_Assert(convCopy, "must allocate copies before copying for convergence");
        xp->averageWithCopy();
        isi->averageWithCopy();
    }

    virtual IsiModel* link_copy();

    virtual void setMexInput(const mxArray *in);
    virtual mxArray* getMexOutput();

/**** for generating spike trains ******/
    virtual double evalAtIsi(int timeindex, double isi);
};

} // namespace AFilter

#endif   // IsiModel_H
