#ifndef CSPLINES2D_H
#define CSPLINES2D_H

/* $Id: CSplines2d.h,v 1.7 2008/08/24 19:47:03 chengs Exp $
   File name   : /bach/AdaptFilter/model/CSplines2d.h
   authors     : Sen Cheng
   created     : Fri 08 Oct 2004 04:17:25 PM PDT
  
 */

#include "../aux/defs.h"
#include "VDynamics2d.h"
#include "CSStorage.h"

class TFilter;

/* ************************************************************
                          class CSplines2d
   ************************************************************ */
namespace AFilter {
    
//class CardSplines2d;

class CSplines2d : public VDynamics2d, public CSFunction {
private: // data structures
    /* at each timestep, store only 16 control points that change
        0:xyval11  1:xyval12  2:xyval13  3:xyval14
        4:xyval21  5:xyval22  6:xyval23  7:xyval24
        8:xyval31  9:xyval32 10:xyval33 11:xyval34
       12:xyval41 13:xyval42 14:xyval43 15:xyval44
		   */

    vector<CSplinesFct2d *> fctCurr;
    int     *ncpx, *ncpy;   // the numbers of spatial control points
    double  **cpx, **cpy; // the lists of spatial control points. cpx[traj][c]

    int *csegx, *csegy;

    bool periodicY;

    struct {
        bool controls;
    } allocated;

protected:
//    void alloc(int NTraj, int NSteps);
    void dealloc();
    void initFctCurr();

    void zero_out();    
    void copyInit(CSplines2d *v) const;
public:
    CSplines2d();
    CSplines2d(CSplines2d &);
    virtual ~CSplines2d();

    virtual CSplines2d* alloc_new() const { return new CSplines2d; }
    virtual CSplines2d* link_copy();

    virtual const char* getName() const { return "CSplines2d"; };

    virtual void init(int startindex, int endindex, int _T, double *x, double *y, int nid, 
            traj_type *id);    


    virtual double getMinX() const;
    virtual double getMaxX() const;
    virtual double getMinY() const;
    virtual double getMaxY() const;


//    virtual int getNParam() const { return nUpdate; };
//    virtual int getNAllParam() const { return  totalncp; }; 

    virtual void getActiveIndices(int t, int &ntraj, traj_type trajs[], 
            int ind[]) const;

//    virtual traj_type getParam(int t, double [], int a[]) const;
//    virtual void getAllParam(int t, double z[]) const;

//    virtual double eval(int t, double x, double y, traj_type traj) const ;

    virtual double eval(double x, double y, traj_type traj) const ;

    virtual void eval(double x, double y, traj_type traj, Grad &g) const;

    /** calculate f(x) and d(log f)/dz, where z are the parameters
     */
    virtual double evalGrad(double x, double y,
            traj_type traj, double partial[]) const ;


//    virtual void setInitial(double x);
//    virtual void setInitialToCurrent();
//    virtual void setInitialEst(int t= -1) const;
//    virtual void setInitialFrom(VDynamics2d *);

//    virtual void updateDynamics(int tStart, int tEnd, double diff[]);

//    virtual bool hasConverged(VDynamics *);
//    virtual CSplines2d* allocForConverge();
//    virtual void copyForConverge(VDynamics *);

    virtual void setMexInput(const mxArray *in);
    virtual mxArray* getMexOutput();

//    friend class CardSplines2d;
};

}; //namespace AFilter

#endif   // CSPLINES2D_H
