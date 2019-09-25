/* $Id: VStat.h,v 1.4 2008/08/20 21:47:14 chengs Exp $
   
   Sen Cheng, Fri Sep  3 12:58:58 PDT 2004

   program description
*/

#ifndef VSTAT_H
#define VSTAT_H

#include "StatLimits.h"
#include "../aux/defs.h"
#include "../aux/TData.h"
#include "../model/AdaptModel.h"

/* ************************************************************
                          class VStat
   ************************************************************ */

namespace AFilter {

/** Abstract base class for all statistics.
   @author Sen Cheng
 */
class VStat 
{
private:
    // storage for requested statistics, even if stat is n-dim.
    double **stat;

    char *label;
protected:
    // # of analyses
    int nana;
    
    // # of outputs
    int ndim;
    
    // # of steps, at which to compute stats
    int *nsteps;

    /// pointer to object with all experimental data
    const TData* data;

    const AdaptModel *model;

    /// 
    const StatLimits *lim;

private:
    void dealloc();
protected:
    VStat(const TData *d=0, const AdaptModel *m=0, const StatLimits *lim= 0);

    int getNDim() const { return ndim;}
    int getNSteps(int ana) const { return nsteps[ana]; }

    double get(int ana, int dim, int steps) const;
    void set(int ana, int dim, int steps, double x);
    
    virtual void init(int ana, int ndim, int *nsteps= 0);

public:
    virtual ~VStat();

    virtual const char* getName() const =0;

    virtual void init_obj() =0;

    virtual void run() =0; 
    
    virtual void setMexInput(const mxArray *in);
    virtual mxArray* getMexOutput();
    
    friend mxArray* collectMexOutput(int n, VStat *stats[]);
};

mxArray* collectMexOutput(int n, VStat *stats[]);

} // namespace AFilter

#endif   // VSTAT_H
