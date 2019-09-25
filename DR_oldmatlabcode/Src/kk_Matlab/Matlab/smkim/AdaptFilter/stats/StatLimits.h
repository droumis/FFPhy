/* $Id: StatLimits.h,v 1.4 2008/08/19 21:51:35 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#ifndef STATLIMITS_H
#define STATLIMITS_H

#include "../aux/defs.h"
#include "../aux/hippoIO.h"
#include "../aux/mexAux.h"

/* ************************************************************
                          class StatLimits
   ************************************************************ */

namespace AFilter {


/** Handles selection criteria for analyzing data and/or estimated models.
 *
 * Each analysis can be run several times (nana) with different criteria.
 * Each selection criterion can consist of several variables (nvar[ana]) with
 * either a single number or a range of valid numbers.
 * In addition, each analysis can be run for multiple values (nsteps[ana]) of an
 * independent variable, e.g., time.
 * 
 * @author Sen Cheng
 * @version 
 */
class StatLimits 
{
private:
    bool valid;     // whether object was initialized and is good to use
    mwSize nana;       // number of analyses to compute

    // For the constraints
    int *nvar;      // number of constraints for each analysis, nvar[ana]
    const char ***varnames; // (*name)[ana][var]
    double **varmin; // min[ana][var]
    double **varmax; // max[ana][var]

    // for the independent variable
    int *nsteps;    // number of time steps for each analysis, nsteps[ana] 
    const char **xname;
    double **xmin;  // min[ana][timestep]
    double **xmax;  // max[ana][timestep]

    // special pointers (shortcuts) for steps
    double **stepBegin;
    double **stepEnd;
private:
    //void allocate(int n);
    void dealloc();
    void zeroOut();
    int name2index(const char varname[], int ana) const;
public:
    StatLimits(const mxArray *in= 0);
    StatLimits(StatLimits &);
    ~StatLimits();

    int getNAnalyses() const { 
        hippo_Assert(valid, "initialize StatLimit object properly");
        return nana; 
    }
    int getNSteps(int ana) const { 
        hippo_Assert(valid, "initialize StatLimit object properly");
        return nsteps[ana]; 
    }
    int* getNSteps() const { 
        hippo_Assert(valid, "initialize StatLimit object properly");
        return nsteps;
    }

    bool isConstraint(const char varname[], int ana) const;

    // Return start time (sec) for analysis referenced ana and step
    double getXMin(mwSize ana, int step) const;

    // Return end time (sec) for analysis referenced ana and step
    double getXMax(mwSize ana, int step) const;

    // if a variable has only one value, it will be assigned to max[var],
    // min[var]= 0
    double getVarMin(const char varname[], int ana) const;
    double getVarMax(const char varname[], int ana) const;

    double getStepBegin(int ana, int step) const {
        hippo_Assert(valid, "initialize StatLimit object properly");
        return stepBegin[ana][step];
    }
    double getStepEnd(int ana, int step) const {
        hippo_Assert(valid, "initialize StatLimit object properly");
        return stepEnd[ana][step];
    }

    void setMexInput(const mxArray *in);
    mxArray* getMexOutput();
};

} // namespace AFilter

#endif   // STATLIMITS_H
