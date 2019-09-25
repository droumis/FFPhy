/*
 * MATLAB Compiler: 3.0
 * Date: Wed Dec  6 10:10:35 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "debugVis" 
 */
#include "debugVis.h"
#include "libmatlbm.h"
#include "libmmfile.h"
#include "visModel_mex_interface.h"

extern mxArray * adaptest;
extern mxArray * cmap;

static mxChar _array1_[19] = { '/', 'h', 'o', 'm', 'e', '/', 'c', 'h', 'e', 'n',
                               'g', 's', '/', 't', 'h', 'e', 't', 'a', '/' };
static mxArray * _mxarray0_;

static mxChar _array3_[16] = { '/', 'd', 'a', 't', 'a', '2', '/', 'b',
                               'e', 'h', 'a', 'v', 'd', 'a', 't', 'a' };
static mxArray * _mxarray2_;

static mxChar _array5_[4] = { '%', '.', '2', 'd' };
static mxArray * _mxarray4_;

static mxChar _array7_[16] = { '/', 'd', 'a', 't', 'a', '2', '/', 's',
                               'p', 'i', 'k', 'e', 'd', 'a', 't', 'a' };
static mxArray * _mxarray6_;

static mxChar _array9_[8] = { 'a', 'd', 'a', 'p', 't', 'e', 's', 't' };
static mxArray * _mxarray8_;

static mxChar _array11_[34] = { '/', 'b', 'a', 'c', 'h', '/', 'A', 'd', 'a',
                                'p', 't', 'F', 'i', 'l', 't', 'e', 'r', '/',
                                'n', 'e', 'w', 'v', 'i', 's', '/', 'c', 'o',
                                'l', 'o', 'r', 'm', 'a', 'p', '2' };
static mxArray * _mxarray10_;
static mxArray * _mxarray12_;

void InitializeModule_debugVis(void) {
    _mxarray0_ = mclInitializeString(19, _array1_);
    _mxarray2_ = mclInitializeString(16, _array3_);
    _mxarray4_ = mclInitializeString(4, _array5_);
    _mxarray6_ = mclInitializeString(16, _array7_);
    _mxarray8_ = mclInitializeString(8, _array9_);
    _mxarray10_ = mclInitializeString(34, _array11_);
    _mxarray12_ = mclInitializeDouble(500.0);
}

void TerminateModule_debugVis(void) {
    mxDestroyArray(_mxarray12_);
    mxDestroyArray(_mxarray10_);
    mxDestroyArray(_mxarray8_);
    mxDestroyArray(_mxarray6_);
    mxDestroyArray(_mxarray4_);
    mxDestroyArray(_mxarray2_);
    mxDestroyArray(_mxarray0_);
}

static void MdebugVis(mxArray * rat,
                      mxArray * din,
                      mxArray * ein,
                      mxArray * tin,
                      mxArray * cin,
                      mxArray * ratio);

_mexLocalFunctionTable _local_function_table_debugVis
  = { 0, (mexFunctionTableEntry *)NULL };

/*
 * The function "mlfDebugVis" contains the normal interface for the "debugVis"
 * M-function from file "/home/chengs/theta/code/debugVis.m" (lines 1-27). This
 * function processes any input arguments and passes them to the implementation
 * version of the function, appearing above.
 */
void mlfDebugVis(mxArray * rat,
                 mxArray * din,
                 mxArray * ein,
                 mxArray * tin,
                 mxArray * cin,
                 mxArray * ratio) {
    mlfEnterNewContext(0, 6, rat, din, ein, tin, cin, ratio);
    MdebugVis(rat, din, ein, tin, cin, ratio);
    mlfRestorePreviousContext(0, 6, rat, din, ein, tin, cin, ratio);
}

/*
 * The function "mlxDebugVis" contains the feval interface for the "debugVis"
 * M-function from file "/home/chengs/theta/code/debugVis.m" (lines 1-27). The
 * feval function calls the implementation version of debugVis through this
 * function. This function processes any input arguments and passes them to the
 * implementation version of the function, appearing above.
 */
void mlxDebugVis(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]) {
    mxArray * mprhs[6];
    int i;
    if (nlhs > 0) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: debugVis Line: 1 Column:"
            " 1 The function \"debugVis\" was called with m"
            "ore than the declared number of outputs (0)."),
          NULL);
    }
    if (nrhs > 6) {
        mlfError(
          mxCreateString(
            "Run-time Error: File: debugVis Line: 1 Column:"
            " 1 The function \"debugVis\" was called with m"
            "ore than the declared number of inputs (6)."),
          NULL);
    }
    for (i = 0; i < 6 && i < nrhs; ++i) {
        mprhs[i] = prhs[i];
    }
    for (; i < 6; ++i) {
        mprhs[i] = NULL;
    }
    mlfEnterNewContext(
      0, 6, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4], mprhs[5]);
    MdebugVis(mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4], mprhs[5]);
    mlfRestorePreviousContext(
      0, 6, mprhs[0], mprhs[1], mprhs[2], mprhs[3], mprhs[4], mprhs[5]);
}

/*
 * The function "MdebugVis" is the implementation version of the "debugVis"
 * M-function from file "/home/chengs/theta/code/debugVis.m" (lines 1-27). It
 * contains the actual compiled code for that M-function. It is a static
 * function and must only be called from one of the interface functions,
 * appearing below.
 */
/*
 * function debugVis(rat, din, ein, tin, cin, ratio)
 */
static void MdebugVis(mxArray * rat,
                      mxArray * din,
                      mxArray * ein,
                      mxArray * tin,
                      mxArray * cin,
                      mxArray * ratio) {
    mexLocalFunctionTable save_local_function_table_
      = mclSetCurrentLocalFunctionTable(&_local_function_table_debugVis);
    mxArray * opt = NULL;
    mxArray * model = NULL;
    mxArray * spikedata = NULL;
    mxArray * data = NULL;
    mxArray * behavdata = NULL;
    mxArray * ans = NULL;
    mxArray * c = NULL;
    mxArray * t = NULL;
    mxArray * e = NULL;
    mxArray * d = NULL;
    mclCopyArray(&rat);
    mclCopyArray(&din);
    mclCopyArray(&ein);
    mclCopyArray(&tin);
    mclCopyArray(&cin);
    mclCopyArray(&ratio);
    /*
     * % ratio: movietime/ realtime  [only works for 2ms timesteps]
     * % 
     * 
     * d=str2num(din); e= str2num(ein); t= str2num(tin); c= str2num(cin);
     */
    mlfAssign(&d, mlfStr2num(NULL, mclVa(din, "din")));
    mlfAssign(&e, mlfStr2num(NULL, mclVa(ein, "ein")));
    mlfAssign(&t, mlfStr2num(NULL, mclVa(tin, "tin")));
    mlfAssign(&c, mlfStr2num(NULL, mclVa(cin, "cin")));
    /*
     * 
     * load(['/home/chengs/theta/' rat '/data2/behavdata' sprintf('%.2d', d)]);
     */
    {
        mxArray * name_
          = mclInitialize(
              mlfHorzcat(
                _mxarray0_,
                mclVa(rat, "rat"),
                _mxarray2_,
                mlfSprintf(NULL, _mxarray4_, mclVv(d, "d"), NULL),
                NULL));
        mclLoadConditional(
          name_,
          "adaptest",
          mclPrepareGlobal(&adaptest),
          "ans",
          &ans,
          "behavdata",
          &behavdata,
          "c",
          &c,
          "cin",
          &cin,
          "cmap",
          mclPrepareGlobal(&cmap),
          "d",
          &d,
          "data",
          &data,
          "din",
          &din,
          "e",
          &e,
          "ein",
          &ein,
          "model",
          &model,
          "opt",
          &opt,
          "rat",
          &rat,
          "ratio",
          &ratio,
          NULL);
        mclLoadConditional(
          name_, "spikedata", &spikedata, "t", &t, "tin", &tin, NULL);
        mxDestroyArray(name_);
    }
    /*
     * load(['/home/chengs/theta/' rat '/data2/spikedata' sprintf('%.2d', d)]);
     */
    {
        mxArray * name_
          = mclInitialize(
              mlfHorzcat(
                _mxarray0_,
                mclVa(rat, "rat"),
                _mxarray6_,
                mlfSprintf(NULL, _mxarray4_, mclVv(d, "d"), NULL),
                NULL));
        mclLoadConditional(
          name_,
          "adaptest",
          mclPrepareGlobal(&adaptest),
          "ans",
          &ans,
          "behavdata",
          &behavdata,
          "c",
          &c,
          "cin",
          &cin,
          "cmap",
          mclPrepareGlobal(&cmap),
          "d",
          &d,
          "data",
          &data,
          "din",
          &din,
          "e",
          &e,
          "ein",
          &ein,
          "model",
          &model,
          "opt",
          &opt,
          "rat",
          &rat,
          "ratio",
          &ratio,
          NULL);
        mclLoadConditional(
          name_, "spikedata", &spikedata, "t", &t, "tin", &tin, NULL);
        mxDestroyArray(name_);
    }
    /*
     * 
     * 
     * data= behavdata{d}{e};
     */
    mlfAssign(
      &data,
      mlfIndexRef(
        mclVv(behavdata, "behavdata"), "{?}{?}", mclVv(d, "d"), mclVv(e, "e")));
    /*
     * data.spiketimes= spikedata{d}{e}{t}{c}.time;
     */
    mlfIndexAssign(
      &data,
      ".spiketimes",
      mlfIndexRef(
        mclVv(spikedata, "spikedata"),
        "{?}{?}{?}{?}.time",
        mclVv(d, "d"),
        mclVv(e, "e"),
        mclVv(t, "t"),
        mclVv(c, "c")));
    /*
     * data.spikeindex= spikedata{d}{e}{t}{c}.index;
     */
    mlfIndexAssign(
      &data,
      ".spikeindex",
      mlfIndexRef(
        mclVv(spikedata, "spikedata"),
        "{?}{?}{?}{?}.index",
        mclVv(d, "d"),
        mclVv(e, "e"),
        mclVv(t, "t"),
        mclVv(c, "c")));
    /*
     * 
     * 
     * load(['adaptest' sprintf('%.2d', d)]);
     */
    {
        mxArray * name_
          = mclInitialize(
              mlfHorzcat(
                _mxarray8_,
                mlfSprintf(NULL, _mxarray4_, mclVv(d, "d"), NULL),
                NULL));
        mclLoadConditional(
          name_,
          "adaptest",
          mclPrepareGlobal(&adaptest),
          "ans",
          &ans,
          "behavdata",
          &behavdata,
          "c",
          &c,
          "cin",
          &cin,
          "cmap",
          mclPrepareGlobal(&cmap),
          "d",
          &d,
          "data",
          &data,
          "din",
          &din,
          "e",
          &e,
          "ein",
          &ein,
          "model",
          &model,
          "opt",
          &opt,
          "rat",
          &rat,
          "ratio",
          &ratio,
          NULL);
        mclLoadConditional(
          name_, "spikedata", &spikedata, "t", &t, "tin", &tin, NULL);
        mxDestroyArray(name_);
    }
    /*
     * model= adaptest{d}{e}{t}{c}.model;
     */
    mlfAssign(
      &model,
      mlfIndexRef(
        mclVg(&adaptest, "adaptest"),
        "{?}{?}{?}{?}.model",
        mclVv(d, "d"),
        mclVv(e, "e"),
        mclVv(t, "t"),
        mclVv(c, "c")));
    /*
     * 
     * global adaptest cmap
     * %load /bach/AdaptFilter/newvis/colormap;
     * load /bach/AdaptFilter/newvis/colormap2;
     */
    {
        mxArray * name_ = mclInitialize(_mxarray10_);
        mclLoadConditional(
          name_,
          "adaptest",
          mclPrepareGlobal(&adaptest),
          "ans",
          &ans,
          "behavdata",
          &behavdata,
          "c",
          &c,
          "cin",
          &cin,
          "cmap",
          mclPrepareGlobal(&cmap),
          "d",
          &d,
          "data",
          &data,
          "din",
          &din,
          "e",
          &e,
          "ein",
          &ein,
          "model",
          &model,
          "opt",
          &opt,
          "rat",
          &rat,
          "ratio",
          &ratio,
          NULL);
        mclLoadConditional(
          name_, "spikedata", &spikedata, "t", &t, "tin", &tin, NULL);
        mxDestroyArray(name_);
    }
    /*
     * opt.playrate= str2num(ratio)*500;
     */
    mlfIndexAssign(
      &opt,
      ".playrate",
      mclMtimes(mlfStr2num(NULL, mclVa(ratio, "ratio")), _mxarray12_));
    /*
     * opt.cmap= cmap;
     */
    mlfIndexAssign(&opt, ".cmap", mclVg(&cmap, "cmap"));
    /*
     * %opt.cmap= jet(1024);
     * 
     * visModel(data, model, opt)
     */
    mclPrintAns(
      &ans,
      mlfNVisModel(
        0,
        mclAnsVarargout(),
        mclVv(data, "data"),
        mclVv(model, "model"),
        mclVv(opt, "opt"),
        NULL));
    mxDestroyArray(d);
    mxDestroyArray(e);
    mxDestroyArray(t);
    mxDestroyArray(c);
    mxDestroyArray(ans);
    mxDestroyArray(behavdata);
    mxDestroyArray(data);
    mxDestroyArray(spikedata);
    mxDestroyArray(model);
    mxDestroyArray(opt);
    mxDestroyArray(ratio);
    mxDestroyArray(cin);
    mxDestroyArray(tin);
    mxDestroyArray(ein);
    mxDestroyArray(din);
    mxDestroyArray(rat);
    mclSetCurrentLocalFunctionTable(save_local_function_table_);
}
