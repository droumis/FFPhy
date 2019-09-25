/*
 * MATLAB Compiler: 3.0
 * Date: Wed Dec  6 10:10:35 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "debugVis" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __debugVis_h
#define __debugVis_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_debugVis(void);
extern void TerminateModule_debugVis(void);
extern _mexLocalFunctionTable _local_function_table_debugVis;

extern void mlfDebugVis(mxArray * rat,
                        mxArray * din,
                        mxArray * ein,
                        mxArray * tin,
                        mxArray * cin,
                        mxArray * ratio);
extern void mlxDebugVis(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
