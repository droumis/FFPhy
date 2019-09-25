/*
 * MATLAB Compiler: 3.0
 * Date: Wed Dec  6 10:10:36 2006
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-m" "-W" "main" "-L"
 * "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "debugVis" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#include "libmatlb.h"
#include "debugVis.h"
#include "visModel_mex_interface.h"
#include "libmmfile.h"

extern _mex_information _main_info;

mxArray * adaptest = NULL;

mxArray * cmap = NULL;

static mexGlobalTableEntry global_table[2]
  = { { "adaptest", &adaptest }, { "cmap", &cmap } };

static mexFunctionTableEntry function_table[2]
  = { { "debugVis", mlxDebugVis, 6, 0, &_local_function_table_debugVis },
      { "visModel", mlxVisModel, -1, -1, &_local_function_table_visModel } };

static const char * path_list_[1] = { "/home/chengs/AdaptFilter/main" };

static _mexInitTermTableEntry init_term_table[3]
  = { { libmmfileInitialize, libmmfileTerminate },
      { InitializeModule_debugVis, TerminateModule_debugVis },
      { InitializeModule_visModel_mex_interface,
        TerminateModule_visModel_mex_interface } };

_mex_information _main_info
  = { 1, 2, function_table, 2, global_table, 1,
      path_list_, 3, init_term_table };

/*
 * The function "main" is a Compiler-generated main wrapper, suitable for
 * building a stand-alone application.  It calls a library function to perform
 * initialization, call the main function, and perform library termination.
 */
int main(int argc, const char * * argv) {
    return mclMain(argc, argv, mlxDebugVis, 0, &_main_info);
}
