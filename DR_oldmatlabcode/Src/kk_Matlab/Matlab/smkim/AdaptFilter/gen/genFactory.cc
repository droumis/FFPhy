/* $Id: genFactory.cc,v 1.4 2008/08/24 19:47:02 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "genFactory.h"
#include "../aux/hippoIO.h"
#include "../aux/mexAux.h"

#include "../gen/RescalingGen.h"
#include "../gen/PoissonGen.h"
 

using namespace AFilter;

VGenerator* AFilter::allocGenerator(const mxArray *opts, TData *data, AdaptModel* model)
{
    VGenerator *gen;
    char *name;
    mxArray *mxTmp= mxGetField(opts,0,"name");
    hippo_Assert(mxTmp, "gen.name was not defined");
    MX_StringAllocAssign(mxTmp,name);     

    if (strcmp(name,"RescalingGen")==0) {
        gen= new RescalingGen;
    } else if (strcmp(name,"PoissonGen")==0) {
        gen= new PoissonGen;
    } else {
        gen= 0;
        MX_Error("unknown generator: %s ", name);
    }
    gen->setMexInput(opts);

    mxFree(name);
    return gen;
}
