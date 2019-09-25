/* $Id: VStatAM.cc,v 1.2 2008/08/20 21:47:14 chengs Exp $
   
   Sen Cheng, Fri Sep  3 12:58:42 PDT 2004
   
   program description
*/

#include "VStatAM.h"

/* ************************************************************
                          class VStatAM
   ************************************************************ */

using namespace AFilter;

VStatAM::VStatAM(const TData *d, const AdaptModel *m, const StatLimits *l)
    : VStat(d,m,l), model(m)
{
    //do something;
}

VStatAM::VStatAM(VStatAM &)
{
    //error;
}

VStatAM::~VStatAM()
{
    //empty;
}

void VStatAM::init_obj()
{
}

void VStatAM::run(const TData *data, const AdaptModel *model, 
		const StatLimits *lim)
{
}

