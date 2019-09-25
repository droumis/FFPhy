#ifndef STATSFACTORY_H
#define STATSFACTORY_H

/* $Id: statsFactory.h,v 1.1.1.1 2006/04/24 14:21:47 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/


#include "../aux/defs.h"
#include "../aux/mexAux.h"
#include "StatLimits.h"
#include "VStatAM.h"
#include "VStat1d.h"
#include "VStat2d.h"

/* ************************************************************
   ************************************************************ */

namespace AFilter {
	class VStatAM;
	class VStat1d;
	class VStat2d;

	void allocStatAMList(const mxArray *opts, VStatAM **&amstat, int &n);
	void allocStat1dList(const mxArray *opts, VStat1d **&vstat, int &n);
	void allocStat2dList(const mxArray *opts, VStat2d **&vstat, int &n);

	VStatAM* allocStatAM(const mxArray *opts, bool warn= true);
	VStat1d* allocStat1d(const mxArray *opts, bool warn= true);
	VStat2d* allocStat2d(const mxArray *opts, bool warn= true);

    VStat* allocStat(const mxArray *opts, TData * data, AdaptModel *model,
        StatLimits *lim, bool warn= true);
} // namespace AFilter

#endif   // STATSFACTORY_H
