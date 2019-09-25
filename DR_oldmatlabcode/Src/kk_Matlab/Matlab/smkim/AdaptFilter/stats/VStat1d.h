#ifndef VSTAT1D_H
#define VSTAT1D_H

/* $Id: VStat1d.h,v 1.3 2008/08/20 21:47:14 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "VStat.h"
#include "../model/VDynamics1d.h"

/* ************************************************************
                          class VStat1d
   ************************************************************ */

namespace AFilter {

class VStat1d : public VStat
{
private:
protected:
    const VDynamics1d *dyn1;
protected:
    VStat1d(const TData*, const AdaptModel*, const StatLimits*);
    VStat1d(VStat1d &);
    VStat1d();
public:
    virtual ~VStat1d();

    virtual void init_obj();
};

} // namespace AFilter

#endif   // VSTAT1D_H
