#ifndef VSTAT2D_H
#define VSTAT2D_H

/* $Id: VStat2d.h,v 1.3 2008/08/20 21:47:14 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "VStat.h"
#include "../model/VDynamics2d.h"

/* ************************************************************
                          class VStat2d
   ************************************************************ */

namespace AFilter {

class VStat2d : public VStat
{
private:
protected:
    const VDynamics2d *dyn2;
protected:
    VStat2d(const TData*, const AdaptModel*, const StatLimits*);
    VStat2d(VStat2d &);
public:
    virtual ~VStat2d();

    virtual void init_obj();
};

} // namespace AFilter

#endif   // VSTAT2D_H
