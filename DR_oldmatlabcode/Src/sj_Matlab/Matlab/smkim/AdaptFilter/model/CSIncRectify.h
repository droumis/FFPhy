#ifndef CSIncRectify_H
#define CSIncRectify_H

/* $Id: CSIncRectify.h,v 1.2 2008/02/01 19:46:49 chengs Exp $

   Sen Cheng, Tue Jun  6 14:06:42 PDT 2006
   definition of class CSIncRectify
*/

#include "CSIncrement.h"

/* ************************************************************
                          class CSIncRectify
   ************************************************************ */
namespace AFilter {
    
class CSIncRectify : public CSIncrement {
protected: // methods
public:
    CSIncRectify();
    virtual ~CSIncRectify();
    virtual CSIncRectify* alloc_new() const { return new CSIncRectify; }

    virtual const char* getName() const {return "CSIncRectify"; };
//    virtual void setAllParam(int t, double z[]) const;
    virtual void moveTo(int t) const;


}; // class CSIncRectify

}; //namespace AFilter

#endif   // CSIncRectify_H
