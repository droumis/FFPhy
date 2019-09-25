// $Id: globals.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#include "afViewUI.h"

extern AFViewUI *avui;


void movie_idle(void *v) {
  avui->viewWindow->getDataModel()->doFrame();
}

#endif
