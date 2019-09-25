// $Id: axes.h,v 1.2 2008/08/24 19:47:11 chengs Exp $
#ifndef __AXES_H__
#define __AXES_H__

#include <FL/gl.h>

// these require that GL_COLOR_MATERIAL be enabled
// and that a font be set
// the x's are given in log units
void drawAxis(float xMin, float xMax,
		float xIncMinor, float xIncMajor,
		float tickSizeMinor, float tickSizeMajor,
		float numOffset);

void drawAxisNew(float xMin, float xMax,
		float xIncMinor, float xIncMajor,
		float tickSizeMinor, float tickSizeMajor,
		float numOffset, float tickOffset,
        bool logflag= false);

void drawLogAxis(float xMin, float xMax,
		float xIncMinor, float xIncMajor,
		float tickSizeMinor, float tickSizeMajor,
		float numOffset);

void drawGrid(float xMin, float xMax,
		float xIncMinor, float xIncMajor,
		float xNumOffset,
		float yMin, float yMax,
		float yIncMinor, float yIncMajor,
		float yNumOffset);

#endif
