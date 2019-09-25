// $Id: axes.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#ifndef __AXES_H__
#define __AXES_H__

#include <FL/gl.h>

// these require that GL_COLOR_MATERIAL be enabled
// and that a font be set
void drawAxis(float xMin, float xMax,
		float xIncMinor, float xIncMajor,
		float tickSizeMinor, float tickSizeMajor,
		float numOffset);

// the x's are given in log units
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
