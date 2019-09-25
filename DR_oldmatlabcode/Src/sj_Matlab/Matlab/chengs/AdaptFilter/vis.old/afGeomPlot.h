// $Id: afGeomPlot.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#ifndef __AFGEOMPLOT_H__
#define __AFGEOMPLOT_H__

#include <FL/gl.h>
#include <GL/glu.h>
#include <assert.h>

#include "glUtils.h"
#include "trajInfo.h"

#include "texture.h"

using namespace std;

class AFGeomPlot {

	public:
		AFGeomPlot(TrajInfo *trajInfo);
		~AFGeomPlot();

		void init();
		void draw();
		void compute(int timestep);
		
	private:
		bool initialized;
		
		TrajInfo *trajInfo;

		// every point that the rat visited is drawn on this texture
		Texture posTex;

		void computePosTex(); // draw the locus of rat pts on the texture
		
		vertex2f ratPos;
};

#endif
