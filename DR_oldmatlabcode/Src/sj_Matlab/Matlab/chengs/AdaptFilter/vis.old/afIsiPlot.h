// $Id: afIsiPlot.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#ifndef __AFISITEXTURE_H__
#define __AFISITEXTURE_H__

// OpenGL includes
#include <FL/gl.h>
#include <GL/glu.h>
#include <assert.h>

// vertex_tmesh, vertex3f, ManageMatrix etc. are defined in here
#include "glUtils.h"
#include "trajInfo.h"

using namespace std;

class AFIsiPlot {
	public:
		AFIsiPlot(TrajInfo *trajInfo);
		~AFIsiPlot();

		/*
		* The following require a current OpenGL context.
		*/
		void init();
		void draw();
		void compute(int time);
		
	private:
		void createVertexArray();
		
	private:
		bool initialized;
		TrajInfo *trajInfo;

		int width;
		int height;

		int nVerts;
		vertex3f *vertices;

        float xMin, xMax, xRes;
        float vMax;
};

#endif
