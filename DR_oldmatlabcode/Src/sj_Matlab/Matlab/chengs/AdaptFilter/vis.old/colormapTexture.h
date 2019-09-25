// $Id: colormapTexture.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#ifndef __COLORMAPTEXTURE_H__ 
#define __COLORMAPTEXTURE_H__

// OpenGL includes
#include <FL/gl.h>
#include <GL/glu.h>
#include <assert.h>

// vertex_tmesh, vertex3f, ManageMatrix etc. are defined in here
#include "glUtils.h"

using namespace std;

class ColormapTexture {
	public:
		ColormapTexture(); 
		~ColormapTexture();
	
		// we assume that the colormap has 2^n, n<=8, elements 
		void setColormap(int cLen, color3f *c);
		
		/* 
		*  The following require a current OpenGL context 
		*  */
		/* Get a texture ID, set the texture interpolation modes, etc. */
		void init();
		
		/* Set this as the texture for OpenGL to render. */
		void bind() const { glBindTexture(GL_TEXTURE_2D,texId); }

		/* Call bind and then update the texture */	
		void subImage() const;
			
	private:
		color3f *colormap;
		int colormapLength;
		
		/* these should be powers of 2 */
		int width;
		int height;
	
		bool initialized;
		
		GLuint texId; // the openGL assigned texture identifier
		color3f *pixels; // the (row-major) image data
};

 #endif
