// $Id: texture.h,v 1.2 2008/08/24 19:47:13 chengs Exp $
#ifndef __TEXTURE_H__
#define __TEXTURE_H__

#include <FL/gl.h>
//#include <GL/glu.h>
#include "glUtils.h"

#include <assert.h>
#include <iostream>

class Texture 
{
	public:
		Texture(int w=256, int h=256) : width(w), height(h)
		{ 	
			texData = new color3f[width*height]; 	
			memset(texData,0,width*height*sizeof(color3f));
		}

		virtual ~Texture() { delete[] texData; }

		void reset()
		{		
			memset(texData,0,width*height*sizeof(color3f));
			subImage();
		}

		void bind() { glBindTexture(GL_TEXTURE_2D,texId); }

		void subImage() 
		{ 
			bind(); 	
			glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height,
					GL_RGB, GL_FLOAT, texData); 
		} 
		
		void init()
		{	
			glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
			
			glGenTextures(1, &texId);
//            std::cerr << "Texture assigned texId " << texId << "\n";
			
			bind();
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_FLOAT, texData);
		}

		// no-op for now. Will not be a no-op if using 
		// textures stored in video ram
		void dataBegin() { }

		// May only be called between a call to
		// dataBegin() and dataEnd()
		// x and y are column, row respectively from lower left
		void setPixel(int x, int y, float r, float g, float b)
		{
			int i = pToI(x,y);
			
			// std::cerr << "setting pixel x=" << x << "   y=" << y << "   i=" << i << "...";
			assert(i<width*height);
			
			texData[i].r=r;
			texData[i].g=g;
			texData[i].b=b;
		}
		
		// here x and y are normalized to range [0,1], lower left is (0,0) 
		void setPixelF(float x, float y, float r, float g, float b)
		{ setPixel((int)(x*(width-1)), (int)(y*(height-1)), r, g, b); }

		void dataEnd() { subImage(); }

	private:
		// convert pixel coords (x,y from lower left) to index into texData
		int pToI(int x, int y)
		{ return y*width+x; }

	private:

		int width, height;
		GLuint texId;
		color3f *texData;		// holds our texture data
};

#endif
