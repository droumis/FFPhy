// $Id: colormapTexture.cc,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#include "colormapTexture.h"

ColormapTexture::ColormapTexture()
{
	width = 8;
	height = 256;

	pixels = new color3f[width * height];	
	
	memset(pixels, 0, width*height*sizeof(color3f)); // zero image memory

	texId = 0; // this will appear as default texture until initialized
	initialized = false;
}

ColormapTexture::~ColormapTexture()
{
	delete[] pixels;
}

void ColormapTexture::setColormap(int cLen, color3f *cMap)
{
	this->colormapLength = cLen;
	this->colormap = cMap;

	for (int r=0; r<height; r++)
		for (int c=0; c<width; c++)
		{
			int pi = r*width+c;
			int ci = (int) (r/(height/(double)colormapLength));
			
			pixels[pi].r = colormap[ci].r;
			pixels[pi].g = colormap[ci].g;
			pixels[pi].b = colormap[ci].b;
		}

	if (initialized) subImage();
}

void ColormapTexture::init()
{
	if (initialized) return;
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glGenTextures(1, &texId);
	bind();
	
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, 
			GL_RGB, GL_FLOAT, pixels);

	// int err = glGetError();
	// assert(err!=0);
	// if (err!=0)
		// cerr << "OpenGL error " << err << " occurred.\n";

	initialized = true;
}

void ColormapTexture::subImage() const
{
	assert(initialized);
	bind();
	glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, 
			GL_RGB, GL_FLOAT, pixels);
}
