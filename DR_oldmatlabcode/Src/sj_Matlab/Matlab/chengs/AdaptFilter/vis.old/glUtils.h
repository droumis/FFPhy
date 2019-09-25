// $Id: glUtils.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#ifndef __GLUTILS_H__
#define __GLUTILS_H__

#include <Fl/gl.h>

#include <iostream>

/* pack these struct so that we can use it as OpenGL vertices
*/

#if !defined (__GCC__)
#pragma pack()
#endif
struct vertex2f {
	vertex2f()
		:x(0),y(0)
		{}
	
	vertex2f(float x, float y)
		: x(x), y(y)
		{}

	vertex2f(const vertex2f &v)
		:x(v.x), y(v.y)
		{}
	
	vertex2f& operator=(vertex2f &v)
	{ this->x=v.x; this->y=v.y; return *this; }
	
	float x;
	float y;
#if defined (__GCC__)
} __attribute__ ((packed));
#else
};
#endif


#if !defined (__GCC__)
#pragma pack()
#endif
struct vertex3f {
	float x;
	float y;
	float z;
#if defined (__GCC__)
} __attribute__ ((packed));
#else
};
#endif

#if !defined (__GCC__)
#pragma pack()
#endif
struct color3f {
	float r;
	float g;
	float b;
#if defined (__GCC__)
} __attribute__ ((packed));
#else
};
#endif


#if !defined(__GCC__)
#pragma pack()
#endif
struct vertex_tmesh {
	float r,g,b;
	float x,y,z;

	friend std::ostream& operator<<(std::ostream& o, const vertex_tmesh& v);
#if defined (__GCC__)
} __attribute__ ((packed));
#else
};
#endif

/* If you are drawing several textured quads sequentially,
*  don't use this fn. This fn. calls glBegin(GL_QUADS), glEnd().
*  If you are drawing n quads, you don't want to repeat those calls
*  n times.
*/
/* Make sure that 2d textures are enabled when calling this! */
class Texture;
void drawTexturedQuad(Texture &tex, 
		float left, float bottom, 
		float right, float top);

#endif
