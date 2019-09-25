// $Id: glUtils.cc,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#include "glUtils.h"
#include "texture.h"

std::ostream& operator<<(std::ostream& o, const vertex_tmesh& v)
{
	return o << "(" << v.r << ", " << v.g << ", " << v.b 
		<< ", " << v.x << ", " << v.y << ", " << v.z << ")";
}

void drawTexturedQuad(Texture &tex, 
		float l, float b, 
		float r, float t)
{
	tex.bind();
	glBegin(GL_QUADS);
	glTexCoord2f(0.0f,0.0f);
	glVertex3f(l,b,0);
	glTexCoord2f(0.0f,1.0f);
	glVertex3f(l,t,0);
	glTexCoord2f(1.0f,1.0f); 
	glVertex3f(r,t,0);
	glTexCoord2f(1.0f,0.0f); 
	glVertex3f(r,b,0);
	glEnd();
}
