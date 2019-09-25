// $Id: Geometry.cc,v 1.1 2008/08/24 19:47:09 chengs Exp $
#include "Geometry.h"
#include "../aux/numerics.h"

using namespace AFilter;

Geometry::Geometry(DataCache *c)
	: cache(c), posTex(256,256)
{
	initialized = false;
}

Geometry::~Geometry()
{

}

void Geometry::init()
{
	if (!initialized) {
		posTex.init();
        computePosTex();
	}
	initialized = true;
}

void Geometry::draw()
{

	// draw outline
	// glEnable(GL_COLOR_MATERIAL);
	// glColor3f(1.0f,1.0f,1.0f);
// 
	// glBegin(GL_LINE_LOOP);
	// glVertex3f(0.0f,0.0f,0.0f);
	// glVertex3f(0.0f,1.0f,0.0f);
	// glVertex3f(1.0f,1.0f,0.0f);
	// glVertex3f(1.0f,0.0f,0.0f);
	// glEnd();
// 
	// glDisable(GL_COLOR_MATERIAL);
	
	// draw posTex
    glEnable(GL_TEXTURE_2D);
	glColor3f(0.3f,0.3f,0.7f);
    drawTexturedQuad(posTex, 0.0f, 0.0f, 1.0f, 1.0f);
    glDisable(GL_TEXTURE_2D);
	
	// draw rat's position
	glPointSize(5.0f);
	float _x = normalize(ratPos.x, cache->getMinXPos(),
			cache->getMaxXPos());
	float _y = normalize(ratPos.y, cache->getMinYPos(),
			cache->getMaxYPos());	
	glBegin(GL_POINTS);
	glColor3f(1.0f,1.0f,1.0f);
	glVertex3f(_x,_y,0.0f);
	glEnd();
}

void Geometry::compute(int timestep)
{
	ratPos.x = cache->getData()->getXPos(timestep);
	ratPos.y = cache->getData()->getYPos(timestep);
//        printf("%d: %.2f, %.2f\n", timestep, ratPos.x, ratPos.y);
}

void Geometry::computePosTex()
{

//    cerr << "Computing track geometry texture...";
	
	const TData *data = cache->getData();

	double xPosMin = cache->getMinXPos();
	double xPosMax = cache->getMaxXPos();
	double yPosMin = cache->getMinYPos();
	double yPosMax = cache->getMaxYPos();
	
	posTex.dataBegin();
	for (int t=0; t<data->getNTimesteps(); t++) {

		// TODO this computation may be incorrect
		
		float _x = normalize(data->getXPos(t),xPosMin,xPosMax);
		float _y = normalize(data->getYPos(t),yPosMin,yPosMax);
		
		posTex.setPixelF(_x,_y,0.3f,0.3f,0.7f);
//        if(t % 100==0) printf("%d: %.2f, %.2f\n", t, _x, _y);
	}
	
	posTex.dataEnd();

//    cerr << " complete.\n";
}
