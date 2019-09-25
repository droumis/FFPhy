// $Id: afGeomPlot.cc,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#include "afGeomPlot.h"
#include "../aux/numerics.h"

AFGeomPlot::AFGeomPlot(TrajInfo *ti)
	: trajInfo(ti), posTex(256,256)
{
	initialized = false;
}

AFGeomPlot::~AFGeomPlot()
{
}

void AFGeomPlot::init()
{
	if (!initialized) {
		posTex.init();
		computePosTex();
	}
	initialized = true;
}

void AFGeomPlot::draw()
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
	drawTexturedQuad(posTex, 0.0f, 0.0f, 1.0f, 1.0f);
	glDisable(GL_TEXTURE_2D);
	
	// draw rat's position
	glPointSize(3.3f);
	float _x = normalize(ratPos.x, trajInfo->getMinXPos(),
			trajInfo->getMaxXPos());
	float _y = normalize(ratPos.y, trajInfo->getMinYPos(),
			trajInfo->getMaxYPos());	
	glBegin(GL_POINTS);
	glColor3f(1.0f,1.0f,1.0f);
	glVertex3f(_x,_y,0.0f);
	glEnd();
}

void AFGeomPlot::compute(int timestep)
{
	ratPos.x = trajInfo->getData()->getXPos(timestep);
	ratPos.y = trajInfo->getData()->getYPos(timestep);
}

void AFGeomPlot::computePosTex()
{

	cerr << "Computing track geometry texture...";
	
	const TData *data = trajInfo->getData();

	double xPosMin = trajInfo->getMinXPos();
	double xPosMax = trajInfo->getMaxXPos();
	double yPosMin = trajInfo->getMinYPos();
	double yPosMax = trajInfo->getMaxYPos();
	
	posTex.dataBegin();
	for (int t=0; t<data->getNTimesteps(); t++) {

		// TODO this computation may be incorrect
		
		float _x = normalize(data->getXPos(t),xPosMin,xPosMax);
		float _y = normalize(data->getYPos(t),yPosMin,yPosMax);
		
		posTex.setPixelF(_x,_y,0.3f,0.3f,0.7f);
	}
	
	posTex.dataEnd();

	cerr << " complete.\n";
}
