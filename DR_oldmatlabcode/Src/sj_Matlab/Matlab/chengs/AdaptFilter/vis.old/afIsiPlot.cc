// $Id: afIsiPlot.cc,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#include "afIsiPlot.h"
#include "axes.h"

AFIsiPlot::AFIsiPlot(TrajInfo *ti)
: trajInfo(ti)
{
	nVerts = 0;
	
	initialized = false;	
}

AFIsiPlot::~AFIsiPlot()
{}

void AFIsiPlot::init()
{
	xMin = (float) log(0.001);
	xMax = (float) log(trajInfo->getMaxTIsi());
    vMax = trajInfo->approxV1Max() *1.2;
//    vMax = 10.0f;
    nVerts = 200;

    if (!initialized) {
        createVertexArray();
    }
	xRes = (xMax-xMin)/(nVerts-1);

    initialized = true;
}

void AFIsiPlot::compute(int time)
{
	assert(initialized);

	float x = xMin;
	for (int i=0; i<nVerts; i++) {
		vertices[i].y = trajInfo->evalIsi(time, exp(x), 0)/vMax;	
		x+=xRes;
	}
}

//void AFIsiPlot::compute(int time)
//{
//    assert(initialized);

//    float xMin = trajInfo->getMinTIsi();
//    float xMax = trajInfo->getMaxTIsi();
//    float xRes = (xMax-xMin)/(nVerts-1);
//    
//    const float VE = 1.1f;
//        
//    float vMax = trajInfo->approxV1Max()*VE;

//    float x = xMin;
//    for (int i=0; i<nVerts; i++) {
//        vertices[i].y = trajInfo->evalIsi(time, x, 0)/vMax;	
//        x+=xRes;
//    }
//}

void AFIsiPlot::draw()
{
    assert(initialized);

    glDisable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    glColor3f(0.2,0.2,0.7);

    glLineWidth(1.0f);

    glEnableClientState(GL_VERTEX_ARRAY);

    glVertexPointer(3,GL_FLOAT,sizeof(vertex3f),vertices);

    glDrawArrays(GL_LINE_STRIP,0,nVerts);

    glPushMatrix();
    glTranslatef(0.0f,1.0f,0.0f);
    glScalef(1.0f,-1.0f,1.0f);
    drawLogAxis(xMin, xMax, (xMax-xMin)/10, (xMax-xMin)/5,
            0.02f,0.07f,0.02f);

    glPopMatrix();
}

void AFIsiPlot::createVertexArray()
{
    vertices = new vertex3f[nVerts];

    /* curve vertices have x,y coordinates in range [0,1]
     *  z coordinate is 0
     */
    float xExtent = 1.0f;
    float xRes = xExtent/nVerts;

    for (int i=0; i<nVerts; i++) {
        vertices[i].x = xRes*i;
        vertices[i].y = (i%2==0) ? 0.0f : 1.0f;
        vertices[i].z = 0.0f;
    }
}
