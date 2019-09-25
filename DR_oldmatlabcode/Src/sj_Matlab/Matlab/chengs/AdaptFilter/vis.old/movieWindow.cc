// $Id: movieWindow.cc,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#include "movieWindow.h"

MovieWindow::MovieWindow(int x, int y, int w, int h, const char *l)
	: ViewWindow(x,y,w,h,l)
{}

void MovieWindow::glInit()
{
	glClearColor (0.0f, 0.0f, 0.0f, 0.0f);	
	glClearDepth (1.0f);

	glDisable(GL_LIGHTING);

	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	gl_font(FL_COURIER,10);
	
	dataModel->glInit();
}

void MovieWindow::glDraw()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0,1,0,6,-100,100);
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	int nTraj = dataModel->getData()->getNTraj();

	/* draw posphase movies */
	for (int i=0; i<nTraj; i++) 
	{
		float l = 0.0f; 
		float b = 6-i-1;

		glPushMatrix();
		glTranslatef(l, b, 0.0f);
	
		dataModel->getTexture(i)->draw();

		glPopMatrix();
	}

	/* draw geometry plot */
	{
		glPushMatrix();
        glTranslatef(0.5f,1.0f,0.0f);
        glScalef(0.5f,-1.0f,1.0f);

		dataModel->getGeomPlot()->draw();

		glPopMatrix();
	}

	/* draw isi plot */
	glPushMatrix();
	glTranslatef(0.0f,0.0f,0);
	glScalef(0.5f,1.0f,1.0f);
	dataModel->getIsiPlot()->draw();
	glPopMatrix();
	
}
