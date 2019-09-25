// $Id: meshWindow.cc,v 1.2 2006/11/07 18:02:14 chengs Exp $
#include "meshWindow.h"
#include "utils.h"

MeshWindow::MeshWindow(int x, int y, int w, int h, const char *l)
	: ViewWindow(x,y,w,h,l)
{
	traj = 0;
}

void MeshWindow::glInit()
{
	glClearColor (0.0f, 0.0f, 0.0f, 0.0f);	
	glClearDepth (1.0f);				

	glDepthFunc (GL_LEQUAL);		/* The Type Of Depth Test To Do */
	glEnable (GL_DEPTH_TEST);		/* Enable Depth Testing */

	glShadeModel (GL_SMOOTH);		/* Enables Smooth Color Shading */
	glDisable (GL_LINE_SMOOTH);		/* Initially Disable Line Smoothing */

	dataModel->glInit();

	glDisable(GL_LIGHTING);

	gl_font(FL_COURIER,10);
	
	// enable backface culling
	//	glEnable(GL_CULL_FACE);
	//	glCullFace(GL_BACK);
}

void MeshWindow::glDraw()
{
	// fnTrace(__PRETTY_FUNCTION__);
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);    

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-5,5,-5,5,-100,100);
	
	glPushMatrix();

	camera.transform();

	glPushMatrix();
	// draw trajectory on which the rat is currently running
	dataModel->getMesh(traj)->draw();
	glPopMatrix();
	
	glPopMatrix();

	// draw colormap texture
	// TODO 
	// doesn't work since I changed to dsiplaying meshes in seperate windows.
	// The problem is that the OpenGL context is created and the texId generated
	// in the main movie window, so that texId is invalid in the window in which 
	// the mesh is displayed.
	// glDisable(GL_COLOR_MATERIAL);
	// glEnable(GL_TEXTURE_2D);
	// glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	// dataModel->getColormapTexture()->bind();

	// float l=-4.8f, r=-4.5f, t=4.8f, b=3.5f;
	// 
	// glBegin(GL_QUADS);
	// glTexCoord2f(0.0f,0.0f);
	// glVertex3f(l,b,0);
	// glTexCoord2f(0.0f,1.0f);
	// glVertex3f(l,t,0);
	// glTexCoord2f(1.0f,1.0f); 
	// glVertex3f(r,t,0);
	// glTexCoord2f(1.0f,0.0f); 
	// glVertex3f(r,b,0);
	// glEnd();
// 
	// glDisable(GL_TEXTURE_2D);

	// glEnable(GL_COLOR_MATERIAL);
	// glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
	// glColor3f(1.0,1.0,1.0);	

	// draw box around colormap texture
	// glLineWidth(2.0f);
// 
	// glBegin(GL_LINE_LOOP);
	// glVertex3f(l,b,0);
	// glVertex3f(l,t,0);
	// glVertex3f(r,t,0);
	// glVertex3f(r,b,0);
	// glEnd();

	// display trajectory numbers
	char temp[256];
	glRasterPos3f(-4.8f,-4.5f,0.0f);
	snprintf(temp,"Trajectory: %d",traj);
	gl_draw(temp);
}

int MeshWindow::handle(int event)
{
	switch (event) {
		case FL_PUSH:
			camera.mouseDown(Fl::event_button(), Fl::event_x(), Fl::event_y());
			return 1;
		case FL_DRAG:
			camera.mouseDrag(Fl::event_state(), Fl::event_x(), Fl::event_y());
			redraw();
			return 1;
		case FL_RELEASE:
			camera.mouseUp(Fl::event_button());
			return 1;
		case FL_MOUSEWHEEL:
			camera.mouseWheel(Fl::event_dx(), Fl::event_dy());
			redraw();
			return 1;
		default:
			return Fl_Gl_Window::handle(event);
	}
}
