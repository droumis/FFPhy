// $Id: viewWindow.cc,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#include "viewWindow.h"

ViewWindow::ViewWindow(int x, int y, int w, int h, const char *l)
: Fl_Gl_Window(x,y,w,h,l)
{
	initialized = false;
}

void ViewWindow::draw() {
	if (!initialized)	
		init();
	else if (!valid())
		glReshape();
	glDraw();
}


