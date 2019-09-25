// $Id: Display.cc,v 1.1 2008/08/24 19:47:08 chengs Exp $
#include "Display.h"

#include "../aux/hippoIO.h"

Display::Display(int x, int y, int w, int h, const char *l)
: Fl_Gl_Window(x,y,w,h,l)
{
	initialized = false;
//    hippo_Mark;
}

Display::~Display() 
{
//    hippo_Mark; 
}

void Display::draw() {
	if (!initialized)	
		init();
	else if (!valid())
		glReshape();
	glDraw();
}


