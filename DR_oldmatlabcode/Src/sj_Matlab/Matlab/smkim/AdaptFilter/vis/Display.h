// $Id: Display.h,v 1.1 2008/08/24 19:47:08 chengs Exp $
#ifndef __Display_H__
#define __Display_H__

#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>

#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

class Display : public Fl_Gl_Window {
	public:
		Display(int x, int y, int w, int h, const char *l);
		virtual ~Display();

//        virtual void glInit()=0;
//        virtual void glDraw()=0;
		virtual void glInit(){};
		virtual void glDraw(){};
		virtual void glReshape() { glViewport(0,0,w(),h()); }
		
		/* these should typically not be overriden */
		virtual void init() { glInit(); glReshape(); initialized = true; }
		virtual void draw();
		
	protected:
		
		bool initialized;
};

#endif
