// $Id: viewWindow.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#ifndef __VIEWWINDOW_H__
#define __VIEWWINDOW_H__

#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>

#include "afDataModel.h"

#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

class ViewWindow : public Fl_Gl_Window {
	public:
		ViewWindow(int x, int y, int w, int h, const char *l);
		virtual ~ViewWindow() {}

		// TODO
		// should this be in ctor?
		virtual void setDataModel(AFDataModel *d) 
		{ this->dataModel = d; }
	
		AFDataModel* getDataModel() { return dataModel; }
		
		virtual void glInit()=0;
		virtual void glReshape() { glViewport(0,0,w(),h()); }
		virtual void glDraw()=0;
		
		/* these should typically not be overriden */
		virtual void init() { glInit(); glReshape(); initialized = true; }
		virtual void draw();
		
	protected:
		AFDataModel *dataModel;
		
		bool initialized;
};

#endif
