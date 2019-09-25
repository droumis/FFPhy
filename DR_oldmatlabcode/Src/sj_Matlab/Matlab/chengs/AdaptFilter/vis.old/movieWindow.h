// $Id: movieWindow.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#ifndef __MOVIEWINDOW_H__
#define __MOVIEWINDOW_H__

#include "viewWindow.h"
#include "eventListeners.h"
#include "utils.h"

using namespace std;

class MovieWindow : public ViewWindow, public MovieEventListener
{
	public:
		MovieWindow::MovieWindow(int x, int y, int w, int h, const char *l);
		~MovieWindow() 
		{ if (dataModel!=NULL) dataModel->removeMovieEventListener(this); }

		void setDataModel(AFDataModel *d)
		{ ViewWindow::setDataModel(d); 
			dataModel->addMovieEventListener(this); }
		
		/* OpenGL related fns. */
		void glInit();
		void glDraw();

		/*event handling */
		void movieEventHandle() { this->redraw(); }
};

#endif
