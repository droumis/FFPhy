// $Id: meshWindow.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#ifndef __MESHWINDOW_H__
#define __MESHWINDOW_H__

#include "meshWindow.h"
#include "eventListeners.h"
#include "camera.h"

class MeshWindow : public ViewWindow, public TrajMeshEventListener
{
	public:
		MeshWindow::MeshWindow(int x, int y, int w, int h, const char *l);
		~MeshWindow() 
		{ if (dataModel!=NULL) dataModel->removeTrajMeshEventListener(traj, this); }

		void setDataModelAndTrajectory(AFDataModel *d, int t) 
		{ 
			this->traj=t; 
			ViewWindow::setDataModel(d); 
			dataModel->addTrajMeshEventListener(traj, this); 
		}

		/* OpenGL related fns. */
		void glInit();
		void glDraw();
	
		/* event handling */
		int handle(int event);
		
		void trajMeshEventHandle() { this->redraw(); }
			
	private:
		int traj;

		Camera camera;
};

#endif


