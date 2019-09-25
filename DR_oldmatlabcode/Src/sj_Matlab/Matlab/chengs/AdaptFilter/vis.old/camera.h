// $Id: camera.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#ifndef __CAMERA_H__
#define __CAMERA_H__

#include <FL/gl.h>

class Camera {
	public:
		Camera();

		// mouse events
		void mouseDown(int button, int x, int y);
		void mouseDrag(int eventState, int x, int y);
		void mouseUp(int button);
		void mouseWheel(int dx, int dy);

		/* Modifies the current OpenGL matrix. */
		void transform() 
		{
			glTranslatef(xShift, yShift, 0);
			glRotatef(dec,1,0,0); glRotatef(rot,0,1,0);
			glScalef(float(zm),float(zm),float(zm));
		}
		
	private:
		// last mouse position
		int mouseX, mouseY;

		// camera positioning
		// float hAng, vAng;
		double rot, dec;
		double xShift, yShift;
		float zm;
};

#endif
