// $Id: camera.cc,v 1.2 2008/08/24 19:47:11 chengs Exp $
#include "camera.h"

Camera::Camera()
{
	mouseX = 0;
	mouseY = 0;
	
	rot=0.0;
	dec=0.0;
	zm=10.0;
	xShift=0.0;
	yShift=0.0;
}

void Camera::mouseDown(int button, int x, int y)
{
	mouseX = x;
	mouseY = y;
}

void Camera::mouseDrag(int eventState, int x, int y)
{
	int dx = x - mouseX;
	int dy = y - mouseY;

	float mPan = 0.10f;
	float mAngle = 2.0f;

	if (eventState & FL_BUTTON1) {
		xShift += mPan*dx;
		yShift += -mPan*dy;
	}
	else if (eventState & FL_BUTTON3) {
		rot += mAngle*dx;
		dec += mAngle*dy;
	}    

	mouseX = x;
	mouseY = y;
}

void Camera::mouseUp(int button)
{}

void Camera::mouseWheel(int dx, int dy)
{
	zm += dy;
}
