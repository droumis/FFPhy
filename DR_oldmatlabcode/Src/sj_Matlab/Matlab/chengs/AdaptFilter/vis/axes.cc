// $Id: axes.cc,v 1.3 2008/08/24 19:47:11 chengs Exp $
#include "axes.h"
#include "../aux/numerics.h"

#include <stdlib.h>
#include <stdio.h>

using namespace AFilter;

void drawAxisNew(float xMin, float xMax,
		float xIncMinor, float xIncMajor,
		float tickSizeMinor, float tickSizeMajor,
		float numOffset, float tickOffset,
        bool logflag)
{
//    double xMinMinor = getCeilingMultiple(xMin, xIncMinor);
//    double xMinMajor = getCeilingMultiple(xMin, xIncMajor);

    glLineWidth(1.0f);
	glBegin(GL_LINES);

	/* minor ticks  */

	// draw minor ticks on x-axis
//    double xMinor = xMinMinor;
    double xMinor = xMin;
	while (xMinor < xMax) {
		double _x = normalize(xMinor, xMin, xMax); 
		glVertex3f(_x, tickOffset, 0.0f);
		glVertex3f(_x, tickSizeMinor+tickOffset, 0.0f);

		xMinor += xIncMinor;
	}

	/* major ticks */

	// draw major ticks on x-axis
//    double xMajor = xMinMajor;
	double xMajor = xMin;
	while (xMajor < xMax) {
		double _x = normalize(xMajor, xMin, xMax);
		glVertex3f(_x, tickOffset, 0.0f);
		glVertex3f(_x, tickSizeMajor+tickOffset, 0.0f);

		xMajor += xIncMajor;
	}

	glEnd();

	/* numbering */

	char temp[32];

    int precision, length;
    precision= (int) -floor(log(xMax-xMin)/log(10))+1;
    length= fabs(xMin) < fabs(xMax) ? (int) floor(log10(fabs(xMax))) :
        (int) floor(log10(fabs(xMin)));
    if(length < 0) length= 0;
    if(precision<0) precision= 0;
	// number major grid lines on x-axis
//    xMajor = xMinMajor;
	xMajor = xMin;
	while (xMajor < xMax) {
		double _x = normalize(xMajor, xMin, xMax);
		if(logflag) snprintf(temp,32,"%.1e",pow(10,xMajor));
        else snprintf(temp,32,"%*.*f",length,precision,xMajor);
		glRasterPos3f(_x, numOffset, 0.0f);
//        glBitmap( 0, 0, 0, 0, -20, 0, NULL );
		gl_draw(temp);

		xMajor += xIncMajor;
	}
}

void drawAxis(float xMin, float xMax,
		float xIncMinor, float xIncMajor,
		float tickSizeMinor, float tickSizeMajor,
		float numOffset)
{
	double xMinMinor = getCeilingMultiple(xMin, xIncMinor);
	double xMinMajor = getCeilingMultiple(xMin, xIncMajor);

	glBegin(GL_LINES);

	/* minor ticks  */

	// draw minor ticks on x-axis
	double xMinor = xMinMinor;
	while (xMinor < xMax) {
		double _x = normalize(xMinor, xMin, xMax); 
		glVertex3f(_x,          0.0f, 0.0f);
		glVertex3f(_x, tickSizeMinor, 0.0f);

		xMinor += xIncMinor;
	}

	/* major ticks */

	// draw major ticks on x-axis
	double xMajor = xMinMajor;
	while (xMajor < xMax)
	{
		double _x = normalize(xMajor, xMin, xMax);
		glVertex3f(_x,            0.0f, 0.0f);
		glVertex3f(_x,   tickSizeMajor, 0.0f);

		xMajor += xIncMajor;
	}

	glEnd();

	/* numbering */

	char temp[32];

	// number major grid lines on x-axis
	xMajor = xMinMajor;
	while (xMajor < xMax)
	{
		double _x = normalize(xMajor, xMin, xMax);
		snprintf(temp,32,"%.2f",xMajor);
		glRasterPos3f(_x, -numOffset, 0.0f);
//        glBitmap( 0, 0, 0, 0, -20, 0, NULL );
		gl_draw(temp);

		xMajor += xIncMajor;
	}
}

void drawLogAxis(float xMin, float xMax,
		float xIncMinor, float xIncMajor,
		float tickSizeMinor, float tickSizeMajor,
		float numOffset)
{
	double xMinMinor = getCeilingMultiple(xMin, xIncMinor);
	double xMinMajor = getCeilingMultiple(xMin, xIncMajor);

	glBegin(GL_LINES);

	/* minor ticks  */

	// draw minor ticks on x-axis
	double xMinor = xMinMinor;
	while (xMinor < xMax)
	{
		double _x = normalize(xMinor, xMin, xMax); 
		glVertex3f(_x,          0.0f, 0.0f);
		glVertex3f(_x, tickSizeMinor, 0.0f);

		xMinor += xIncMinor;
	}

	/* major ticks */

	// draw major ticks on x-axis
	double xMajor = xMinMajor;
	while (xMajor < xMax)
	{
		double _x = normalize(xMajor, xMin, xMax);
		glVertex3f(_x,            0.0f, 0.0f);
		glVertex3f(_x,   tickSizeMajor, 0.0f);

		xMajor += xIncMajor;
	}

	glEnd();

	/* numbering */

	char temp[32];

	// number major grid lines on x-axis
	xMajor = xMinMajor;
	while (xMajor < xMax)
	{
		double _x = normalize(xMajor, xMin, xMax);
		snprintf(temp,32,"%.2f",exp(xMajor));
		glRasterPos3f(_x, -numOffset, 0.0f);
		gl_draw(temp);

		xMajor += xIncMajor;
	}
}

void drawGrid(float xMin, float xMax,
		float xIncMinor, float xIncMajor,
		float xNumOffset,
		float yMin, float yMax,
		float yIncMinor, float yIncMajor,
		float yNumOffset)
{


	double xMinMinor = getCeilingMultiple(xMin, xIncMinor);
	double xMinMajor = getCeilingMultiple(xMin, xIncMajor);

	double yMinMinor = getCeilingMultiple(yMin, yIncMinor);
	double yMinMajor = getCeilingMultiple(yMin, yIncMajor);

	glLineWidth(0.5f);
	glEnable(GL_COLOR_MATERIAL);
	glBegin(GL_LINES);

	/* minor grid lines */

	// set color for minor grid lines
	glColor3f(0.0f,1.0f,0.0f);

	// draw minor grid lines on x-axis
	double xMinor = xMinMinor;
	while (xMinor < xMax)
	{
		double _x = normalize(xMinor, xMin, xMax); 
		glVertex3f(_x, 0.0f, 0.0f);
		glVertex3f(_x, 1.0f, 0.0f);

		xMinor += xIncMinor;
	}

	// draw minor grid lines on y-axis
	double yMinor = yMinMinor;
	while (yMinor < yMax)
	{
		double _y = normalize(yMinor, yMin, yMax);
		glVertex3f(0.0f, _y, 0.0f);
		glVertex3f(1.0f, _y, 0.0f); 

		yMinor += yIncMinor;
	}

	/* major grid lines */

	// set color for major grid lines
	glColor3f(0.2f,0.8f,0.2f);

	// draw major grid lines on x-axis
	double xMajor = xMinMajor;
	while (xMajor < xMax)
	{
		double _x = normalize(xMajor, xMin, xMax);
		glVertex3f(_x, 0.0f, 0.0f);
		glVertex3f(_x, 1.0f, 0.0f);

		xMajor += xIncMajor;
	}

	// draw major grid lines on y-axis
	double yMajor = yMinMajor;
	while (yMajor < yMax)
	{
		double _y = normalize(yMajor, yMin, yMax);
		glVertex3f(0.0f, _y, 0.0f);
		glVertex3f(1.0f, _y, 0.0f);

		yMajor += yIncMajor;
	}

	glEnd();

	/* numbering */

	char temp[32];

	// number major grid lines on x-axis
	xMajor = xMinMajor;
	while (xMajor < xMax)
	{
		double _x = normalize(xMajor, xMin, xMax);
		snprintf(temp,32,"%.2f",xMajor);
		glRasterPos3f(_x, 1.0f+xNumOffset, 0.0f);
		gl_draw(temp);

		xMajor += xIncMajor;
	}

	// number major grid lines on y-axis
	yMajor = yMinMajor;
	while (yMajor < yMax)
	{
		double _y = normalize(yMajor, yMin, yMax);
		snprintf(temp,32,"%.2f",yMajor);
		glRasterPos3f(1.0f+yNumOffset, _y, 0.0f);
		gl_draw(temp);

		yMajor += yIncMajor;
	}

	glDisable(GL_COLOR_MATERIAL);


}
