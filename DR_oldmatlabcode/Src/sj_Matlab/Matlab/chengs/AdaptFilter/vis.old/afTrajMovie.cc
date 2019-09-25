// $Id: afTrajMovie.cc,v 1.2 2006/11/07 18:02:14 chengs Exp $
#include "afTrajMovie.h"
#include "axes.h"
#include "../aux/numerics.h"

#include <iostream>


const int AFTrajMovie::SPIKE_EXIST = 10000;

// the coordinates of the texture rect's corners
const float AFTrajMovie::TX_L = 0.08f;
const float AFTrajMovie::TX_R = 0.95f;
const float AFTrajMovie::TX_B = 0.12f;
const float AFTrajMovie::TX_T = 0.90f;

// the offset of the axes from the texture
const float AFTrajMovie::AO_L = 0.03f;
const float AFTrajMovie::AO_B = 0.08f;

AFTrajMovie::AFTrajMovie(int n, TrajInfo *ti)
:	traj(n), trajInfo(ti)
{
	cerr << traj << "\t";

	width = 128;
	height = 64;

	/* setup image data and OpenGL texture */
	pixels = new color3f[width * height];	

//	cerr << "stride of pixels is " << &pixels[1] << " - " << &pixels[0] << " = " 
//		<< (char*)&pixels[1] - (char*)&pixels[0] << "\n";
//	cerr << "sizeof ( color3f ) = " << sizeof(color3f) << "\n";

	memset(pixels, 0, width*height*sizeof(color3f)); // zero image memory

	texId = 0; // this will appear as default texture until initialized
	initialized = false;

	active = false;
	ratPos = 0.0f;
}

AFTrajMovie::~AFTrajMovie()
{
	delete[] pixels;
}

void AFTrajMovie::setColormap(int cLen, color3f *c)
{
	this->colormapLength = cLen;
	this->colormap = c;
}

void AFTrajMovie::init()
{
	if (!initialized) {
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glGenTextures(1, &texId);
	cerr << "Texture for traj " << traj << " assigned texId " << texId << "\n";
	bind();
	
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, 
			GL_RGB, GL_FLOAT, pixels);

	// int err = glGetError();
	// assert(err!=0);
	// if (err!=0)
		// cerr << "OpenGL error " << err << " occurred.\n";
	}
	
	initialized = true;
}

void AFTrajMovie::subImage() const
{ 	
	assert(initialized);
	bind();
	glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height,
			GL_RGB, GL_FLOAT, pixels); 
}

/* Note that this will compute the texture for the given time whether or not the rat is on this
*  trajectory at the given time. This is desirable for various reasons. For example, when
*  we are jumping around in time (not moving continuously through time, computing every
*  time step) the texture needs to be updated b/c it may have changed between the current 
*  time and the last computed time.
*/
void AFTrajMovie::compute(int time)
{
	// TODO: be smarter about updating the image when the rat's not on this trajectory
	assert(initialized);
	assert(colormap != NULL);
	assert(trajInfo != NULL);

	// determine if active and get the rat's position if necessary
	if (trajInfo->getTrajectory(time) == traj) {
		active = true;
		ratPos = trajInfo->getRatLinPos(time);
	} else {
		active = false;
	}

	computeSpikes(time-SPIKE_EXIST >= 0 ? time-SPIKE_EXIST : 0 ,time);
	
	float posMin = trajInfo->getMinPos();
	float posMax = trajInfo->getMaxPos();
	float phaseMin = trajInfo->getMinPhase();
	float phaseMax = trajInfo->getMaxPhase();
	float posRes = (posMax - posMin) / (width-1);
	float phaseRes = (phaseMax - phaseMin) / (height-1);
	float vMax = trajInfo->approxV2Max();
    int pi= 0;
	
	for (int r=0; r<height; r++) {
		for (int c=0; c<width; c++) {
			float pos = posMin + posRes*c;
			float phase = phaseMin + phaseRes*r;
			float v = trajInfo->evalPosPhase(time,pos,phase,traj)/vMax;
			if (v>1.0f) v=1.0f; // clamp to range [0,1]
 
			// color the pixel using the colormap
			int ci = (int)(v*(colormapLength-1.0)); 
			pixels[pi].r = colormap[ci].r;
			pixels[pi].g = colormap[ci].g;
			pixels[pi].b = colormap[ci].b;
            pi++;
		}
	}

	// update the texture
	subImage();
}

void AFTrajMovie::computeSpikes(int t1, int t2)
{
	// do nothing if animal wasn't on this traj at t1 or t2
	if (trajInfo->getTrajectory(t1)!=traj 
			&& trajInfo->getTrajectory(t2)!=traj) 
		return;

	// erase the vector's contents
	spikes.erase(spikes.begin(),spikes.end());

	// find the index into TData::spiketimes of first spike in timerange [t1,t2]
	int i = trajInfo->getSpikeIndex(t1);
	const TData *data = trajInfo->getData();
	int endspike= data->getEndspike();

	double endTime = data->timearray[t2];
	
	for (;i<=endspike;i++) {
		if (data->getSpikeTraj(i) != traj) continue;	
			
		float t = data->getSpikeTime(i);
		if (t > endTime) break;

		// add this spike to our list of spikes to draw
		spikes.push_back(vertex2f(data->getSpikePos(i),data->getSpikePhase(i)));
	}
}

void AFTrajMovie::draw()
{
	// draw movie texture
	glEnable(GL_TEXTURE_2D);
	this->bind();

	glBegin(GL_QUADS);
	glTexCoord2f(0.0f,0.0f);
	glVertex3f(TX_L,TX_B,0);
	glTexCoord2f(0.0f,1.0f);
	glVertex3f(TX_L,TX_T,0);
	glTexCoord2f(1.0f,1.0f); 
	glVertex3f(TX_R,TX_T,0);
	glTexCoord2f(1.0f,0.0f); 
	glVertex3f(TX_R,TX_B,0);
	glEnd();

	glDisable(GL_TEXTURE_2D);
	
	glEnable(GL_COLOR_MATERIAL);

	// TODO can I just put this in glInit ? 
	glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE); 
	
	drawPosAxis();	
	drawPhaseAxis();
    drawSpikes();
	
	// display trajectory number
	char temp[256];
	glRasterPos3f(TX_L,TX_T+0.01f,0.0f);
	snprintf(temp,"Trajectory: %d",traj);
	gl_draw(temp);

	glDisable(GL_COLOR_MATERIAL);
}

// draw spikes' in position-phase grid
void AFTrajMovie::drawSpikes()
{
	glPushMatrix();
	glTranslatef(TX_L, TX_B, 0.0f);
	glScalef(TX_R-TX_L,TX_T-TX_B,1.0f);

	// draw "+" to represent a spike
	glColor3f(1.0f,1.0f,1.0f);
	glBegin(GL_LINES);
    float x, y;
	for (vector<vertex2f>::iterator s=spikes.begin(); s<spikes.end(); s++) {
		x = normalize(s->x, trajInfo->getMinPos(), trajInfo->getMaxPos()); 
		y = normalize(s->y, trajInfo->getMinPhase(), trajInfo->getMaxPhase()); 
		glVertex3f(x-0.005, y, 0.0f);
		glVertex3f(x+0.005, y, 0.0f);
		glVertex3f(x, y-0.02, 0.0f);
		glVertex3f(x, y+0.02, 0.0f);
	}
	glEnd();
	glPopMatrix();
}

void AFTrajMovie::drawSpikesPos()
{
	glPushMatrix();
	glTranslatef(TX_L, TX_B, 0.0f);
	glScalef(TX_R-TX_L,1.0f,1.0f);

	glColor3f(1.0f,1.0f,1.0f);
	glBegin(GL_LINES);
	for (vector<vertex2f>::iterator s=spikes.begin(); s<spikes.end(); s++) {
		float _x = normalize(s->x, trajInfo->getMinPos(), trajInfo->getMaxPos());
		glVertex3f(_x, -0.03f, 0.0f);
		glVertex3f(_x,  0.14f, 0.0f);
	}
	glEnd();
	glPopMatrix();
}

void AFTrajMovie::drawPosAxis()
{
	glPushMatrix();
	glTranslatef(TX_L, TX_B, 0.0f);
	glScalef(TX_R-TX_L,1.0f,1.0f);

	glColor3f(0.7f,0.7f,1.0f);
	drawAxis(trajInfo->getMinPos(), trajInfo->getMaxPos(),
			10.0f, 50.0f, 0.02f, 0.07f, 0.1f);

    // draw white bar at rat's current location
	if (active) {
		float _x = normalize(ratPos, trajInfo->getMinPos(), trajInfo->getMaxPos()); 
		float l = _x-0.02f, r = _x+0.02f, b = 0.0f, t = 0.03f;
		glColor3f(1.0f,1.0f,0.0f);
		glBegin(GL_QUADS);
		glVertex3f(l,b,0.0f);
		glVertex3f(l,t,0.0f);
		glVertex3f(r,t,0.0f);
		glVertex3f(r,b,0.0f);
		glEnd();
	}
	glPopMatrix();
}

void AFTrajMovie::drawPhaseAxis()
{
	glPushMatrix();
	glTranslatef(TX_L, TX_B, 0.0f);
	glRotatef(0.0f, 1, 0, 0);
	glRotatef(0.0f, 0, 1, 0);
	glRotatef(90, 0, 0, 1);
	glScalef(TX_T-TX_B,-1.0f,1.0f);

	glColor3f(0.7f,0.7f,1.0f);
	drawAxis(trajInfo->getMinPhase(), trajInfo->getMaxPhase(),
			M_PI_4, M_PI_2, 0.005f, 0.015f, TX_L/1.4f);

	glPopMatrix();
}
