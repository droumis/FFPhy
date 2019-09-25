// $Id: afTrajMesh.cc,v 1.2 2006/11/07 18:02:14 chengs Exp $
#include "afTrajMesh.h"
#include "utils.h"
#include "../aux/numerics.h"
#include "axes.h"
#include <assert.h>

const float AFTrajMesh::AO = 0.1f;
const float AFTrajMesh::LO = 0.1f;
const float AFTrajMesh::VE = 1.3f;

AFTrajMesh::AFTrajMesh(int t, TrajInfo *ti)
{
	cerr << t << "\t";

	initialized = false;
	
	trans.x = -0.5f;		trans.y = 0.0f; 	trans.z = 0.5f;
	rot.x = 0.0f; 		rot.y = 0.0f; 		rot.z = 0.0f;
	scale.x = 1.0f; 	scale.y = 1.0f;		scale.z = -1.0f;

	nVertsX = 0;
	nVertsZ = 0;
	
	colormap = NULL;
	
	setTrajInfo(t, ti);
}

/* This is to be called after the OpenGL context is setup. 
* This function sets up the geometry. traj and trajInfo may be NULL
* when this is called.
*/
void AFTrajMesh::init()
{
	if (!initialized) {
	createVertexArray();
	createMeshIndexArray();
	createGridIndexArray();
	createSpikeArray();
	}
	
	initialized = true;
}

void AFTrajMesh::setTrajInfo(int t, TrajInfo *ti)
{
	this->traj = t;
	this->trajInfo = ti;
}

void AFTrajMesh::setColormap(int cLen, color3f* c)
{
	this->colormapLength = cLen;
	this->colormap = c;
}

AFTrajMesh::~AFTrajMesh() {
	delete[] vertices;
	delete[] indices_mesh;
	delete[] indices_grid;
}

void AFTrajMesh::draw() 
{
	assert(trajInfo!=NULL);
	assert(initialized);

	// glEnableClientState(GL_VERTEX_ARRAY);
	// glEnableClientState(GL_COLOR_ARRAY);
// 
	// glColorPointer(3,GL_FLOAT,sizeof(vertex_tmesh),&vertices[0]);	
	// glVertexPointer(3,GL_FLOAT,sizeof(vertex_tmesh),&vertices[3]);
// 
	// setup OpenGL to render vertex_tmesh arrays
	//	short stride = (short)((char*)(vertices+1) - (char*)(vertices));
	//	int stride = sizeof(vertex_tmesh);
	glInterleavedArrays(GL_C3F_V3F, 0, vertices);
	
	{
		glPushMatrix();
		
		// setup transformation
		glTranslatef(trans.x, trans.y, trans.z);
		glRotatef(rot.y, 0, 1, 0);
		glRotatef(rot.x, 1, 0, 0);
		glRotatef(rot.z, 0, 0, 1);
		glScalef(scale.x, scale.y, scale.z);

		drawMesh();
		// drawGrid();
		drawAxes();

		glEnable(GL_COLOR_MATERIAL);
		drawPosAxisGrid();
		drawPhaseAxisGrid();
		//drawSpikesProjected();
		//drawSpikesOnMesh();

		glPopMatrix();
	}
}

void AFTrajMesh::computeMesh(int time)
{
	if (trajInfo == NULL) return;
	assert(colormap!=NULL);
	
	float posMin = trajInfo->getMinPos();
	float posMax = trajInfo->getMaxPos();
	float phaseMin = trajInfo->getMinPhase();
	float phaseMax = trajInfo->getMaxPhase();
	float vMax = trajInfo->approxV2Max();
	
	for (int i=0; i<nVertsZ*nVertsX; i++) {
		float pos = posMin + vertices[i].x * (posMax - posMin);
		float phase = phaseMin + vertices[i].z * (phaseMax - phaseMin);

		// set y-value
		vertices[i].y = trajInfo->evalPosPhase(time,pos,phase,traj)/vMax;

		// TODO
		// compute normal
    
		// set color
		float v = vertices[i].y;
		v = v<1.0f ? v : 1.0f;
		int ci = (int)(vertices[i].y * (colormapLength-1)); 

		vertices[i].r = (float)colormap[ci].r;
		vertices[i].g = (float)colormap[ci].g;
		vertices[i].b = (float)colormap[ci].b;
	}
}

void AFTrajMesh::computeSpikes(int t1, int t2)
{
	// TODO
}

void AFTrajMesh::drawMesh() 
{
	// glPointSize(5.0f);
	// glDrawArrays(GL_POINTS,0,nVertsZ*nVertsX);
	
	for (int r=0; r<nVertsZ-1; r++) {
		glDrawElements(GL_TRIANGLE_STRIP,2*nVertsX,GL_UNSIGNED_INT,&indices_mesh[r*2*nVertsX]);    
	}
}

/*
void AFTrajMesh::drawGrid()
{
	// this may not work with interleaved arrays
	glDisableClientState(GL_COLOR_ARRAY);

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_DIFFUSE);
	glColor3f(0.0,0.0,1.0);

	// draw grid lines parallel to pos axis
	for (int r=0; r<nVertsZ; r+=4) {
		glDrawElements(GL_LINE_STRIP,nVertsX,GL_UNSIGNED_INT,&indices_grid[r*nVertsX]);    
	}

	// draw grid lines parallel to phase axis
	for (int c=0; c<nVertsX; c+=25) {
		glDrawElements(GL_LINE_STRIP,nVertsZ,GL_UNSIGNED_INT,&indices_grid[nVertsX*nVertsZ+c*nVertsZ]);        
	}  

	glDisable(GL_COLOR_MATERIAL);
}*/

void AFTrajMesh::drawAxes()
{
	assert(trajInfo != NULL);
	
	//	float posMin = trajInfo->getMinPos();
	//	float posMax = trajInfo->getMaxPos();
	//	float phaseMin = trajInfo->getMinPhase();
	//	float phaseMax = trajInfo->getMaxPhase();
	//	float vMax = trajInfo->approxVMax();

	glDisable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
	glColor3f(0.0,1.0,0.0);

	glBegin(GL_LINES);
	// x axis
	glVertex3f(-AO,0.0f,-AO);
	glVertex3f(1.0f,0.0f,-AO);
	// z axis
	glVertex3f(-AO,0.0f,-AO);
	glVertex3f(-AO,0.0f,1.0f);
	// y axis
	glVertex3f(-AO,0.0f,-AO);
	glVertex3f(-AO,VE,-AO); 
	glEnd();

	glRasterPos3f(1.0f,-LO,-AO);
	gl_draw("position"); 

	glRasterPos3f(-AO,-LO,0.0f);
	gl_draw("phase");

	glRasterPos3f(-AO,1.0f,-AO);
	gl_draw("eval");   
}


void AFTrajMesh::drawPosAxisGrid()
{
	// z position at which to locate grid
	float z = -AO;
	
	float xMin = trajInfo->getMinPos();
	float xMax = trajInfo->getMaxPos();

	float yMin = trajInfo->getV2Min();
	float yMax = VE*trajInfo->approxV2Max();
	
	const double XINC_MINOR = 20.0f;
	const double XINC_MAJOR = 100.0f;
	
	double xMinMinor = getCeilingMultiple(xMin, XINC_MINOR);
	double xMinMajor = getCeilingMultiple(xMin, XINC_MAJOR);
	
	const double YINC_MINOR = 5.0f; 
	const double YINC_MAJOR = 10.0f;
	
	double yMinMinor = getCeilingMultiple(yMin, YINC_MINOR);
	double yMinMajor = getCeilingMultiple(yMin, YINC_MAJOR);
	

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
		glVertex3f(_x, 0.0f, z);
		glVertex3f(_x,   VE, z);
		
		xMinor += XINC_MINOR;
	}

	// draw minor grid lines on y-axis
	double yMinor = yMinMinor;
	while (yMinor < yMax)
	{
		double _y = rerange(yMinor, yMin, yMax, 0.0, VE);
		glVertex3f(0.0f, _y, z);
		glVertex3f(1.0f, _y, z); 

		yMinor += YINC_MINOR;
	}

	/* major grid lines */
	
	// set color for major grid lines
	glColor3f(0.2f,0.8f,0.2f);
	
	// draw major grid lines on x-axis
	double xMajor = xMinMajor;
	while (xMajor < xMax)
	{
		double _x = normalize(xMajor, xMin, xMax);
		glVertex3f(_x, 0.0f, z);
		glVertex3f(_x,   VE, z);

		xMajor += XINC_MAJOR;
	}
	
	// draw major grid lines on y-axis
	double yMajor = yMinMajor;
	while (yMajor < yMax)
	{
		double _y = rerange(yMajor, yMin, yMax, 0.0, VE);
		glVertex3f(0.0f, _y, z);
		glVertex3f(1.0f, _y, z);

		yMajor += YINC_MAJOR;
	}
	
	glEnd();

	/* numbering */
	
	char temp[32];
	
	// number major grid lines on x-axis
	xMajor = xMinMajor;
	while (xMajor < xMax)
	{
		double _x = normalize(xMajor, xMin, xMax);
		snprintf(temp,"%.2f",xMajor);
		glRasterPos3f(_x, VE+0.05f, z);
		gl_draw(temp);
		
		xMajor += XINC_MAJOR;
	}

	// number major grid lines on y-axis
	yMajor = yMinMajor;
	while (yMajor < yMax)
	{
		double _y = rerange(yMajor, yMin, yMax, 0.0f, VE);
		snprintf(temp,"%.2f",yMajor);
		glRasterPos3f(1.15f, _y, z);
		gl_draw(temp);
		
		yMajor += YINC_MAJOR;
	}
	
	glDisable(GL_COLOR_MATERIAL);
}

void AFTrajMesh::drawPhaseAxisGrid() 
{
	glPushMatrix();
	glTranslatef(-AO,0.0f,0.0f);
	glRotatef(-90,0,1,0);
	glScalef(1.0f,VE,1.0f);

	drawGrid(trajInfo->getMinPhase(), trajInfo->getMaxPhase(),
			0.7854f, 1.5708f, 0.1f,
			trajInfo->getV2Min(), VE*trajInfo->approxV2Max(),
			5.0f, 10.0f, 0.1f);
	
	glPopMatrix();
}

void AFTrajMesh::drawSpikesOnMesh() 
{
	// TODO
}

void AFTrajMesh::drawSpikesProjected() 
{
	// TODO
}

void AFTrajMesh::createVertexArray()
{
	nVertsX = 160;
	nVertsZ = 50;
	vertices = new vertex_tmesh[nVertsX * nVertsZ];

	// cerr << "stride of vertices is " << &vertices[1] << " - " << &vertices[0] << " = " << (char*)&vertices[1] - (char*)&vertices[0] << "\n";
	// cerr << "sizeof ( vertex_tmesh) = " << sizeof(vertex_tmesh) << "\n";

	/* mesh vertices have x,y coordinates in range [0,1]
	* (z coordinates are normalized into this range, 
	* but because we can only imperfectly compute the range of
	* z, this normalization is only approximate)
	 */
	float xExtent = 1.0f;
	float zExtent = 1.0f;
	float xRes = xExtent/nVertsX;
	float zRes = zExtent/nVertsZ;

	// vertex array is row (strip parallel to x-axis) major
	for (int i=0; i<nVertsX*nVertsZ; i++) {
		//		vertex_tmesh& v = vertices[i];
		// position
		vertices[i].x = xRes*(i%nVertsX);
		vertices[i].y = 0.0f;
		vertices[i].z = zExtent - zRes*(i/nVertsX);
		// color
		vertices[i].r = 1.0f;
		vertices[i].g = 0.0f;
		vertices[i].b = 0.0f;

		//cerr << " vertices[" << i << "] = " << vertices[i] << "\n";
	}
}

/* The mesh index array consists of nVertsZ-1 strips parallel
*  to the x-axis.
*/
void AFTrajMesh::createMeshIndexArray()
{
	indices_mesh = new unsigned int[2*nVertsX*(nVertsZ-1)];

	// we create nVertsZ strips 
	// (this is optimal when nVertsZ < nVertsX)
	for (int r=0; r<nVertsZ-1; r++)
		for (int i=0; i<2*nVertsX; i++) {
			int j = i+2*r*nVertsX;
			if (i%2==0)
				indices_mesh[j]=i/2+r*nVertsX;
			else
				indices_mesh[j]=(i-1)/2+(r+1)*nVertsX;
		}   
}

void AFTrajMesh::createGridIndexArray() 
{
	indices_grid = new unsigned int[2*nVertsX*nVertsZ];

	// create horizontal grid
	for (int r=0; r<nVertsZ; r++)
		for (int c=0; c<nVertsX; c++) {
			int j = r*nVertsX + c;
			indices_grid[j] = j;
		}  

	// create vertical grid
	for (int r=0; r<nVertsZ; r++)
		for (int c=0; c<nVertsX; c++) {
			indices_grid[nVertsX*nVertsZ + c*nVertsZ + r] = r*nVertsX + c; 
		} 
}

void AFTrajMesh::createSpikeArray() {
	nUsedSpikes = 0;
}
