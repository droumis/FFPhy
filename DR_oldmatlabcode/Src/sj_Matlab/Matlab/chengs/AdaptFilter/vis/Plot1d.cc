// $Id: Plot1d.cc,v 1.1 2008/08/24 19:47:09 chengs Exp $
//
#include "../model/VDynamics1d.h"
#include "axes.h"
#include "Plot1d.h"

using namespace AFilter;

// the coordinates of the texture rect's corners
const float Plot1d::TX_L = 0.06f;
const float Plot1d::TX_R = 0.97f;
const float Plot1d::TX_B = 0.12f;
const float Plot1d::TX_T = 0.90f;

// the offset of the axes from the texture
const float Plot1d::AO_L = 0.03f;
const float Plot1d::AO_B = 0.08f;

Plot1d::Plot1d(VDynamics1d *f)
: fct(f), logx(false), logy(false)
{
    initialized= false;
    nVerts = 200;
    nTraj= fct->getNFct();
	
    vertices = new float*[nTraj];
    for(int n= 0; n < nTraj; n++) vertices[n]= new float[nVerts*2];
}

Plot1d::~Plot1d()
{
    for(int n= 0; n < nTraj; n++) delete[] vertices[n];
    delete[] vertices;
    vertices= 0;
}

void Plot1d::init()
{
    if(initialized) return;
    /* curve vertices have x,y coordinates in range [0,1]
     *  z coordinate is 0
     */
    const double small= 1e-5;
    xMin = (float) fct->getMinX();
    xMax = (float) fct->getMaxX();
    if(logx) {
        xMin= xMin > tiny ? log10(xMin) : log10(tiny);
        xMax= xMax > tiny ? log10(xMax) : log10(tiny);
    }
    xMin+= small;
    xMax-= small;

    yMin = (float) fct->getMinZ();
    yMax = (float) fct->getMaxZ();
    hippo_Print(yMin);
    hippo_Print(yMax);
    yMin= (yMin < 0) ? 1.2f*yMin : yMin/1.2f;
    yMax= (yMax < 0) ? yMax/1.2f : yMax*1.2f;

    xRes = (xMax-xMin)/(nVerts-1)-small;

    assert(xMin+xRes*(nVerts-1) < xMax);

    
    // init x and z with permanent coordinates
    // init y initially with zig-zag line
    float xStep= 1.0f/(float)nVerts;
    for(int n= 0; n < nTraj; n++) {
        for (int i=0; i<nVerts; i++) {
            vertices[n][2*i] = (float) xStep*i;
            vertices[n][2*i+1] = (i%2==0) ? 0.0f : 1.0f;
        }

    }
    //        for(int j=0; j<2*nVerts; j++) printf("%4.2f\t", vertices[0][j]);
    //        printf("\n");
    initialized= true;
}


void Plot1d::compute(int time)
{
    if(!logx) 
        for(int traj= 0; traj < nTraj; traj++) {
            float x = xMin;
            for (int i=0; i<nVerts; i++) {
                vertices[traj][2*i+1]= (fct->evalAtTime(time, x, traj)-yMin)/(yMax-yMin);	
                x+=xRes;
            }
        }
    else if(!logy)
        for(int traj= 0; traj < nTraj; traj++) {
            float x = xMin;
            for (int i=0; i<nVerts; i++) {
                vertices[traj][2*i+1] = (fct->evalAtTime(time, pow(10.0f,x), traj)-yMin)/(yMax-yMin);	
                x+=xRes;
            }
        }
    else
        for(int traj= 0; traj < nTraj; traj++) {
            float x = xMin;
            for (int i=0; i<nVerts; i++) {
                vertices[traj][2*i+1] = (log10(fct->evalAtTime(time, pow(10.0f,x), traj))-log10(yMax))/(log10(yMax)-log10(yMin));	
                x+=xRes;
            }
        }
}

void Plot1d::draw()
{
    assert(initialized);
    for (int traj=0; traj<nTraj; traj++) {
        glPushMatrix();
        glTranslatef(TX_L,(float)TX_B-traj,0.0f);
        glScalef(TX_R-TX_L,TX_T-TX_B,1.0f);

        glDisable(GL_LIGHTING);
        glEnable(GL_COLOR_MATERIAL);
        glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
//        glColor3f(1,0,0);
//        fillPlotSpace();

        if(nTraj== 1)
            glColor3f(0.2,0.2,0.7);
        else if (traj % 2)
            glColor3f(0.7,0.7,0.2);
        else
            glColor3f(0.7,0.2,0.2);

        glLineWidth(4.0f);

        glEnableClientState(GL_VERTEX_ARRAY);

//        for(int j=0; j<2*nVerts; j++) printf("%.2f\t", vertices[0][j]);
//        printf("\n");

        glVertexPointer(2,GL_FLOAT,0,vertices[traj]);
        glDrawArrays(GL_LINE_STRIP,0,nVerts);

        if(nTraj > 1) {
            glColor3f(1.0,1.0,1.0);
            // display trajectory number
            char temp[256];
            glRasterPos3f(0.02f,0.85f,0.0f);
            snprintf(temp, 256, "Trajectory: %d",traj);
            gl_draw(temp);
        }

//        glPushMatrix();
//        glTranslatef(0.0f,1.0f-traj,0.0f);
//        glScalef(1.0f,-1.0f,1.0f);
        glPopMatrix();

        drawXAxis(traj);
        drawYAxis(traj);

    }

}

void Plot1d::drawXAxis(int traj)
{
	glPushMatrix();
	glTranslatef(TX_L, TX_B-traj, 0.0f);
	glScalef(TX_R-TX_L,1.0f,1.0f);

    glColor3f(0.7f,0.7f,1.0f);
    drawAxisNew(xMin, xMax, (xMax-xMin)/10, (xMax-xMin)/5,
            -0.02f,-0.07f,-0.25,-0.04, logx);
//    hippo_Print(xMin);
//    hippo_Print(xMax);
    glPopMatrix();
}

void Plot1d::drawYAxis(int traj)
{
    glPushMatrix();
    glTranslatef(TX_L, TX_B-traj, 0.0f);
    glRotatef(0.0f, 1, 0, 0);
    glRotatef(0.0f, 0, 1, 0);
    glRotatef(90, 0, 0, 1);
    glScalef(TX_T-TX_B,-1.0f,1.0f);

    glColor3f(0.7f,0.7f,1.0f);
    drawAxisNew(yMin, yMax, (yMax-yMin)/10, (yMax-yMin)/5,
            -0.007f,-0.014f,-0.07f,-0.007f,logy);

//    hippo_Print(yMin);
//    hippo_Print(yMax);

    glPopMatrix();
}


