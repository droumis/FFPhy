// $Id: Plot2d.cc,v 1.1 2008/08/24 19:47:09 chengs Exp $
#include "Plot2d.h"
#include "Plot1d.h"
#include "axes.h"
#include "../aux/numerics.h"
#include "../aux/hippoIO.h"

#include <iostream>


using namespace AFilter;

Plot2d::Plot2d(VDynamics2d *vd)
:   fct(vd)
{
    initialized= false;
    nTraj= fct->getNFct();
	
	width = 128;
	height = 64;

    texId= new GLuint[nTraj];
    pixels= new color3f*[nTraj];
	/* setup image data and OpenGL texture */
    for(int n= 0; n < nTraj; n++) { 
        pixels[n] = new color3f[width * height];	

//	cerr << "stride of pixels is " << &pixels[1] << " - " << &pixels[0] << " = " 
//		<< (char*)&pixels[1] - (char*)&pixels[0] << "\n";
//	cerr << "sizeof ( color3f ) = " << sizeof(color3f) << "\n";

        memset(pixels[n], 0, width*height*sizeof(color3f)); // zero image memory
        texId[n]= 0;
    }

	initialized = false;
}

Plot2d::~Plot2d()
{
    for(int n= 0; n < nTraj; n++) delete[] pixels[n];
	delete[] pixels;
}

void Plot2d::setColormap(int cLen, color3f *c)
{
	this->colormapLength = cLen;
	this->colormap = c;
}

void Plot2d::init()
{
    if (!initialized) {
        const double small= 1e-5;
        for(int n= 0; n < nTraj; n++) {
            glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
            glGenTextures(1, texId+n);
//            cerr << "Texture for traj " << n << " assigned texId " << texId[n] << "\n";
            bind(n);

            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, 
                    GL_RGB, GL_FLOAT, pixels[n]);

            // int err = glGetError();
            // assert(err!=0);
            // if (err!=0)
            // cerr << "OpenGL error " << err << " occurred.\n";
        }
        xMin = fct->getMinX()+small;
        xMax = fct->getMaxX()-small;
        yMin = fct->getMinY()+small;
        yMax = fct->getMaxY()-small;
        xRes = (xMax - xMin) / (width-1) -small;
        yRes = (yMax - yMin) / (height-1) -small;

        zMin = (float) fct->getMinZ();
        zMax = (float) fct->getMaxZ();
        hippo_Print(zMin);
        hippo_Print(zMax);
        zMin= (zMin < 0) ? 1.2f*zMin : zMin/1.2f;
        zMax= (zMax < 0) ? zMax/1.2f : zMax*1.2f;

        assert(xMin+xRes*(width-1) < xMax);
        assert(yMin+yRes*(height-1) < yMax);
    }

    initialized = true;
}

void Plot2d::subImage() const
{ 	
    assert(initialized);
    for(int n= 0; n < nTraj; n++) {
        bind(n);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height,
                GL_RGB, GL_FLOAT, pixels[n]); 
    }
}

/* Note that this will compute the texture for the given time whether or not the rat is on this
 *  trajectory at the given time. This is desirable for various reasons. For example, when
 *  we are jumping around in time (not moving continuously through time, computing every
 *  time step) the texture needs to be updated b/c it may have changed between the current 
 *  time and the last computed time.
 */
void Plot2d::compute(int time)
{
//    hippo_Print(time);
    // TODO: be smarter about updating the image when the rat's not on this trajectory
    hippo_Assert(initialized, "Plot2d not initialized");
    hippo_Assert(colormap, "colormap not set");
    hippo_Assert(fct, "function not set");

    for(traj_type traj= 0; traj < nTraj; traj++) {
        int pi= 0;
        for (int r=0; r<height; r++) {
            float y = yMin + yRes*r;
//            if(!(y > yMin && y < yMax))  {
//                printf("%.6f \t%.6f \t%.6f\n",yMin,y,yMax);
//            printf("%.6f \t%.6f \t%.6f\n",fct->getMinY(),y,fct->getMaxY());
//            }
//                printf("y=%.8f\n",y);
            for (int c=0; c<width; c++) {
                float x = xMin + xRes*c;
//                assert(x > xMin && x < xMax);
//            if(!(x > xMin && x < xMax))  {
//                printf("%.6f \t%.6f \t%.6f\n",xMin,x,xMax);
//            }
//            printf("%.6f \t%.6f \t%.6f\n",fct->getMinX(),x,fct->getMaxX());
//                printf("x=%.8f\n",x);
                float v = (fct->evalAtTime(time,x,y,traj)-zMin)/(zMax-zMin);
                if (v<0.0f) v= 0.0f; // clamp to range [0,1]
                if (v>1.0f) v= 1.0f; // clamp to range [0,1]

                // color the pixel using the colormap
                int ci = (int)(v*(colormapLength-1.0)); 
                pixels[traj][pi].r = colormap[ci].r;
                pixels[traj][pi].g = colormap[ci].g;
                pixels[traj][pi].b = colormap[ci].b;
                pi++;
            }
        }
    }

    // update the texture
    subImage();
}

void Plot2d::draw()
{
    // draw movie texture
    for(traj_type traj= 0; traj < nTraj; traj++) {
        glPushMatrix();
        glTranslatef(Plot1d::TX_L,Plot1d::TX_B-traj,0.0f);
        glScalef(Plot1d::TX_R-Plot1d::TX_L,Plot1d::TX_T-Plot1d::TX_B,1.0f);
        glEnable(GL_TEXTURE_2D);
        this->bind(traj);

        glBegin(GL_QUADS);
        glTexCoord2f(0.0f,0.0f);
//        glVertex3f(Plot1d::TX_L,Plot1d::TX_B-traj,0);
        glVertex3f(0,0,0);
        glTexCoord2f(0.0f,1.0f);
//        glVertex3f(Plot1d::TX_L,Plot1d::TX_T-traj,0);
        glVertex3f(0,1,0);
        glTexCoord2f(1.0f,1.0f); 
//        glVertex3f(Plot1d::TX_R,Plot1d::TX_T-traj,0);
        glVertex3f(1,1,0);
        glTexCoord2f(1.0f,0.0f); 
//        glVertex3f(Plot1d::TX_R,Plot1d::TX_B-traj,0);
        glVertex3f(1,0,0);
        glEnd();

        glDisable(GL_TEXTURE_2D);

        glEnable(GL_COLOR_MATERIAL);

        // TODO can I just put this in glInit ? 
        glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE); 

        glColor3f(1.0,1.0,1.0);
        // display trajectory number
        char temp[256];
        glRasterPos3f(0.0f,1.05f,0.0f);
        snprintf(temp, 256, "Trajectory: %d",traj);
        gl_draw(temp);

        glDisable(GL_COLOR_MATERIAL);
        glPopMatrix();

        drawXAxis(traj);	
        drawYAxis(traj);
    }
}


void Plot2d::drawXAxis(traj_type traj)
{
    glPushMatrix();
    glTranslatef(Plot1d::TX_L, Plot1d::TX_B-traj, 0.0f);
    glScalef(Plot1d::TX_R-Plot1d::TX_L,1.0f,1.0f);

    glColor3f(0.7f,0.7f,1.0f);
    drawAxis(xMin, xMax, 10.0f, 50.0f, 0.02f, 0.07f, 0.15f);

    glPopMatrix();
}

void Plot2d::drawYAxis(traj_type traj)
{
    glPushMatrix();
    glTranslatef(Plot1d::TX_L, Plot1d::TX_B-traj, 0.0f);
    glRotatef(0.0f, 1, 0, 0);
    glRotatef(0.0f, 0, 1, 0);
    glRotatef(90, 0, 0, 1);
    glScalef(Plot1d::TX_T-Plot1d::TX_B,-1.0f,1.0f);

    glColor3f(0.7f,0.7f,1.0f);
    drawAxis(yMin, yMax, M_PI_4, M_PI_2, 0.005f, 0.015f, Plot1d::TX_L/1.4f);

    glPopMatrix();
}
