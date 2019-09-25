/* $Id: Plot2d.h,v 1.1 2008/08/24 19:47:09 chengs Exp $
 *
 * Draw 2-d color density plot.
 */

#ifndef AF_Plot2d_H
#define AF_Plot2d_H

// OpenGL includes
#include <FL/gl.h>
//#include <GL/glu.h>
#include <assert.h>
#include <vector>

// vertex_tmesh, vertex3f, ManageMatrix etc. are defined in here
#include "glUtils.h"
#include "../model/VDynamics2d.h"
#include "VPlot.h"

using namespace std;

namespace AFilter {

/** Plot two dimensional model.
 */
class Plot2d : public VPlot {
private:
    bool initialized;
    VDynamics2d *fct;

    GLuint *texId; // the openGL assigned texture identifier
    color3f **pixels; // the (row-major) image data [traj][i]

    color3f *colormap;
    int colormapLength;

    /* these should be powers of 2 */
    int width;
    int height;

    int nTraj;
    int nVerts;

    float xMin, xMax, xRes;
    float yMin, yMax, yRes;
    float zMin, zMax;

    bool logx, logy;

private:

    // axis drawing
    void drawXAxis(traj_type traj);
    void drawYAxis(traj_type traj);

public:
    Plot2d(VDynamics2d *vd);
    virtual ~Plot2d();

    void setColormap(int cLen, color3f *c);

    /* 
     *  The following require a current OpenGL context 
     *  */
    /* Get a texture ID, set the texture interpolation modes, etc. */
    void init();

    /* Set this as the texture for OpenGL to render. */
    void bind(int n) const { glBindTexture(GL_TEXTURE_2D,texId[n]); }

    /* Call bind and then update the texture */	
    void subImage() const;

    void setLogX(bool val= true) { logx= val; }
    void setLogY(bool val= true) { logy= val; }

    void compute(int time);

    /* draws the texture, axes, and label into the rectangle with corners
     *  (0.0, 0.0) and (1.0, 1.0)
     */
    void draw();


}; // class Plot2d

}; // namespace AFilter 

#endif // #ifdef AF_Plot2d_H
