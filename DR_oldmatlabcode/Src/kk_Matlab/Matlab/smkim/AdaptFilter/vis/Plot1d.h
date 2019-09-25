// $Id: Plot1d.h,v 1.1 2008/08/24 19:47:09 chengs Exp $
#ifndef AF_Plot1d_h
#define AF_Plot1d_h

// OpenGL includes
#include <FL/gl.h>
//#include <GL/glu.h>
#include <assert.h>

// vertex_tmesh, vertex3f, ManageMatrix etc. are defined in here
#include "glUtils.h"
#include "VPlot.h"

using namespace std;

namespace AFilter {

class VDynamics1d;

/** Plot one dimensional model.
 */
class Plot1d : public VPlot {
private:
    bool initialized;
    VDynamics1d *fct;

    int width;
    int height;

    int nTraj;
    int nVerts;
    float **vertices;  // [traj][(x,y) X nVerts]

    float xMin, xMax, xRes;
    float yMin, yMax;

    bool logx, logy;

public:
    // the coordinates of the texture rect's corners
    static const float TX_L;
    static const float TX_R;
    static const float TX_B;
    static const float TX_T;

    // the offset of the axes from the texture
    static const float AO_L;
    static const float AO_B;

public:
    Plot1d(VDynamics1d *fct);
    virtual ~Plot1d();

    /*
     * The following require a current OpenGL context.
     */
    void init();
    void draw();
    void compute(int time);

    void setLogX(bool val= true) { logx= val; }
    void setLogY(bool val= true) { logy= val; }

    bool getLogX() const { return logx; }
    bool getLogY() const { return logy; }

    float getXMin() const { return xMin; }
    float getXMax() const { return xMax; }
    float getYMin() const { return yMin; }
    float getYMax() const { return yMax; }
private:
    void createVertexArray();
    void drawXAxis(int traj);
    void drawYAxis(int traj);

};

}; // namespace AFilter 

#endif // #ifdef AF_Plot1d_h
