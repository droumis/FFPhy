// $Id: Geometry.h,v 1.1 2008/08/24 19:47:09 chengs Exp $
#ifndef __Geometry_H__
#define __Geometry_H__

#include <FL/gl.h>
//#include <GL/glu.h>
#include <assert.h>

#include "glUtils.h"
#include "DataCache.h"

#include "texture.h"

using namespace std;

namespace AFilter {

class Geometry {

private:
    bool initialized;

    DataCache *cache;

    // every point that the rat visited is drawn on this texture
    Texture posTex;


    vertex2f ratPos;
private: 
    void computePosTex(); // draw the locus of rat pts on the texture
public:
    Geometry(DataCache *cache);
    ~Geometry();

    void init();
    void draw();
    void compute(int timestep);

}; // class
}; // namespace

#endif
