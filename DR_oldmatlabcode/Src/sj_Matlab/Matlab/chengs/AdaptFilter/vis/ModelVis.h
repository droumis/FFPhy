/* $Id: ModelVis.h,v 1.1 2008/08/24 19:47:09 chengs Exp $
   File name   : vis/ModelVis.h
   authors     : Sen Cheng
   created     : Sun Feb 20 10:57:17 PST 2005
  
 */
#ifndef AF_ModelVis_H
#define AF_ModelVis_H

#include "../aux/defs.h"
#include "../aux/hippoIO.h"

#include <assert.h>
#include <vector>

#include "Geometry.h"
#include "colormapTexture.h"
#include "eventListeners.h"
#include "Activity.h"

using namespace std;


namespace AFilter {

class DataCache;
class EventHandler;
class VPlot;
class Plot1d;

class IsiModel;

/** Show estimated model.
 * 
 * @author Sen Cheng
 * @date   2005/02/20
 */
class ModelVis 
{
private:

protected:
    bool initialized;
    TData *data;
    DataCache *cache;
    IsiModel *model;  //@@ set model 
    int nTraj;

    ColormapTexture colormapTexture;
    color3f *colormap; 
    int colormapLength;

    EventHandler *events;

    VPlot    *trajPlot; //@@ define

    Geometry *geomPlot;
    Activity *actPlot;
    Plot1d   *isiPlot; 
protected:
    inline void setColormap(const mxArray *cmap);

    inline color3f* getColormap() const { return colormap; }
    inline int getColormapLength() const { return colormapLength; }
public:
    ModelVis(TData *d= 0);
    virtual ~ModelVis();

    virtual const char* getName() const { return "ModelVis"; }

    virtual void glInit();

    virtual void show();

    virtual void close() {};

    inline virtual EventHandler* getEventHandler() const { return events;}

    inline TData* getData() const { return data; }

    inline IsiModel* getModel() const { return model; }
    

    virtual void plotTraj(int time);

    virtual void plotGeometry(int time);

    virtual void plotTrajActivity(int time, int dir);

    virtual void plotIsiActivity(int time, int dir) ;

    virtual void plotIsi(int time);

    virtual void setMexInput(const mxArray *model, const mxArray *opts);

}; // class

}; // namespace

#endif   // AF_ModelVis_H
