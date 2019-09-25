/* $Id: ModelVis.cc,v 1.2 2008/09/01 18:24:23 chengs Exp $
   authors     : Sen Cheng
   created     : 2005/02/20
 */

#include <FL/Fl.H>

#include "../model/modelFactory.h"
#include "../model/IsiModel.h"
#include "EventHandler.h"
#include "DataCache.h"
#include "Plot2d.h"
#include "Plot1d.h"

#include "ModelVis.h"

using namespace AFilter;

ModelVis::ModelVis(TData *d) 
{
    initialized= false;
    data= d;
    colormap= 0;
    colormapLength= 0;
    model =0;   
    nTraj= data->getNTraj();

    trajPlot= 0;
    isiPlot =0; 
    geomPlot= 0;
    actPlot= 0;

    cache= new AFilter::DataCache(data);
    events= new EventHandler(nTraj);
};

ModelVis::~ModelVis() 
{
    if(cache) delete cache;
    if(events) delete events;
    if(colormap) delete[] colormap;

    if(actPlot) delete actPlot;
    if(geomPlot) delete geomPlot;
    if(trajPlot) delete trajPlot;
    if(isiPlot) delete isiPlot;
    if(model) delete model;
};

void ModelVis::setColormap(const mxArray *c)
{
	double *cd = mxGetPr(c); 
	this->colormapLength = mxGetNumberOfElements(c)/3; 

	this->colormap = new color3f[colormapLength];
	// cerr << "Building colormap of length " << colormapLength << "\n";
	for (int i=0; i<colormapLength; i++) {
		// Matlab arrays are row major
		colormap[i].r = cd[i];
		colormap[i].g = cd[i + colormapLength];
		colormap[i].b = cd[i + 2*colormapLength];

		// cerr << "colormap[" << i << "] = (" << colormap[i].r << ", " << colormap[i].g << ", " << colormap[i].b << ")\n";
	}

	colormapTexture.setColormap(colormapLength, colormap);
}


void ModelVis::glInit()
{
    if (initialized) return;
	
    VDynamics *dyn= model->getSpatialModel();
    int dim= dyn->getDim();
    Plot2d *ip;
    switch (dim) {
        case 0:
            trajPlot= new Plot1d(model->get1dSpatialModel());
            trajPlot->init();
            break;
        case 1:
            trajPlot= new Plot1d(model->get1dSpatialModel());
            trajPlot->init();
            break;
        case 2:
            ip= new Plot2d(model->get2dSpatialModel());
            ip->init();
            ip->setColormap(getColormapLength(),getColormap());
            trajPlot= ip;
            break;
        default:
            hippo_Error("Do not know how to plot dynamics model ", 
                    dyn->getName());
    }

	// allocate objects that represent the model
    isiPlot= new Plot1d(model->getIsiModel());

	// initialize all visible objects
    isiPlot->setLogX();
    isiPlot->init();

    /* construct geometry */
//    cerr << "Constructing geometry plot.\n";
	geomPlot = new Geometry(cache);
	geomPlot->init();
    actPlot= new Activity(cache, isiPlot);
    actPlot->init();
    colormapTexture.init();

    initialized= true;
}

void ModelVis::show()
{
    close();
}


void ModelVis::plotTraj(int time) { trajPlot->compute(time); trajPlot->draw();};

void ModelVis::plotGeometry(int t)
{
    geomPlot->compute(t);
    geomPlot->draw();
}


void ModelVis::plotTrajActivity(int time, int dir) 
    { 
        actPlot->compute(time, dir); actPlot->drawTraj(); 
    }

void ModelVis::plotIsiActivity(int time, int dir) 
    { 
        actPlot->compute(time, dir); actPlot->drawIsi(); 
    }

void ModelVis::plotIsi(int time) 
    { 
        isiPlot->compute(time); isiPlot->draw(); 
    }

void ModelVis::setMexInput(const mxArray *model_param, const mxArray *opts)
{
    model= dynamic_cast<IsiModel*> ( allocModel(model_param, data) );

    hippo_Assert(model, "Can only plot models derived from IsiModel.");

    const mxArray *cmap= mxGetField(opts,0,"cmap");
    hippo_Assert(cmap, "cmap not defined in options");
    setColormap(cmap);
}
