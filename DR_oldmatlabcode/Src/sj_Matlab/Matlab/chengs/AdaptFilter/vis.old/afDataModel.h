// $Id: afDataModel.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#ifndef __AFDATAMODEL_H__
#define __AFDATAMODEL_H__

#include <assert.h>
#include <vector>

#include "afTrajMesh.h"
#include "afTrajMovie.h"
#include "afIsiPlot.h"
#include "afGeomPlot.h"
#include "colormapTexture.h"
#include "eventListeners.h"
#include "movieContext.h"

using namespace std;

class AFDataModel
{
	public:
		AFDataModel(TData *data, PosPhase_Isi *am);
		~AFDataModel();

		void setColormap(const mxArray *c);
		
		/* Initialize all of the member objects that need
		*  initialization
		*/
		void glInit();

		void doFrame();
		
		void compute(int time); // clever about computing only things that have changed
		void computeAll(int time); // computes everything
	
		const TData* getData() { return trajInfo->getData(); }

		TrajInfo* getTrajInfo() { return trajInfo; }
		
		int getMTime() { return movieContext.getMTime(); }
		
		AFTrajMesh* getMesh(int traj)
		{ assert(traj<nTraj); return meshes[traj]; }
		
		AFTrajMovie* getTexture(int traj) 
		{ assert(traj<nTraj); return movies[traj]; } 
		
		AFIsiPlot* getIsiPlot() { return isiPlot; }
		
		AFGeomPlot* getGeomPlot() { return geomPlot; }
		
		ColormapTexture* getColormapTexture()
		{ return &colormapTexture; }

		// open a new meshWindow
		void openSingleMeshWindow(int trajectory);	
		
		// movie playing
		void setMovieTime(int t) { movieContext.setMTime(t); doFrame(); }
		void playMovie() { movieContext.play(); }
		void pauseMovie() { movieContext.pause(); }
		void stopMovie() { movieContext.stop(); }
		void rewMovie() { movieContext.rew(); }
		void ffMovie() { movieContext.ff(); }

		/* add and remove event listeners */
		void addTrajMeshEventListener(int traj, TrajMeshEventListener* l);
		void removeTrajMeshEventListener(int traj, TrajMeshEventListener *l);
		void addMovieEventListener(MovieEventListener *l);
		void removeMovieEventListener(MovieEventListener *l);
		// route these to MovieContext object
		void addMTimeEventListener(MTimeEventListener *l) 
		{ movieContext.addMTimeEventListener(l); }
		void removeMTimeEventListener(MTimeEventListener *l) 
		{ movieContext.removeMTimeEventListener(l); }
		
	private:
		bool initialized;
		
		MovieContext movieContext;
		
		TrajInfo *trajInfo; 
		int nTraj; // total number of trajs (data->getNTraj())
		
		AFIsiPlot *isiPlot; // interspike interval 
		AFGeomPlot *geomPlot; // track geometry
		
		/* these both have nTraj elems */
		AFTrajMesh **meshes;	// array of ptrs to the trajectories' meshes
		AFTrajMovie **movies;	// array of ptrs to the trajectories' movies

		ColormapTexture colormapTexture;

		color3f *colormap;
		int colormapLength;

		// event listeners
		vector<TrajMeshEventListener*> *trajMeshEventListeners; // an nTraj element array
		vector<MovieEventListener*> movieEventListeners;

		void notifyTrajMeshEventListeners(int traj);
		void notifyMovieEventListeners();
};

#endif
