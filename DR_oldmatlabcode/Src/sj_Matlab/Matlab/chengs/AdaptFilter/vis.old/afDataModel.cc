// $Id: afDataModel.cc,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#include <matrix.h> // matlab API 

#include "afDataModel.h"
#include "afMeshViewUI.h" // for opening mesh windows

AFDataModel::AFDataModel(TData *data, PosPhase_Isi *am)
{
	initialized = false;

	cerr << "Setting TData and VDynamics objects.\n" ;
	trajInfo = new TrajInfo(data, am);
	nTraj = data->getNTraj();
	cerr << "...with " << nTraj << " trajectories.\n";
	
	/* construct meshes */
	cerr << "Constructing meshes: \t";
	meshes = new (AFTrajMesh*)[nTraj];	
	for (int i=0; i<nTraj; i++) { 
		meshes[i] = new AFTrajMesh(i, trajInfo); 
	}
	cerr << "\n";

	/* construct movies */
	cerr << "Constructing movies: \t";
	movies = new (AFTrajMovie*)[nTraj];
	for (int i=0; i<nTraj; i++) {
		movies[i] = new AFTrajMovie(i, trajInfo);
	}
	cerr << "\n";

	/* construct isi plot */
	cerr << "Constructing isi plot.\n";
	isiPlot = new AFIsiPlot(trajInfo);

	/* construct geometry */
	cerr << "Constructing geometry plot.\n";
	geomPlot = new AFGeomPlot(trajInfo);
	
	/* construct array of mesh event listeners */
	trajMeshEventListeners = new vector<TrajMeshEventListener*>[nTraj];

	movieContext.setDataModel(this);
}

AFDataModel::~AFDataModel()
{
	delete[] meshes;
	delete[] movies;
	delete isiPlot;
	delete geomPlot;
	delete[] trajMeshEventListeners;
}

void AFDataModel::setColormap(const mxArray *c)
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

	for (int i=0; i<nTraj; i++)
		meshes[i]->setColormap(colormapLength, colormap);
	for (int i=0; i<nTraj; i++)
		movies[i]->setColormap(colormapLength, colormap);
	colormapTexture.setColormap(colormapLength, colormap);
}

void AFDataModel::glInit()
{
	if (initialized) return;
	
	// initialize all visible objects
	
	for (int i=0; i<nTraj; i++)
		meshes[i]->init();

	for (int i=0; i<nTraj; i++)
		movies[i]->init();

	isiPlot->init();

	geomPlot->init();

	colormapTexture.init();
}

void AFDataModel::doFrame()
{
	/* determine the time for which we need to compute 
	*  the objects.
	*/
	movieContext.updateTime();
	int mTime = movieContext.getMTime();
	
	/* compute objects */
	for (int i=0; i<nTraj; i++)
	{
		if (!trajMeshEventListeners[i].empty()) {
			meshes[i]->compute(mTime);
			notifyTrajMeshEventListeners(i);
		}
		movies[i]->compute(mTime);
	}
	
	geomPlot->compute(mTime);
	
	isiPlot->compute(mTime);
	
	notifyMovieEventListeners();
}

void AFDataModel::compute(int time)
{	
	// TODO
	// there is room for optimization here
	// for (int i=0; i<nTraj; i++) {
		// movies[i]->compute(time);
	// }
	
	// TODO
	// may need to compute other meshes if they are displayed on screen
	// need to go through windows, see which meshes are being displayed
	// meshes[trajInfo->getTrajectory(time)]->computeMesh(time);
}

void AFDataModel::computeAll(int time)
{
	// for (int i=0; i<nTraj; i++) {
		// meshes[i]->computeMesh(time);
		// movies[i]->compute(time);
	// }
}

void AFDataModel::openSingleMeshWindow(int traj)
{
	if (traj >= getData()->getNTraj()) return;

	AFMeshViewUI *mvui = new AFMeshViewUI();
	mvui->setDataModelAndTrajectory(this, traj);

	mvui->show();
}

void AFDataModel::addTrajMeshEventListener(int traj, TrajMeshEventListener *l)
{
	assert(traj < getData()->getNTraj());

	/* We need to compute the mesh if it isn't up to date. It won't be up 
	*  to date if the trajMeshEventListenets vector was previously empty
	*/
	if (trajMeshEventListeners[traj].empty()) meshes[traj]->compute(movieContext.getMTime());
	
	trajMeshEventListeners[traj].push_back(l);
}

void AFDataModel::removeTrajMeshEventListener(int traj, TrajMeshEventListener *l)
{
	for (vector<TrajMeshEventListener*>::iterator p = trajMeshEventListeners[traj].begin();
			p < trajMeshEventListeners[traj].end();
			p++)
	{
		if (*p == l)
			trajMeshEventListeners[traj].erase(p);
	}
}

void AFDataModel::addMovieEventListener(MovieEventListener *l)
{
	movieEventListeners.push_back(l);
}

void AFDataModel::removeMovieEventListener(MovieEventListener *l)
{
	for (vector<MovieEventListener*>::iterator p = movieEventListeners.begin();
			p < movieEventListeners.end();
			p++)
	{
		if (*p == l)
			movieEventListeners.erase(p);
	}
}

void AFDataModel::notifyTrajMeshEventListeners(int traj)
{	
	for (vector<TrajMeshEventListener*>::iterator p = trajMeshEventListeners[traj].begin();
			p < trajMeshEventListeners[traj].end();
			p++)
	{
		(*p)->trajMeshEventHandle();
	}
}

void AFDataModel::notifyMovieEventListeners()
{
		for (vector<MovieEventListener*>::iterator p = movieEventListeners.begin();
			p < movieEventListeners.end();
			p++)
	{
		(*p)->movieEventHandle();
	}
}


