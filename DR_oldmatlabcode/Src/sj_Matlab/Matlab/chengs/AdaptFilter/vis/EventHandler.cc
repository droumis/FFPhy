// $Id: EventHandler.cc,v 1.1 2008/08/24 19:47:08 chengs Exp $
#include <assert.h>
#include "EventHandler.h"


EventHandler::EventHandler(int n)
    : nTraj(n)
{
    trajMeshEventListeners= new vector<TrajMeshEventListener*>[nTraj];
}

EventHandler::~EventHandler()
{
    delete[] trajMeshEventListeners;
}

void EventHandler::addTrajMeshEventListener(int traj, TrajMeshEventListener *l)
{
	assert(traj < nTraj);

	/* We need to compute the mesh if it isn't up to date. It won't be up 
	*  to date if the trajMeshEventListenets vector was previously empty
	*/
//@@    if (trajMeshEventListeners[traj].empty()) meshes[traj]->compute(movieContext.getMTime());
	
	trajMeshEventListeners[traj].push_back(l);
}

void EventHandler::removeTrajMeshEventListener(int traj, TrajMeshEventListener *l)
{
	for (vector<TrajMeshEventListener*>::iterator p = trajMeshEventListeners[traj].begin();
			p < trajMeshEventListeners[traj].end();
			p++)
	{
		if (*p == l)
			trajMeshEventListeners[traj].erase(p);
	}
}

void EventHandler::addMovieEventListener(MovieEventListener *l)
{
	movieEventListeners.push_back(l);
}

void EventHandler::removeMovieEventListener(MovieEventListener *l)
{
	for (vector<MovieEventListener*>::iterator p = movieEventListeners.begin();
			p < movieEventListeners.end();
			p++)
	{
		if (*p == l)
			movieEventListeners.erase(p);
	}
}

void EventHandler::notifyTrajMeshEventListeners(int traj)
{	
	for (vector<TrajMeshEventListener*>::iterator p = trajMeshEventListeners[traj].begin();
			p < trajMeshEventListeners[traj].end();
			p++)
	{
		(*p)->trajMeshEventHandle();
	}
}

void EventHandler::notifyMovieEventListeners()
{
		for (vector<MovieEventListener*>::iterator p = movieEventListeners.begin();
			p < movieEventListeners.end();
			p++)
	{
		(*p)->movieEventHandle();
	}
}


void EventHandler::addMTimeEventListener(MTimeEventListener *l)
{
	mTimeEventListeners.push_back(l);
}

void EventHandler::removeMTimeEventListener(MTimeEventListener *l)
{
	for (vector<MTimeEventListener*>::iterator p = mTimeEventListeners.begin();
			p < mTimeEventListeners.end();
			p++)
	{
		if (*p == l)
			mTimeEventListeners.erase(p);
	}
}

void EventHandler::notifyMTimeEventListeners()
{	
	for (vector<MTimeEventListener*>::iterator p = mTimeEventListeners.begin();
			p < mTimeEventListeners.end();
			p++)
	{
		(*p)->mTimeEventHandle();
	}
}
