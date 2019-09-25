// $Id: EventHandler.h,v 1.1 2008/08/24 19:47:09 chengs Exp $
#ifndef __EventHandler_H__
#define __EventHandler_H__

#include <vector>
#include "eventListeners.h"

using namespace std;

namespace AFilter {

class EventHandler 
{
    private:
        int nTraj;
    private:
        // event listeners
        vector<TrajMeshEventListener*> *trajMeshEventListeners; // an nTraj element array
        vector<MovieEventListener*> movieEventListeners;

        vector<MTimeEventListener*> mTimeEventListeners;
    public:
        EventHandler(int n);
        ~EventHandler();

        /* add and remove event listeners */
        void addTrajMeshEventListener(int traj, TrajMeshEventListener* l);
        void removeTrajMeshEventListener(int traj, TrajMeshEventListener *l);
        void addMovieEventListener(MovieEventListener *l);
        void removeMovieEventListener(MovieEventListener *l);

        // route these to MovieContext object
        void addMTimeEventListener(MTimeEventListener *l) ;
        //        { movieContext.addMTimeEventListener(l); }
        void removeMTimeEventListener(MTimeEventListener *l) ;
        //        { movieContext.removeMTimeEventListener(l); }


        void notifyTrajMeshEventListeners(int traj);
        void notifyMovieEventListeners();
        void notifyMTimeEventListeners();

}; // class

}; // namespace

using namespace AFilter;

#endif
