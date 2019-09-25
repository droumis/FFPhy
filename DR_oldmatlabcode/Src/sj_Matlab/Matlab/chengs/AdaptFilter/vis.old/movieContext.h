// $Id: movieContext.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#ifndef __MOVIECONTEXT_H__
#define __MOVIECONTEXT_H__

#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#include <vector>

#include "eventListeners.h"

class AFDataModel;

using namespace std;

/** A state machine responsible for computing the current time for the movie and
 * handling play, stop, rew, etc. messages.
 */
class MovieContext {
public:
    MovieContext();
	
    void setDataModel(AFDataModel *dataModel);

    /************* Movie beginning, ending, etc. *************/
    void play();
    void pause();
    void stop();
    void rew();
    void ff();

    /************* Movie parameters *************/
    void setPlayBegin(int mt) { 
        this->playBeginMTime = mt; gettimeofday(&this->playBeginPTime,NULL); 
        notifyMTimeEventListeners();
    }

    void updateTime(); // compute the mTime when the movie is playing

    /* In the next "set" functions we call setPlayBegin.
     * This is so that these functions can be called while the movie is playing.
     * When the movie is playing we have to do more than set mTime - we have to
     * set the playBegin times as well. 
     */
    void setMTime(int t) { this->mTime = t; setPlayBegin(t); }
    int getMTime() { return mTime; }
    void setPlayRate(float r) { this->playRate = r; setPlayBegin(mTime); }
    float getPlayRate() { return playRate; }

    /* add and remove event listeners */
    void addMTimeEventListener(MTimeEventListener *l);
    void removeMTimeEventListener(MTimeEventListener *l);

private:
    void _play(); // utility function used for play, rew, ff

    /************* Playback Times *************/
    /* There are two types of time: playback time and model time. */
    bool moviePlaying;
    timeval playBeginPTime; // system time (in milliseconds)  at which play began
    int playBeginMTime; // model time (in timesteps) at which play began
    float playRate; // the ratio: model time/playback time
    int mTime; // model time that is currently evaluated
    int prevMTime; // used to calculate a time delta btw current and previous frame 
    int frameAccum; // count number of frames drawn 
    int frameRate; // the frame rate within the last uiDt microseconds
    int mTimeMin; // min. model time that can be evaluated (0 for now)
    int mTimeMax; // max. model time that can be evaluated 
    /* This converts playback time to model time */
    int time_PtoM(int pt) { return (int)(playRate * pt); }

    /************* GUI updating *************/
    /* We allow some changes to accumulate while playing the movie before 
     *  updating the GUI. 
     *  For example, the time slider does not update every frame.
     */
    // all of these time are in microseconds
    int _dt, _lastDt; // used in doFrame
    int dtAccum; // used in doFrame
    int uiDt;

    /* event listeners */
    vector<MTimeEventListener*> mTimeEventListeners;
    void notifyMTimeEventListeners();
};

#endif
