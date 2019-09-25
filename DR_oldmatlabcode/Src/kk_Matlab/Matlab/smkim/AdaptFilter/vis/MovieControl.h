// $Id: MovieControl.h,v 1.1 2008/08/24 19:47:09 chengs Exp $
#ifndef __MovieControl_H__
#define __MovieControl_H__

#include "Display.h"
#include "EventHandler.h"

#include "../aux/hippoIO.h"
#include "../aux/defs.h"
#include <png.h>

using namespace std;
using namespace AFilter;

namespace AFilter {

class TData;
class ModelVis;

/** Handles the display of the movie like play, stop, pause, ff, rw, etc.
 *  The actual drawing of the functions is left to methods of ModelVis.
 */

class MovieControl : public Display, public MovieEventListener
{
private:
    ModelVis *model;
    TData *data;
    EventHandler *events;

    int nTraj;

    /************* Playback Times *************/
    /* There are two types of time: playback time and model time. */
    bool moviePlaying;
    timeval playBeginPTime; // system time at which play began
    int playBeginMTime; // model time (in timesteps) at which play began
    double playRate; // current play rate, model time/playback time
                    // [frames/ sec]
    double playRateNormal; // regular play rate
    int mTime; // model time that is currently evaluated
    int mTimeMin; // min. model time that can be evaluated (0 for now)
    int mTimeMax; // max. model time that can be evaluated 

    /************* GUI updating *************/
    double _dt; // [usec]

    /********** dumping output to files ************/
    bool dumpFrames;
    int outputFormat; // 0: none, 1: png, 2: ppm 
    unsigned char *pixels;  // pixels of the current frame
    void (MovieControl::*writeFrameFct)();

    void writeppm();        // get current frame and write ppm to disk

    png_bytep* row_ptr;
    void writepng();        // get current frame and write png to disk


public:

    MovieControl(int x, int y, int w, int h, const char *l);
    ~MovieControl();

    void setModel(ModelVis *d);

    /* OpenGL related fns. */
    void glInit();
    // Main methods for graphing the adaptive filter results.
    // Called automatically from within OpenGL.
    void glDraw();

    /*event handling */
    void movieEventHandle() { this->redraw(); }

    // open a new meshWindow
    void openSpecialWindow(int trajectory) { hippo_Empty;	}

    // movie playing
    void play();
    void pause();
    void stop();
    void rew();
    void ff();

    void quit();

    // display current frame (at mTime)
    void updateFrame();


    /************* Movie parameters *************/
    void setPlayBegin(int mt) { 
        this->playBeginMTime = mt; gettimeofday(&this->playBeginPTime,NULL); 
        events->notifyMTimeEventListeners();
    }

    /* In the next "set" functions we call setPlayBegin.
     * This is so that these functions can be called while the movie is playing.
     * When the movie is playing we have to do more than set mTime - we have to
     * set the playBegin times as well. 
     */
//    void setMTime(int t) { this->mTime = t; setPlayBegin(t); }
    void setMTime(int t);
    int getMTime() const { return mTime; }
    int getMTimeMin() const { return mTimeMin; }
    int getMTimeMax() const { return mTimeMax; }

    void setPlayRate(double r) { this->playRate = r; setPlayBegin(mTime); }
    double getPlayRate() { return playRate; }
    short int getPlayDir() { return playRate> 0 ? 1 : -1; }

    /********** dumping output to .ppm files ************/
    void setDumpFrames(bool b);

    void setMexInput(const mxArray *opts);


    static void cb_movie_idle(void *v);
private:
    void _play(); // utility function used for play, rew, ff


}; // class


}; // namespace

#endif
