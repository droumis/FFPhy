// $Id: movieContext.cc,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#include "movieContext.h"
#include "globals.h" // for movie_idle
#include <FL/Fl.H>

const double PLAYRATE_NORMAL = 0.01;
//const double PLAYRATE_FAST = 0.1;

MovieContext::MovieContext()
{
	moviePlaying = false;
	
	playBeginPTime.tv_sec = 0;
	playBeginPTime.tv_usec = 0;
	playBeginMTime = 0;
	mTime = 0;
	prevMTime = 0;
	frameRate = 0;
	frameAccum = 0;
	mTimeMin = 0;
	mTimeMax = 0; 
	playRate = PLAYRATE_NORMAL;
}

void MovieContext::setDataModel(AFDataModel *dataModel)
{ 
    mTimeMin = dataModel->getTrajInfo()->getMinTimeIndex();  
	mTimeMax = dataModel->getTrajInfo()->getMaxTimeIndex();  
}

void MovieContext::play() {
	setPlayRate(PLAYRATE_NORMAL);
	_play();
}

void MovieContext::pause() {
	Fl::remove_idle(movie_idle);
	moviePlaying = false;
}

void MovieContext::stop() {
	pause();
    setMTime(mTimeMin);
}

void MovieContext::ff() 
{
    double r= getPlayRate();
    if (r < 0) r= PLAYRATE_NORMAL;
    else r+= PLAYRATE_NORMAL;
	setPlayRate(r);
	_play();
}

void MovieContext::rew()
{
    double r= getPlayRate();
    if (r > 0) r= -PLAYRATE_NORMAL;
    else r-= PLAYRATE_NORMAL;
	setPlayRate(r);
	_play();
}

void MovieContext::_play()
{
	if (!moviePlaying) {
		// TODO: this probably needs to be changed!!
		Fl::add_idle(movie_idle);
		moviePlaying = true;
	}
}

void MovieContext::updateTime()
{
	timeval t;
	gettimeofday(&t, NULL);

	/* update mTime if the movie is playing */
	if (moviePlaying) {
		_lastDt = _dt;
		_dt = 1000000*(t.tv_sec - playBeginPTime.tv_sec) + (t.tv_usec - playBeginPTime.tv_usec);

		prevMTime = mTime;
		mTime = playBeginMTime + time_PtoM(_dt);

		// accumulate changes then refresh GUI
		dtAccum += (_dt-_lastDt);
		if (dtAccum > uiDt) {
			// cerr << "accumulated mTime changes\n";
			notifyMTimeEventListeners();
		}    
	}

	// handle out-of-bounds mTime
	if (mTime < mTimeMin || mTime > mTimeMax) {
		// Stop the movie before the time gets out of bounds.
		// (in terms of the number of timesteps that the adaptive model can eval)
		if (moviePlaying) {
			stop(); // stop movie 
        }
		if (mTime >= mTimeMax)
            mTime= mTimeMax;
        else
            mTime= mTimeMin;
	}

	// cerr << "mTime = " << mTime <<"\n";
	// cerr << "FLTK "  << (Fl::has_idle(movie_idle) ? "has" : "does not have ") <<  " movie_idle\n"; 
	// cerr << "mTimeMax = " << mTimeMax << "\n";
}

void MovieContext::addMTimeEventListener(MTimeEventListener *l)
{
	mTimeEventListeners.push_back(l);
}

void MovieContext::removeMTimeEventListener(MTimeEventListener *l)
{
	for (vector<MTimeEventListener*>::iterator p = mTimeEventListeners.begin();
			p < mTimeEventListeners.end();
			p++)
	{
		if (*p == l)
			mTimeEventListeners.erase(p);
	}
}

void MovieContext::notifyMTimeEventListeners()
{	
	for (vector<MTimeEventListener*>::iterator p = mTimeEventListeners.begin();
			p < mTimeEventListeners.end();
			p++)
	{
		(*p)->mTimeEventHandle();
	}
}
