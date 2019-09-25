// $Id: DataCache.cc,v 1.1 2008/08/24 19:47:08 chengs Exp $


#include "../aux/TData.h"
#include "DataCache.h"

using namespace AFilter;

DataCache::DataCache(TData *d)
    : data(d)
{
	int n=data->ntraj;
	posMinT 	= new lazyVar<double>[n];
	posMaxT 	= new lazyVar<double>[n];
	phaseMaxT 	= new lazyVar<double>[n];
	phaseMinT 	= new lazyVar<double>[n];
}

DataCache::~DataCache()
{
	delete[] posMinT;
	delete[] posMaxT;
	delete[] phaseMinT;
	delete[] phaseMaxT;
}


int DataCache::getSpikeIndex(int timestep) const
{
    double time = data->timearray[timestep];
    int endspike= data->getEndspike();
    for (int i= data->getStartspike(); i<= endspike; i++)
        if (data->spiketimes[i] > time)
            return i;
    return endspike;
}

double DataCache::getMinPos()
{
	if (!posMin.valid) {
		minDouble min;
		forEachPos(min);
		posMin = min.m;
	}
	return posMin;
}

double DataCache::getMinPos(traj_type t)
{
	if (!posMinT[t].valid) {
		minDouble min;
		forEachPos(min);
		posMinT[t] = min.m;		
	}
	return posMinT[t];
}

double DataCache::getMaxPos()
{
	if (!posMax.valid) {
		maxDouble max;
		forEachPos(max);
		posMax = max.m;
	}
	return posMax;
}

double DataCache::getMaxPos(traj_type t) 
{
	if (!posMaxT[t].valid)
	{
		maxDouble max;
		forEachPos(max);
		posMaxT[t] = max.m;
	}
	return posMaxT[t];
}

double DataCache::getMinPhase()
{
	if (!phaseMin.valid) {
		minDouble min;
		forEachPhase(min);
		phaseMin = min.m;
	}
	return phaseMin;
}

double DataCache::getMinPhase(traj_type t)
{
	if (!phaseMinT[t].valid) {
		minDouble min;
		forEachPhase(min);
		phaseMinT[t] = min.m;		
	}
	return phaseMinT[t];
}

double DataCache::getMaxPhase()
{
	if (!phaseMax.valid) {
		maxDouble max;
		forEachPhase(max);
		phaseMax = max.m;
	}
	return phaseMax;
}

double DataCache::getMaxPhase(traj_type t)
{
	if (!phaseMaxT[t].valid) {
		maxDouble max;
		forEachPhase(max);
		phaseMaxT[t] = max.m;
	}
	return phaseMaxT[t];
}

double DataCache::getMinXPos() 
{
	if (!xPosMin.valid) {
        if(data->xpos) {
            minDouble min;
            forEachXPos(min);	 // @@invalid read?
            xPosMin = min.m;
        } else
            xPosMin= 0;
	}
	return xPosMin;
}

double DataCache::getMaxXPos()
{
	if (!xPosMax.valid) {
        if(data->xpos) {
            maxDouble max;
            forEachXPos(max);	 // @@invalid read?
            xPosMax = max.m;
        } else
            xPosMax= 0;
	}
	return xPosMax;
}

double DataCache::getMinYPos()
{
	if (!yPosMin.valid) {
        if(data->ypos) {
            minDouble min;
            forEachYPos(min);		 // @@invalid read?
            yPosMin = min.m;
        } else 
            yPosMin= 0;
	}
	return yPosMin;
}

double DataCache::getMaxYPos()
{
	if (!yPosMax.valid) {
        if(data->ypos) {
            maxDouble max;
            forEachYPos(max);	 // @@invalid read?
            yPosMax = max.m;
        } else
            yPosMax= 0;
	}
	return yPosMax;
}

double DataCache::getMinIsi()
{
	if (!tIsiMin.valid) {
        if(data->isi && data->isi[0]) {
            minDouble min;
            forEachIsi(min);	
            tIsiMin = min.m;
        } else 
            tIsiMin= 0;
	}
	return tIsiMin;
}

double DataCache::getMaxIsi()
{
	if (!tIsiMax.valid) {
        if(data->isi && data->isi[0]) {
            maxDouble max;
            forEachIsi(max);	
            tIsiMax = max.m;
        } else 
            tIsiMax= 0;
	}
	return tIsiMax;
}

