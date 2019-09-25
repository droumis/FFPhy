// $Id: DataCache.h,v 1.1 2008/08/24 19:47:08 chengs Exp $
#ifndef AF_DataCache_h
#define AF_DataCache_h

#include <float.h>
#include <iostream>
#include "../aux/defs.h"
#include "../aux/TData.h"
#include "../aux/numerics.h"

using namespace std;


namespace AFilter {
/* 
*  This class does evaluation of lots of things that require time-consuming
*  iterations through the arrays in TData. Because the computations performed by
*  these functions are expensive, their results are cached. Of course,
*  by caching, this assumes that TData is not going to change!
*/
class DataCache 
{	
	public:
		DataCache(TData *data);
		~DataCache();

		/* get info for a given time */
//        traj_type getTrajectory(int timestep) const
//        { return data->getTrajectory(timestep); }
//        
//        float getRatLinPos(int timestep) const
//        { return data->getPos(timestep); }
//        
//        float getRatPhase(int timestep) const
//        { return data->getPhase(timestep); }
//        
		const TData* getData() { return data; }
	
		// get first spike after time
		// returns data->nspikes if all spikes come before this time
		int getSpikeIndex(int time) const;
		
		/* Approximate the max. value returned by eval 
		*  for the VDynamics1d and VDynamics2d object respectively
		*/
		float approxV1Max();
		float approxV2Max();
		float getV1Min() const { return 0.0f; }
		float getV2Min() const { return 0.0f; }
        // these min and max are only for data, not necessarily for the model!
        int    getMinTimeIndex() const { return 0; };
        int    getMaxTimeIndex() const { return data->getNTimesteps(); };
		
		/* get min, max over all or over a single trajectory */
		double getMinPos();
		double getMinPos(traj_type t);
		double getMaxPos();
		double getMaxPos(traj_type t);
		double getMinPhase();
		double getMinPhase(traj_type t);
		double getMaxPhase();
		double getMaxPhase(traj_type t);	

		/* get min, max X and Y (track coordinate) positions */
		double getMinXPos();
		double getMaxXPos();
		double getMinYPos();
		double getMaxYPos();
		
		/* get min, max isi */
		double getMinIsi();
		double getMaxIsi();
	private:
		lazyVar<double> tIsiMax, tIsiMin;
		lazyVar<float> v1Max, v2Max;
		lazyVar<double> posMax, posMin, phaseMax, phaseMin;
		lazyVar<double> *posMaxT, *posMinT, *phaseMaxT, *phaseMinT;
		lazyVar<double> xPosMin, xPosMax, yPosMin, yPosMax;
		
	private:
		TData *data;

		template<typename T> void forEachPos(T& t);
		template<typename T> void forEachPhase(T& t);		
		template<typename T> void forEachXPos(T& t);		
		template<typename T> void forEachYPos(T& t);
		template<typename T> void forEachIsi(T& t);
};

template<typename T> void DataCache::forEachPos(T& t)
{
	for (int i=getMinTimeIndex(); i<=getMaxTimeIndex(); i++)
		t(data->posarray[i]);
}

template<typename T> void DataCache::forEachPhase(T& t)
{
    for (int i=getMinTimeIndex(); i<=getMaxTimeIndex(); i++)
        t(data->phasearray[i]);
}

template<typename T> void DataCache::forEachXPos(T& t)
{
    for (int i=getMinTimeIndex(); i<=getMaxTimeIndex(); i++)
        t(data->xpos[i]);
}

template<typename T> void DataCache::forEachYPos(T& t)
{
    for (int i=getMinTimeIndex(); i<=getMaxTimeIndex(); i++)
        t(data->ypos[i]);
}

template<typename T> void DataCache::forEachIsi(T& t)
{
    for (int i=getMinTimeIndex(); i<=getMaxTimeIndex(); i++)
        t(data->isi[i]);
}


struct minDouble {
    double m;
    minDouble() { m=DBL_MAX; }
    void operator()(double d) { m=(d<m?d:m); } 	
};

struct maxDouble {
    double m;
    maxDouble() { m=DBL_MIN; }
    void operator()(double d) { m=(d>m?d:m); }
};

}; // namespace AFilter

#endif // AF_DataCache_h
