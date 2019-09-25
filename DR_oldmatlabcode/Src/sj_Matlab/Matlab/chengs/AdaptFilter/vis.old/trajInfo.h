// $Id: trajInfo.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#ifndef TRAJ_INFO_H
#define TRAJ_INFO_H

#include <float.h>
#include <iostream>
#include "../aux/defs.h"
#include "../model/PosPhase_Isi.h"
#include "../aux/TData.h"
#include "../aux/numerics.h"

using namespace std;
using namespace AFilter;

/* 
*  This class does evaluation of lots of things that require time-consuming
*  iterations through the arrays in TData. Because the computations performed by
*  these functions are expensive, their results are cached. Of course,
*  by caching, this assumes that TData is not going to change!
*/
class TrajInfo 
{	
	public:
		TrajInfo(TData *data, PosPhase_Isi *am);
		~TrajInfo();

		double evalIsi(int time, double isi, traj_type traj)
		{ 
            //@@ hack, this should really be done in PosPhase_Isi
            if(isi >= tIsiMin && isi <= tIsiMax)
                return vd1d->eval(time,isi,traj); 
            else
                return 1;
        }
		
		double evalPosPhase(int time, double pos, double phase, traj_type traj)
		{ return vd2d->eval(time,pos,phase,traj); }

		/* get info for a given time */
		traj_type getTrajectory(int timestep) const
		{ return data->getTrajectory(timestep); }
		
		float getRatLinPos(int timestep) const
		{ return data->getPos(timestep); }
		
		float getRatPhase(int timestep) const
		{ return data->getPhase(timestep); }
		
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
        int    getMinTimeIndex() const;
        int    getMaxTimeIndex() const;
		
		/* get min, max over all or over a single trajectory */
		double getMinTIsi();
		double getMaxTIsi();
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
		
	private:
		lazyVar<double> tIsiMax, tIsiMin;
		lazyVar<float> v1Max, v2Max;
		lazyVar<double> posMax, posMin, phaseMax, phaseMin;
		lazyVar<double> *posMaxT, *posMinT, *phaseMaxT, *phaseMinT;
		lazyVar<double> xPosMin, xPosMax, yPosMin, yPosMax;
		
	private:
		TData *data;
		VDynamics1d *vd1d;
		VDynamics2d *vd2d;
        PosPhase_Isi *amodel;

		template<typename T> void forEachPos(T& t);
		template<typename T> void forEachPhase(T& t);		
		template<typename T> void forEachXPos(T& t);		
		template<typename T> void forEachYPos(T& t);
};

template<typename T> void TrajInfo::forEachPos(T& t)
{
	for (int i=getMinTimeIndex(); i<=getMaxTimeIndex(); i++)
		t(data->posarray[i]);
}

template<typename T> void TrajInfo::forEachPhase(T& t)
{
    for (int i=getMinTimeIndex(); i<=getMaxTimeIndex(); i++)
        t(data->phasearray[i]);
}

template<typename T> void TrajInfo::forEachXPos(T& t)
{
    for (int i=getMinTimeIndex(); i<=getMaxTimeIndex(); i++)
        t(data->xpos[i]);
}

template<typename T> void TrajInfo::forEachYPos(T& t)
{
    for (int i=getMinTimeIndex(); i<=getMaxTimeIndex(); i++)
        t(data->ypos[i]);
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

#endif
