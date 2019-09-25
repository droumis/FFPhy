// $Id: trajInfo.cc,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#include "trajInfo.h"

TrajInfo::TrajInfo(TData *d, PosPhase_Isi *am)
    : data(d), amodel(am)
{
    vd1d= amodel->getIsiModel();
    vd2d= amodel->getXPModel();
	int n=data->ntraj;
	posMinT 	= new lazyVar<double>[n];
	posMaxT 	= new lazyVar<double>[n];
	phaseMaxT 	= new lazyVar<double>[n];
	phaseMinT 	= new lazyVar<double>[n];

    // init. min and max isi values
    getMinTIsi();
    getMaxTIsi();
}

TrajInfo::~TrajInfo()
{
	delete[] posMinT;
	delete[] posMaxT;
	delete[] phaseMinT;
	delete[] phaseMaxT;
}

float TrajInfo::approxV1Max()
{
	const int nTime = 100;
	const int nX = 200;
    const float isiStart= 0.001;

	if (!v1Max.valid) {
		int tRes = (getMaxTimeIndex()-getMinTimeIndex()+1)/(nTime-1);	

		getMinTIsi();
		getMaxTIsi();
		float xRes = (log(tIsiMax) - log(isiStart)) / (nX-1);

		v1Max = 0.0f;

		for (int t= getMinTimeIndex(); t <= getMaxTimeIndex(); t+=tRes) {
			for (float x=log(isiStart); x<=log(tIsiMax); x+= xRes) {					
				float v = evalIsi(t, exp(x), 0);
				if (v > v1Max) v1Max = v;
			}
		}
        hippo_Print(v1Max);
	}

	return v1Max;
}

float TrajInfo::approxV2Max()
{
	// number of datapoints to use in our approximation
	const int nTime = 200;
	const int nPos = 100;
	const int nPhase = 40;

	if (!v2Max.valid) {
		getMinPos();
		getMaxPos();
		float posRange = posMax - posMin;
		getMinPhase();
		getMaxPhase();
		float phaseRange = phaseMax - phaseMin; 
		int tRes = (getMaxTimeIndex()-getMinTimeIndex()+1)/(nTime-1);	
		float posRes = posRange/(nPos-1);
		float phaseRes = phaseRange/(nPhase-1);

		v2Max = 0.0f;

		for (int t=getMinTimeIndex(); t <= getMaxTimeIndex(); t+=tRes) {
			for (float pos=posMin; pos<=posMax;pos+=posRes ) {
				for (float phase=phaseMin;phase<=phaseMax; phase+=phaseRes) {
					traj_type ct = getTrajectory(t);
					if (ct!=-1) {
						float v =  evalPosPhase(t,pos,phase,ct);
						if (v > v2Max) v2Max = v;
					}
				}      
			}
		}
        hippo_Print(v2Max);
	}

	return v2Max;
}

int TrajInfo::getMinTimeIndex() const
{
    if(vd1d->getStartindex() != vd2d->getStartindex()) {
        hippo_Print(vd1d->getStartindex() );
        hippo_Print(vd2d->getStartindex() );
    }
    hippo_Assert(vd1d->getStartindex() == vd2d->getStartindex(), 
            "isi and xp model do not have same startindex");
    return vd2d->getStartindex();
}

int TrajInfo::getMaxTimeIndex() const
{
    if(vd1d->getEndindex() != vd2d->getEndindex()) {
        hippo_Print(vd1d->getEndindex() );
        hippo_Print(vd2d->getEndindex() );
    }
    hippo_Assert(vd1d->getEndindex() == vd2d->getEndindex(), 
            "isi and xp model do not have same endindex");
    return vd2d->getEndindex();
}

int TrajInfo::getSpikeIndex(int timestep) const
{
    double time = data->timearray[timestep];
    int endspike= data->getEndspike();
	for (int i= data->getStartspike(); i<= endspike; i++)
		if (data->spiketimes[i] > time)
			return i;
	return endspike;
}

double TrajInfo::getMinTIsi()
{
	if (!tIsiMin.valid)
		tIsiMin = 0.0f;

	return tIsiMin;
}

double TrajInfo::getMaxTIsi()
{
	if (!tIsiMax.valid) {
		tIsiMax = amodel->getMaxIsi();
	}

	return tIsiMax;
}

double TrajInfo::getMinPos()
{
	if (!posMin.valid) {
		minDouble min;
		forEachPos(min);
		posMin = min.m;
	}
	return posMin;
}

double TrajInfo::getMinPos(traj_type t)
{
	if (!posMinT[t].valid) {
		minDouble min;
		forEachPos(min);
		posMinT[t] = min.m;		
	}
	return posMinT[t];
}

double TrajInfo::getMaxPos()
{
	if (!posMax.valid) {
		maxDouble max;
		forEachPos(max);
		posMax = max.m;
	}
	return posMax;
}

double TrajInfo::getMaxPos(traj_type t) 
{
	if (!posMaxT[t].valid)
	{
		maxDouble max;
		forEachPos(max);
		posMaxT[t] = max.m;
	}
	return posMaxT[t];
}

double TrajInfo::getMinPhase()
{
	if (!phaseMin.valid) {
		minDouble min;
		forEachPhase(min);
		phaseMin = min.m;
	}
	return phaseMin;
}

double TrajInfo::getMinPhase(traj_type t)
{
	if (!phaseMinT[t].valid) {
		minDouble min;
		forEachPhase(min);
		phaseMinT[t] = min.m;		
	}
	return phaseMinT[t];
}

double TrajInfo::getMaxPhase()
{
	if (!phaseMax.valid) {
		maxDouble max;
		forEachPhase(max);
		phaseMax = max.m;
	}
	return phaseMax;
}

double TrajInfo::getMaxPhase(traj_type t)
{
	if (!phaseMaxT[t].valid) {
		maxDouble max;
		forEachPhase(max);
		phaseMaxT[t] = max.m;
	}
	return phaseMaxT[t];
}

double TrajInfo::getMinXPos() 
{
	if (!xPosMin.valid) {
        if(data->xpos) {
            minDouble min;
            forEachXPos(min);	
            xPosMin = min.m;
        } else
            xPosMin= 0;
	}
	return xPosMin;
}

double TrajInfo::getMaxXPos()
{
	if (!xPosMax.valid) {
        if(data->xpos) {
            maxDouble max;
            forEachXPos(max);
            xPosMax = max.m;
        } else
            xPosMax= 0;
	}
	return xPosMax;
}

double TrajInfo::getMinYPos()
{
	if (!yPosMin.valid) {
        if(data->ypos) {
            minDouble min;
            forEachYPos(min);	
            yPosMin = min.m;
        } else 
            yPosMin= 0;
	}
	return yPosMin;
}

double TrajInfo::getMaxYPos()
{
	if (!yPosMax.valid) {
        if(data->ypos) {
            maxDouble max;
            forEachYPos(max);
            yPosMax = max.m;
        } else
            yPosMax= 0;
	}
	return yPosMax;
}
