/* $Id: Activity.cc,v 1.2 2008/10/23 21:24:10 chengs Exp $
   authors     : Sen Cheng
   created     : 2005/02/20
 */

#include "Activity.h"
#include "Plot1d.h"
#include "../aux/numerics.h"

using namespace AFilter;

const int Activity::SPIKE_EXIST = 10000;

Activity::Activity(DataCache *d, Plot1d *isiPlot, int dim)
{
    fDim= dim;
    cache= d;
    data= cache->getData();
    nTraj= data->getNTraj();
    for(traj_type traj=0; traj < nTraj; traj++) {
        spikes.pos.push_back(vector<float>() );
        if (fDim==2) 
            spikes.phase.push_back(vector<float>() );
    }
    ratPos= -1.0f;
    currtraj= -1;

    minLogIsi= isiPlot->getXMin();
    maxLogIsi= isiPlot->getXMax();
    logisi= isiPlot->getLogX();
}

Activity::Activity(Activity &)
{
}

Activity::~Activity()
{
    //empty;
}

void Activity::init()
{
    minPos= data->getMinPos();
    maxPos= data->getMaxPos();
}

void Activity::compute(int time, int dir)
{
    if(dir > 0) 
        computeSpikes(time-SPIKE_EXIST >= 0 ? time-SPIKE_EXIST : 0 ,time);
    else  {
        int N= data->getNTimesteps();
        computeSpikes(time, time+SPIKE_EXIST < N  ? time+SPIKE_EXIST : N-1);
    }

    ratPos = (float) data->getPos(time);
    currtraj= data->getTrajectory(time);
    currIsi= (float) data->getIsi(time);
}

void Activity::computeSpikes(int t1, int t2)
{
    traj_type traj;
	// erase the vector's contents
    for(traj=0; traj < nTraj; traj++) {
        spikes.pos[traj].erase(spikes.pos[traj].begin(),spikes.pos[traj].end());
        if (fDim==2)
            spikes.phase[traj].erase(spikes.phase[traj].begin(),spikes.phase[traj].end());
    }
    spikes.isi.erase(spikes.isi.begin(),spikes.isi.end());

	// find the index into TData::spiketimes of first spike in timerange [t1,t2]
	int i = cache->getSpikeIndex(t1);
	int endspike= data->getEndspike();

	double endTime = data->getTime(t2);
	
	for (;i<=endspike;i++) {
        traj= data->getSpikeTraj(i);
        if(traj < 0 | traj >= nTraj) continue;
			
		double t = data->getSpikeTime(i);
		if (t > endTime) break;
		// add this spike to our list of spikes to draw
		spikes.pos[traj].push_back(data->getSpikePos(i));
        if (fDim==2)
            spikes.phase[traj].push_back(data->getSpikePhase(i));
        spikes.isi.push_back(data->getIsi(0+data->getSpikeTimeIndex(i)));
//        cerr << data->getIsi(data->getSpikeTimeIndex(i)) << endl;
	}
}

void Activity::drawSpikesTraj()
{
    for(traj_type traj=0; traj < nTraj; traj++) {
        glPushMatrix();
        glTranslatef(Plot1d::TX_L,Plot1d::TX_B-traj,0.0f);
        glScalef(Plot1d::TX_R-Plot1d::TX_L,Plot1d::TX_T-Plot1d::TX_B,1.0f);

        glColor3f(1.0f,1.0f,1.0f);
//        glBegin(GL_LINES);
//        for (vector<float>::iterator s=spikes.pos[traj].begin(); s<spikes.pos[traj].end(); s++) {
//            float _x = normalize(*s, minPos, maxPos);
//            glVertex3f(_x, -0.03f, 0.0f);
//            glVertex3f(_x,  0.14f, 0.0f);
//        }
//        glEnd();

        // draw spike positions
        glPointSize(5.0f);
        glBegin(GL_POINTS);
        for(unsigned i=0; i<spikes.pos[traj].size(); ++i) {
            float _x = normalize(spikes.pos[traj][i], minPos, maxPos);
            float _y;
            if (fDim==2)
                _y= normalize(spikes.phase[traj][i], 0, 2*M_PI);
            else
                _y= 0;
            glVertex3f(_x,_y,0.0f);
        }
        glEnd();

        glPopMatrix();
    }
}

void Activity::drawSpikesIsi()
{
    glPushMatrix();
    glTranslatef(Plot1d::TX_L, Plot1d::TX_B, 0.0f);
    glScalef(Plot1d::TX_R-Plot1d::TX_L,1.0f,1.0f);

    glColor3f(1.0f,1.0f,1.0f);
    glBegin(GL_LINES);
    for (vector<float>::iterator s=spikes.isi.begin(); s<spikes.isi.end(); s++) {
        float _x = logisi ? normalize(log10(*s), minLogIsi, maxLogIsi):
                normalize(*s, minLogIsi, maxLogIsi);
        if(_x<0 || _x >1) continue;
        glVertex3f(_x, -0.03f, 0.0f);
        glVertex3f(_x,  0.14f, 0.0f);
    }
    glEnd();
    glPopMatrix();
}

void Activity::drawPos()
{
    if(currtraj < 0) return;
    glPushMatrix();
    glTranslatef(Plot1d::TX_L, Plot1d::TX_B-currtraj, 0.0f);
    glScalef(Plot1d::TX_R-Plot1d::TX_L,1.0f,1.0f);

    float _x = normalize(ratPos, minPos, maxPos);
    float l = _x-0.015f, r = _x+0.015f, b = 0.0f, t = 0.05f;
    glColor3f(1.0f,1.0f,1.0f);
    glBegin(GL_QUADS);
    glVertex3f(l,b,0.0f);
    glVertex3f(l,t,0.0f);
    glVertex3f(r,t,0.0f);
    glVertex3f(r,b,0.0f);
    glEnd();
    glPopMatrix();
}

void Activity::drawCurrIsi()
{
    glPushMatrix();
    glTranslatef(Plot1d::TX_L, Plot1d::TX_B, 0.0f);
    glScalef((Plot1d::TX_R-Plot1d::TX_L),Plot1d::TX_T-Plot1d::TX_B,1.0f);

    float _x= logisi ?  normalize(log10(currIsi), minLogIsi, maxLogIsi) :
        normalize(currIsi, minLogIsi, maxLogIsi);
    if(_x > 1) _x= 1;
    float l = _x-0.02f, r = _x+0.02f, b = 0.0f, t = 0.03f;
    glColor3f(1.0f,1.0f,1.0f);
    glBegin(GL_QUADS);
    glVertex3f(l,b,0.0f);
    glVertex3f(l,t,0.0f);
    glVertex3f(r,t,0.0f);
    glVertex3f(r,b,0.0f);
    glEnd();
    glPopMatrix();
}
