#ifndef TDATA_H
#define TDATA_H

/* $Id: TData.h,v 1.8 2008/08/24 19:46:51 chengs Exp $
   
   Sen Cheng, 2004/../..

   program description
 */

#include "../aux/defs.h"

// struct mxArray_tag;
// typedef struct mxArray_tag mxArray;


/* ************************************************************
   class TData
************************************************************ */
namespace AFilter {

class TData {
public:
    int            ntimesteps;	      // the total number of timesteps
    int		   ntraj;	      // the total number of trajectories, treating each direction separately
    
    traj_type            *fieldID;    // the indicator variable for the number of the trajectory the animal is currently on
    double            *timearray;  // the list of times for estimation
    double            *posarray;   // linearized position [cm]
    double            *xpos;   // x position [cm]
    double            *ypos;   // y position [cm]
    double            *phasearray;   // the theta phase at each time step
    
    double            *x;   // x variable (not necessarily spatial)
    double            *y;   // y variable (not necessarily spatial)

    int            nspikes;	      // the total number of spikes
    double            *spiketimes;  //  spike times [seconds]
    int               *spikeindex;  // time index  of spikes
    count_type *dN;   // dN[t] spike count in [t,t+1)
    double *isi;      // not really isi, but "time since last spike"!
    
    int startspike;
    int endspike;
    int nvalidspikes; // number of valid spikes

    double deltaTime; // size of time interval [units?]
    //    bool    allocSpikeTimes;

    //    void alloc();
    void init_mex(const mxArray *in, bool allocSpikeTimes);
    //    void destroy();
    void destroy_mex();
private:
    double meanRate;  // mean firing rate (over all)

public:
    TData(const mxArray *in=0);
    ~TData();

    void setMexInput(const mxArray *in);
    mxArray* getMexOutput();

public:
    inline int getNTraj() const { return ntraj; }
    inline int getNValidSpikes() const { return nvalidspikes; }
    inline int getStartspike() const { return startspike; }
    inline int getEndspike() const { return endspike; }
    inline double getTimestep() const { return deltaTime; }
    inline int getNTimesteps() const { return ntimesteps; }
    inline double getMeanRate() const { return meanRate; }

    // returns index into posarray, phasearray nearest to given time
    int timeToIndex(double time) const;
    double getMinTime() const { return timearray[0];};
    double getMaxTime() const { return timearray[ntimesteps-1];};
    double getTime(int tindex) const { return timearray[tindex];};

    /* Used to get data at a particular time */
    traj_type getTrajectory(int time) const; // get the trajectory on which the rat was running 
    double getPos(int time) const; // get the rat's linear distance
    double getXPos(int time) const; // get the rat's x position
    double getYPos(int time) const; // get the rat's y position
    double getPhase(int time) const; // get the rat's theta phase 
    double getIsi(int time) const; // get the time since last spike

    double getX(int time) const; // get x variable
    double getY(int time) const; // get y variable

    double* getXPtr() const; // get x variable
    double* getYPtr() const; // get y variable

    /* get information about spikes */
    double getSpikeTime(int spikeIndex) const; // time at which spike occurred
    int getSpikeTimeIndex(int spikeIndex) const; // time index at which spike occurred
    double getSpikePos(int spikeIndex) const; // pos at which spike occurred
    double getSpikeXPos(int spikeIndex) const; // pos at which spike occurred
    double getSpikeYPos(int spikeIndex) const; // pos at which spike occurred
    double getSpikePhase(int spikeIndex) const; // time at which spike occurred
    traj_type getSpikeTraj(int spikeIndex) const; // trajectory on which spike occurred

    double getSpikeX(int spikeIndex) const; // x at which spike occurred
    double getSpikeY(int spikeIndex) const; // y at which spike occurred

    /* get min, max over all or over a single trajectory */
    double getMinPos() const;
    double getMinPos(traj_type t) const;
    double getMaxPos() const;
    double getMaxPos(traj_type t) const;
    double getMinPhase() const;
    double getMinPhase(traj_type t) const;
    double getMaxPhase() const;
    double getMaxPhase(traj_type t) const;
		
private:
    void init();
};



/* ************************************************************
   class AdaptCellData
   @@not used yet, introduce in future
   ************************************************************ */

class AdaptCellData {
public:
    double      *spiketimes;  // the times of the spikes in seconds
    int          nspikes;     // the total number of spikes
    TData   *behav;       // behavior data

    void alloc() {};
    void init_mex(const mxArray *in, bool allocSpikeTimes) {};
    void destroy() {};
};


}; // namespace 


#endif   // TDATA_H
