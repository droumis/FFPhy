/* $Id: Activity.h,v 1.2 2008/10/23 21:24:10 chengs Exp $
   authors     : Sen Cheng
   created     : 2005/02/20
 */

#ifndef AF_Activity_H
#define AF_Activity_H

#include <vector>
#include "glUtils.h"
#include "DataCache.h"

using namespace std;
using namespace AFilter;

/* ************************************************************
                          class Activity
   ************************************************************ */

namespace AFilter {

class Plot1d;

/** Visualize spiking activity and animals's position.
 * @author: Eric Foley, Sen Cheng
 * @date:   
 */
class Activity 
{
private:
    // spikes to display
    struct {
        vector<vector<float> > pos, phase; 
        vector<float> isi; 
    } spikes;

    DataCache *cache;
    const TData *data;
    double minPos, maxPos;
    double minLogIsi, maxLogIsi;
    bool logisi;

    int nTraj;
    float ratPos;
    traj_type currtraj;
    float currIsi;

    int fDim;
private:
    // determines spikes on this traj between t1 and t2
    void computeSpikes(int t1, int t2);
public:
    static const int SPIKE_EXIST;// # timesteps a spike should remain on screen


public:
    Activity(DataCache *, Plot1d *isiPlot, int dim);
    Activity(Activity &);
    ~Activity();

    void compute(int time, int dir);

    void init();
    void drawTraj() { drawSpikesTraj(); drawPos(); }

    /// draw only isi's on the isi axis
    void drawIsi() { /*drawCurrIsi();*/ drawSpikesIsi(); };

    /// draw only isi's on the isi axis
    void drawCurrIsi();

    /// draw the isi of spikes' 
    void drawSpikesIsi();

    /// draw only spikes' position on the position axis
    void drawSpikesTraj();

    /// draw white rectangle at rat's position on x axis
    void drawPos();
}; // class Activity


} // namespace AFilter

#endif   // AF_Activity_H
