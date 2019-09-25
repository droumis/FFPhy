/* $Id: TData.cc,v 1.10 2008/08/24 19:46:51 chengs Exp $
   
Sen Cheng, 2004/../..

program description
*/

#include <stdlib.h>
#include <math.h>
#include "../aux/hippoIO.h" 
#include "../aux/mexAux.h" 
#include "TData.h"

#include <float.h> // for DBL_MAX

using namespace AFilter;

TData::TData(const mxArray *in)
    : ntimesteps(-1), fieldID(0), spiketimes(0), spikeindex(0), dN(0), isi(0)
{
    xpos= ypos= 0;
    x= y= 0;
    startspike= 0;
    endspike= 0;
    nvalidspikes= 0;
    meanRate= 0;

    if(in) setMexInput(in);
}

TData::~TData()
{
    if (isi) { free(isi); }
    if (dN) free(dN);
    if (fieldID) delete[] fieldID;
    if (spikeindex) delete[] spikeindex;
}

void TData::init() 
{

    deltaTime= timearray[1]-timearray[0]; //@@ what if timearray is not equidistant


    hippo_Assert(dN==0, "dN was previously assigned");
    hippo_Assert(isi==0, "isi was previously assigned");

    hippo_CAlloc(isi,ntimesteps, double);

    hippo_CAlloc(dN, ntimesteps, count_type);

    // cut off spiketimes before timearray
    startspike= 0;
    while(startspike < nspikes && spiketimes[startspike] < timearray[0]) 
        startspike++;
    // cut off spiketimes after timearray
    endspike= nspikes-1;
    double maxTime= timearray[ntimesteps-1];
    while (endspike > startspike  && spiketimes[endspike] > maxTime) 
        endspike--;
    nvalidspikes= endspike-startspike+1;

//    hippo_Print(startspike);
//    hippo_Print(endspike);
//    hippo_Print(nspikes);

    // find timestep index of valid spikes
    spikeindex= new int[nspikes];
    for(int j= startspike; j <= endspike; j++)
        spikeindex[j]= timeToIndex(spiketimes[j]);

    int t;
    for (t=0; t < ntimesteps; t++) {
        dN[t]= 0; 
        isi[t]= -1;
    }

    // which "times since last spike" (TLS) are valid?
    if(nvalidspikes > 0) {
        // TLS up to (and incl) first spike are invalid
        t= spikeindex[startspike];
        double  currenttime= timearray[t], lastspiketime= currenttime;
        int nextspike= startspike;
        if(t < ntimesteps) {
            while(nextspike <= endspike && 
                    spiketimes[nextspike] < currenttime+deltaTime) {
                dN[t]++;
                nextspike++;
            }
            t++;
        }
        for (; t < ntimesteps; t++) {
            currenttime= timearray[t];
            //            printf("%.3f\t", isi[t]);
            isi[t]= currenttime - lastspiketime;
            while(nextspike <= endspike && 
                    spiketimes[nextspike] < currenttime+deltaTime) {
                dN[t]++;
                lastspiketime= currenttime;
                nextspike++;
            }
        }
    } 

    // mean firing rate
    meanRate= nvalidspikes/(timearray[ntimesteps-1]-timearray[0]);
//    hippo_Print(this);
//    hippo_Print(meanRate);
//    hippo_Print(getMeanRate());
    
    cout << "Data initialized: ";
    cout << ntimesteps << " timesteps, " ;
    cout << nvalidspikes << " valid spikes" << endl;

    // [histogram isi's]
//    const int N= 180;
//    int h[N], ind, sum= 0;
//    for(int i=0; i<N; i++) h[i]=0;
//    for (t= 0; t < ntimesteps; t++) {
//        ind= (int) round(1000*isi[t]);
//        if(ind >= N|| ind < 0) ind= N-1;
//        h[ind]++;
//    }
//    for(int i=0; i<N; i++) {
//        cerr << i << '\t' << h[i] << '\t';
//        if(i<N-2) cerr << h[i]-h[i+2] ;
//        cerr << '\n';
//        sum+= h[i];
//    }
//    cerr << "sum= " << sum << "\t, ntimesteps= " << ntimesteps << endl;

    // [print time since last spike]
//    printf("%.3f\t%.3f\n%.3f\t%.3f\n", spiketimes[startspike],
//    timearray[spikeindex[startspike]],
//    spiketimes[startspike+1],
//    timearray[spikeindex[startspike+1]]);
//        for (t=spikeindex[startspike]; t < spikeindex[startspike+1]+10; t++) {
//            printf("k= %d, t=%.3f, dN[k]= %d, %.3f\n", t, timearray[t], 
//                    dN[t], isi[t]);
//        }

}


void TData::setMexInput(const mxArray *in)
{
    int i;
    MX_FieldAssign(in, "time", timearray,ntimesteps);
    //    hippo_Print(ntimesteps);

    mxArray *cellPtr;
    cellPtr= mxGetField(in,0,"linpos");
    if(cellPtr) {
        MX_FieldAssign(in, "linpos", posarray,i);
        hippo_Assert(ntimesteps==i,"posarray must have same length as timearray");
    } 

    cellPtr= mxGetField(in,0,"x");
    if(cellPtr) {
        MX_FieldAssign(in, "x", x,i);
        hippo_Assert(ntimesteps==i,"xmust have same length as timearray");
    } 

    cellPtr= mxGetField(in,0,"y");
    if(cellPtr) {
        MX_FieldAssign(in, "y", y,i);
        hippo_Assert(ntimesteps==i,"ymust have same length as timearray");
    } 

    cellPtr= mxGetField(in,0,"xpos");
    if(cellPtr) {
        MX_FieldAssign(in, "xpos", xpos,i);
        hippo_Assert(ntimesteps==i,"xpos must have same length as timearray");
    } 

    cellPtr= mxGetField(in,0,"ypos");
    if(cellPtr) {
        MX_FieldAssign(in, "ypos", ypos,i);
        hippo_Assert(ntimesteps==i,"ypos must have same length as timearray");
    } 

    cellPtr= mxGetField(in,0,"phase");
    if(cellPtr) {
        MX_FieldAssign(in, "phase", phasearray,i);
        hippo_Assert(ntimesteps==i,"phase must have same length as timearray");
    } 

    // fieldId's,  ntraj
    // thetahat is float -> need to convert
    double *tmp;
    short int trajMax=0, trajTmp;
    MX_FieldAssign(in, "traj", tmp, i);
//        hippo_Print(ntimesteps);
    hippo_Assert(ntimesteps==i,"fieldID must have same length as timearray");
    fieldID= new traj_type[ntimesteps];
    for (i = 0; i < ntimesteps; i++) {
        trajTmp= (traj_type) *(tmp++);
        fieldID[i]= trajTmp;
        if (trajTmp > trajMax) trajMax= trajTmp;
    }
    ntraj= trajMax+1;

    cellPtr= mxGetField(in,0,"spiketimes");
    if(cellPtr) {
        MX_FieldAssign(in, "spiketimes", spiketimes,nspikes);
    } else {
        spiketimes= 0;
        nspikes= 0;
//        MX_Warn("Data object loaded with %d spikes", nspikes);
    }

    /*
    cellPtr= mxGetField(in,0,"spikeindex");
    if(cellPtr) {
        MX_FieldAssign(in, "spikeindex", tmp, i);
    } else {
        i= 0;
        spikeindex= 0;
    }
    hippo_Print(i);
    hippo_Print(nspikes);
    hippo_Assert(i== nspikes,"spikeindex must have same length as spiketimes");
    spikeindex= new int[nspikes];
    for (i = 0; i < nspikes; i++) spikeindex[i]= (int) *(tmp++);
    */

    hippo_Assert(ntimesteps >=2, "not enough data given");

    init();
}


mxArray* TData::getMexOutput()
{
    MX_Error("function empty %s","");
    return 0;
}

int TData::timeToIndex(double time) const
{ 
    return (int) floor((time+tiny-timearray[0])/deltaTime); 
}

traj_type TData::getTrajectory(int time) const
{
    hippo_Assert(time>=0 && time<ntimesteps, "time out of bounds");
    return fieldID[time];
}

double TData::getPos(int time) const
{
    hippo_Assert(time>=0 && time<ntimesteps, "time out of bounds");
    return posarray[time];
}

double TData::getXPos(int time) const
{
    hippo_Assert(time>=0 && time<ntimesteps, "time out of bounds");
    return xpos? xpos[time] : 0;
}

double TData::getYPos(int time) const
{
    hippo_Assert(time>=0 && time<ntimesteps, "time out of bounds");
    return ypos ? ypos[time] : 0;
}

double TData::getX(int time) const
{
    hippo_Assert(time>=0 && time<ntimesteps, "time out of bounds");
    return x? x[time] : 0;
}

double TData::getY(int time) const
{
    hippo_Assert(time>=0 && time<ntimesteps, "time out of bounds");
    return y? y[time] : 0;
}

double* TData::getXPtr() const 
{
    hippo_Assert(x, "variable x was not assigned in input");
    return x; 
};

double* TData::getYPtr() const 
{ 
    hippo_Assert(y, "variable y was not assigned in input");
    return y; 
};

double TData::getPhase(int time) const
{
    hippo_Assert(time>=0 && time<ntimesteps, "time out of bounds");
    return phasearray[time];
}

double TData::getIsi(int time) const
{
    hippo_Assert(time>=0 && time<ntimesteps, "time out of bounds");
    return isi[time];
}

// using std iterators and utility classes with operator() would make this a lot neater
double TData::getMinPos() const
{
    double min=DBL_MAX;
    for (int i=0; i<ntimesteps; i++)
        min = posarray[i] < min ? posarray[i] : min;
    return min;
}

double TData::getMinPos(traj_type t) const
{
    double min=DBL_MAX;
    for (int i=0; i<ntimesteps; i++) {
        if (fieldID[i] == t)
            min = posarray[i] < min ? posarray[i] : min;
    }
    return min;
}

double TData::getMaxPos() const
{
    double max= -DBL_MAX;
    for (int i=0; i<ntimesteps; i++)
        max = posarray[i] > max ? posarray[i] : max;
    return max;
}

double TData::getMaxPos(traj_type t) const
{
    double max=-DBL_MAX;
    for (int i=0; i<ntimesteps; i++) {
        if (fieldID[i] == t) 
            max = posarray[i] > max ? posarray[i] : max;
    }
    return max;
}

double TData::getMinPhase() const
{
    double min=DBL_MAX;
    for (int i=0; i<ntimesteps; i++) {
        min = phasearray[i] < min ? phasearray[i] : min;
    }
    return min;
}

double TData::getMinPhase(traj_type t) const
{
    double min=DBL_MAX;
    for (int i=0; i<ntimesteps; i++) {
        if (fieldID[i] != t) 
            min = phasearray[i] < min ? phasearray[i] : min;
    }
    return min;
}

double TData::getMaxPhase() const
{
    double max=-DBL_MAX;
    for (int i=0; i<ntimesteps; i++)
        max = phasearray[i] > max ? phasearray[i] : max;
    return max;
}

double TData::getMaxPhase(traj_type t) const
{
    double max=-DBL_MAX;
    for (int i=0; i<ntimesteps; i++) {
        if (fieldID[i] != t) 
            max = phasearray[i] > max ? phasearray[i] : max;
    }
    return max;
}

double TData::getSpikeTime(int index) const
{
    hippo_Assert(index >= startspike && index  <= endspike, "spike index out of bounds");
    return spiketimes[index];
}

int TData::getSpikeTimeIndex(int index) const
{
//    hippo_Print(startspike);
//    hippo_Print(index);
//    hippo_Print(endspike);
    hippo_Assert(index >= startspike && index  <= endspike, "spike index out of bounds");
//    hippo_Print(spikeindex[index]);
    return spikeindex[index];
}

double TData::getSpikePos(int index) const
{
    hippo_Assert(index >= startspike && index  <= endspike, "spike index out of bounds");
    return posarray[spikeindex[index]];
}

double TData::getSpikeXPos(int index) const
{
    hippo_Assert(index >= startspike && index  <= endspike, "spike index out of bounds");
    return xpos ? xpos[spikeindex[index]] : 0;
}

double TData::getSpikeYPos(int index) const
{
    hippo_Assert(index >= startspike && index  <= endspike, "spike index out of bounds");
    return ypos ? ypos[spikeindex[index]] : 0;
}

double TData::getSpikeX(int index) const
{
    hippo_Assert(index >= startspike && index  <= endspike, "spike index out of bounds");
    return x? x[spikeindex[index]] : 0;
}

double TData::getSpikeY(int index) const
{
    hippo_Assert(index >= startspike && index  <= endspike, "spike index out of bounds");
    return y? y[spikeindex[index]] : 0;
}

double TData::getSpikePhase(int index) const
{
    hippo_Assert(index >= startspike && index  <= endspike, "spike index out of bounds");
    return phasearray[spikeindex[index]];  
}

traj_type TData::getSpikeTraj(int index) const
{
    hippo_Assert(index >= startspike && index  <= endspike, "spike index out of bounds");
    return fieldID[spikeindex[index]];  
}

