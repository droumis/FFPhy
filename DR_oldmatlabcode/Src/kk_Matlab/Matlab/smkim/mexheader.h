#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#define sqr(x) ((x) * (x))
#define absv(x) ((x) > 0) ? (x) : (-1.0 * (x))
#define round(x) (((x) - .5 >= (int) (x)) ? ((int) ceil((x))) : (int) (x))
#define TIMESCALE 1e4
#define TTPOINTS 40	/* 40 points per individual waveform */
#define TTRECSIZE 160  	/* 160 points per four channel waveform */
#define MAXWAVEFORMVAL	SHRT_MAX
#define MINWAVEFORMVAL	SHRT_MIN

 
typedef struct {
   unsigned long        time;
   short        x1, y1, x2, y2;
} PosRec;

typedef struct {
	unsigned long 		timestamp;
	int			numsamples;
	double			sampfreq;
        double                  samplingrate;
        short			*data;
} ContRec;

typedef struct {
	unsigned long	time;
   	short			x[TTRECSIZE]; 
} TTRec; 

typedef unsigned long *SpikeList ;

typedef struct _SpikePos {
    double   timestamp;
#ifdef EMERY
	double 		x[2];
#else
    short   x[2];
#endif
} SpikePos;

typedef struct _PlaceField {
    double   mu[2];
    double   sigmasq[2];
    double   roe;
	double	 rmax;
} PlaceField;

typedef struct {
   unsigned long time;
   char s1, s2;
} Flag;
 
typedef struct {
   float x, y;
} Point;
 
#define     SECTSCALE  10000.0
#define     MSECTSCALE  10.0
#define     NUMPARMS    5
#define     TOLERANCE   .1
#define     EPSILON     .01
#define		TRACKERRATE	20.0	/* assume 20 Hz sampling */
#define 	SCALEFACTOR		10000.0
#define 	MINSPIKES 50

#define TRACKER_XPIXELS     364
#define	TRACKER_YPIXELS 	256
#define	TRACKER_RATE		60 		/* 60 Hz */
