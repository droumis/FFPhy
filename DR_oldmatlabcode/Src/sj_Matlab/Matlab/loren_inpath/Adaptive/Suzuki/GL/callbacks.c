#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>
#include "adapt.h"

#define NTHETAP 10

/* global variables */
float epsilon[3];
float mle[3];
float stdev[3];
int 	presentation;
int 	spedup;
int    rtxi;
int    rtti;
int rti;
int posi;
int    rtxinc;
int    rttinc;
int    timeinc;
int animate, grid;
int currentspeed;
float *fixationtime;
float *timeinfo;
float starttime;
float *esttheta;
float *estthetaptr;
float *cpt;
float *cpisi;
float *thetacurrx;
float *thetacurrt;
float *realcpt;
float *realcpisi;
float *realthetat;
float *realthetaisi;
float *realthetacurrx;
float *realthetacurrt;
int nrealthetat;
int nrealthetaisi;
float *spiketime;
float *spiketimeptr;
float *currentspike;
float *spiketimeloc;
float *spikeisi;
float timestep;
float lastspiketime = 0;
float spiket = 0;

float tscale;

float lambdat[NLAMBDAXP];
float lambdaisix[NLAMBDATP];
float lambdatorig[NLAMBDAXP];
float lambdaisiorig[NLAMBDATP];

int 	dumpppm = 0;
int     ppmfilenum = 10000;
unsigned char    *pixels;
int     width = INIT_WINDOW_WIDTH;
int     height = INIT_WINDOW_HEIGHT;

int ntrials;
int currenttrial;
int *trialinfo;
int *trialID;
int *ntimesteps;
int *correct;
int totaltimesteps;
int maxtimesteps;
int currenttrialtime;
int currenttime;
int nspikes;
int lastpass = 0;
int nspikesused = 0; 
int nspikestmp = 0;
int ncpt;
int ncpisi;
int realncpt;
int realncpisi;
int    interpspline = 0;
int nrthetaparms;
float    leftwall;
float    rightwall;
float    traj2start;
float    peakratex;
float    peakratet;

typedef struct _WindowLocations {
    float bottom;
    float tbottom;
    float ttop;
    float xbottom;
    float xtop;
    float top;
    float left;
    float right;
} WindowLocations; 

typedef struct _TextInfo {
    float timeloc[2];
    float xtitleloc[2];
    char  xtitle[80];
    float xlabelloc[2];
    char  xlabel[80];
    int   nxxtick;
    float **xxtickloc;
    char  **xxtick;
    int   nxytick;
    float **xytickloc;
    char  **xytick;
    float ttitleloc[2];
    char ttitle[80];
    float tlabelloc[2];
    char  tlabel[80];
    int   ntxtick;
    float **txtickloc;
    char  **txtick;
    int   ntytick;
    float **tytickloc;
    char  **tytick;
} TextInfo; 

WindowLocations winloc;
TextInfo    textinfo;

/* get points from files */

void Init(void) {
    
    int i, j;
    int nelem;    


    fprintf(stderr, "\nScanning input files\n");

    /* read in the times of the spikes */
    nspikes = getnelements(info.spikefile, sizeof(float));
    spiketime = (float *)malloc(sizeof(float)*nspikes);
    readdata(info.spikefile, spiketime, nspikes, sizeof(float)); 

    /* read in the trial information */
    nelem = getnelements(info.spikefile, sizeof(int));
    trialinfo = (int *) calloc(nelem, sizeof(int));
    readdata(info.trialinfofile, trialinfo, nelem, sizeof(int)); 
    ntrials = *trialinfo;
    trialID = trialinfo+1;
    ntimesteps = trialinfo + 1 + ntrials;
    correct = trialinfo + 1 + 2 * ntrials;

    /* read in the time information */
    nelem = getnelements(info.spikefile, sizeof(float));
    timeinfo = (float *)malloc(sizeof(float)*nelem);
    readdata(info.timefile, timeinfo, nelem, sizeof(float));
    /* the timestep is the first element of time */
    timestep = *timeinfo;
    /* fixationtime starts at element 1 */
    fixationtime = timeinfo+1;
    /* get the total number of timesteps */
    totaltimesteps = 0;
    for (i = 0; i < ntrials; i++) {
        totaltimesteps += round(ntimesteps[i]);
	if (round(ntimesteps[i]) > maxtimesteps) {
	    maxtimesteps = round(ntimesteps[i]);
        }
    }
    maxtimesteps = 3000;
    

    nelem = getnelements(info.estthetafile, sizeof(float));
    esttheta = (float *) malloc(sizeof(float) * nelem);
    readdata(info.estthetafile, esttheta, nelem, sizeof(float));
    /* the first two values of estthetafile are the number of x control points and the
     * number of t control points */
    ncpt = round(esttheta[0]);
    ncpisi = round(esttheta[1]);
    cpt = esttheta+2;
    cpisi = cpt + ncpt;
    /* the start of the second trajectory is the middle control point */
    traj2start = cpt[(int) ncpt/2]; 
    /* thetacurrx is the current array of x control point heights */
    thetacurrx = (float *) malloc(sizeof(float) * ncpt);
    /* thetacurrt is the current array of t control point heights */
    thetacurrt = (float *) malloc(sizeof(float) * ncpisi);
    /* copy the starting values into cpt and cpisi */
    for (i = 0; i < ncpt; i++) {
        thetacurrx[i] = esttheta[2+ncpt+ncpisi+i];
    }
    for (i = 0; i < ncpisi; i++) {
        thetacurrt[i] = esttheta[2+2*ncpt+ncpisi+i];
    }
    /* advance the pointer past the control point initial values */
    esttheta = esttheta + 2 * (1 + ncpt + ncpisi);    
    estthetaptr = esttheta;
    
    /* if this is to be compared to simulated data, read in the data 
    if (info.simdata) {
        /* start with the x control points 
        nelem = getnelements(info.realthetatfile, sizeof(float));
        realthetat = (float *) malloc(sizeof(float) * nelem);
        readdata(info.realthetatfile, realthetat, nelem, sizeof(float));
        /* the first value is the number of x control points 
        realncpt = round(realthetat[0]);
        realcpt = realthetat+1;
        /* increment realthetat to start at the beginning of the data. Note that this
         * leaves a memory hole of 1 element 
        realthetat = realthetat + (1 + realncpt);
	nrealthetat = round((nelem - (1 + realncpt)) / (realncpt));
        /* check the size of realthetat. If it contains only two columns of data, we need
         * to interpolate between those two columns to get the appropriate values at each
         * time step 

        /* do the t control points 
        nelem = getnelements(info.realthetaisifile, sizeof(float));
        realthetaisi = (float *) malloc(sizeof(float) * nelem);
        readdata(info.realthetaisifile, realthetaisi, nelem);
        /* the first value is the number of x control points 
        realncpisi = round(realthetaisi[0]);
        realcpisi = realthetaisi+1;
        /* increment realthetat to start at the beginning of the data. Note that this
         * leaves a memory hole of 1 element 
        realthetaisi = realthetaisi + (1 + realncpisi);
        nrealthetaisi = round((nelem - (1 + realncpisi)) / (realncpisi));

        /* We're going to interpolate between the elements of the input x and t splines 
	interpspline = 1;
	/* allocate space for the realthetacurrx and realthetacurrt arrays 
	realthetacurrx= (float *) malloc(sizeof(float) * realncpt);
	realthetacurrt = (float *) malloc(sizeof(float) * realncpisi);
	for (i = 0; i < realncpt; i++) {
	    realthetacurrx[i] = realthetat[i];
	}
	for (i = 0; i < realncpisi; i++) {
	    realthetacurrt[i] = realthetaisi[i];
	}
    } */

    spiketimeloc = (float *)malloc(sizeof(float)*nspikes);
    spikeisi = (float *)malloc(sizeof(float)*nspikes);

    printf(" . ");
    printf(" . ");
    printf(" . ");
    printf("Done , nspikes = %d, totaltimesteps = %d\n", nspikes,
	    totaltimesteps); fflush(stdout);

    spiketimeptr = spiketime; 

    /* allocate space for the pixels to be dumped */
    pixels = (unsigned char *) malloc(sizeof(unsigned char) * width * height * 3);


    leftwall = cpt[0];
    rightwall = maxtimesteps;
    peakratex = -1e100;
    peakratet = -1e100;

    estthetaptr = esttheta;
    for (i = 0; i < totaltimesteps; i++) {
        if (*(estthetaptr++) > 0) {
	    for (j = 1; j < 5; j++) {
		if ((*estthetaptr) > peakratex) {
		    peakratex = (*estthetaptr); 
		}
		estthetaptr++;
	    }
	    //if ((*estthetaptr >= 0) && (cpisi[round(*estthetaptr)] <= 1)) {
	    if (*estthetaptr > 0) {
		estthetaptr++;
		for (j = 6; j < 10; j++) {
		    if ((*estthetaptr) > peakratet) {
			peakratet = (*estthetaptr); 
		    }
		    estthetaptr++;
		} 
	    }
	    else {
		estthetaptr += 5;
	    }
	}
	else {
	    /* move on to the next time point */
	    estthetaptr += 9;
	}
    }
    /* reset the pointer */
    estthetaptr = esttheta;
    

    printf("Found peak values: %f, %f\n", peakratex, peakratet); fflush(stdout);

    /* set up the locations for the various sections of the display */
    winloc.bottom = -1 * BELOW_T_HEIGHT;
    winloc.tbottom = 0;
    winloc.ttop = winloc.tbottom + T_HEIGHT;
    winloc.xbottom = winloc.ttop + ABOVE_T_HEIGHT;
    winloc.xtop = winloc.xbottom + X_HEIGHT;
    winloc.top = winloc.xtop + ABOVE_X_HEIGHT;
    winloc.left = leftwall - (rightwall - leftwall) / 8;
    winloc.right = rightwall + (rightwall - leftwall) / 8;
    /* define grid */
    define_grid();
    /* define grid labels */
    definelabels();


    tscale = (rightwall - leftwall) / 3;

    /* define the sets of x positions for the spatial and temporal splines */
    definesplinepos();
    currentspike = spiketime;
    
    rtxi = 0;
    rtti = 0;
    currenttime = 0;
    currenttrialtime = 0;
    animate = 1;
    grid = 1;

    currentspeed = 10;
    rtxinc = realncpt * currentspeed;
    rttinc = realncpisi * currentspeed;
    timeinc = currentspeed;
    presentation = 0;
    spedup = 0;
    glEnable(GL_POINT_SMOOTH);


}


/* display function */
void Display(void) {
    
    int i;
    int cseg;

    glClear(GL_COLOR_BUFFER_BIT);
    
    draw_grid(); 
    drawgridlabels();  
    
    /*if (info.simdata) {
        if (interpspline) {
        calctheta(realncpt, realncpisi, realthetat, realthetaisi, nrealthetat,
                      nrealthetaisi, (float) posi / (float) ntimesteps, realthetacurrx,
                      realthetacurrt);
        }
        drawfieldt(realcpt, realthetacurrx, realncpt, 0.0, 1.0, 0, 1, 1, 0);  
        drawfieldisi(realcpisi, realthetacurrt, realncpisi, 0.0, 1.0, 0);  
    }*/
    drawfieldt(cpt, thetacurrx, ncpt, 1.0, 1.0, .3, 1.0, .5 , 0.2 ); 
    drawfieldisi(cpisi, thetacurrt, ncpisi, .3, 1.0, 1.0); 
    drawtimeloc();
    /* do drawspikes for each timestep */
    //drawspikes();
    drawtime();  

    glutSwapBuffers();

    if (dumpppm) {
        glReadPixels(1, 1, width, height, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid *) pixels);
        /* write the ppm file */
        writeppm(&ppmfilenum, pixels);
    }
    if ((currenttrial < ntrials)&&(animate == 1)) {
        /* check to see if we are in a valid trial */
	if (trialID[currenttrial]) {
	    /* update the current thetahat */
	    for (i = 0; i < timeinc; i++) {
	        /* break out of the loop if we have reached the end of this trial*/
		if (currenttrialtime == ntimesteps[currenttrial]) {
		    currenttrial++;
		    currenttrialtime = 0;
		    nspikestmp = 0;
		    break;
                }
		if ((cseg = *(estthetaptr)) > 0) {
		    thetacurrx[cseg-1] = *(++estthetaptr);
		    thetacurrx[cseg] = *(++estthetaptr);
		    thetacurrx[cseg+1] = *(++estthetaptr);
		    thetacurrx[cseg+2] = *(++estthetaptr);
		    estthetaptr++;
		}
		else {
		    estthetaptr += 5;
		}
		if ((cseg = *(estthetaptr)) > 0) {
		    thetacurrt[cseg-1] = *(++estthetaptr);
		    thetacurrt[cseg] = *(++estthetaptr);
		    thetacurrt[cseg+1] = *(++estthetaptr);
		    thetacurrt[cseg+2] = *(++estthetaptr);
		    estthetaptr++;
		}
		else {
		    estthetaptr += 5;
		}
		currenttrialtime++;
		currenttime++;
	    }
	    //if (!spedup & presentation && (posi > 6000)) {
		/* speed up a bunch */
	//	spedup = 1;
	//	currentspeed *= 64;
	//	rtxinc = realncpt * currentspeed;
	//	rttinc = realncpisi * currentspeed;
	//	timeinc = currentspeed;
	 //   }
	}
	else if (currenttrial < ntrials) {
	    /* skip over this trial by jumping by ntimesteps[currenttrial] */
	    estthetaptr += ntimesteps[currenttrial] * 10;
	    currenttime += ntimesteps[currenttrial];
	    currenttrial++;
	    currenttrialtime = 0;
        }
    }
    else {
        usleep(100000);
    }
}

/* change camera when window is resized */
void Reshape(int w, int h) {
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(winloc.left, winloc.right, winloc.bottom, winloc.top, -1000, 1000);
    glViewport(0, 0, w, h);
    width = w;
    height = h;
    glMatrixMode(GL_MODELVIEW);
    
}

void definesplinepos()
{
    float xstep;
    int i;

    xstep = (rightwall - leftwall)/NLAMBDAXP;


    for(i = 0; i < NLAMBDAXP; i++) {
        lambdatorig[i] = (i+1) * xstep;
        lambdat[i] = leftwall + i * xstep;
    }
    for(i = 0; i < NLAMBDATP; i++) {
        lambdaisiorig[i] = (float) pow(10,(float) (i * 3.0 / (float) NLAMBDATP)); 
        lambdaisix[i] = leftwall + log10(lambdaisiorig[i]) * tscale;
    }
}


void drawfieldt(float *cpt, float *cpy, int ncpt, float r1, float g1, float b1, float r2, 
        float g2, float b2) 
{
    float lambday[NLAMBDAXP];
    float *lt;
    float *ly;
    int i;

    
    glColor3d(r1, g1, b1);
    glLineWidth(3.0f);
    
    cardinal(cpt, cpy, ncpt, lambdatorig, lambday, NLAMBDAXP);
    lt = lambdat;
    ly = lambday;
    
    glBegin(GL_LINE_STRIP);
    for(i = 0; i < NLAMBDAXP; i++) {
        glVertex3f(*(lt++), winloc.xbottom + *(ly++),0.0f);
    }
    glEnd();
    
    glLineWidth(1.0f);
    
}

void drawfieldisi(float *cpisi, float *cpy, int ncpt, float r, float g, float b) {

    float lambday[NLAMBDATP];
    float tmp2;
    float *lt, *ly;
    int i;
    
    /* plot the temporal spline out to 1 second */
    glColor3d(r, g, b);
    glLineWidth(3.0f);
    
    cardinal(cpisi, cpy, ncpisi, lambdaisiorig, lambday, NLAMBDATP); 

    lt = lambdaisix;
    ly = lambday;
    
    tmp2 = peakratex / peakratet;
    /* scale the lambda to be on a log scale */
    glBegin(GL_LINE_STRIP);
    for(i = 1; i < NLAMBDATP; i++) {
        glVertex3f(*(lt++), winloc.tbottom + *(ly++)*tmp2, 0.0f);
    }
    glEnd();
}

void drawgridlabels(void) 
{
    int i, j, len;


    /*glNewList(GRIDLABELS, GL_COMPILE); */
    glColor3f(1.0,1.0,1.0);

    glRasterPos2f(textinfo.xtitleloc[0], textinfo.xtitleloc[1]);
    len = strlen(textinfo.xtitle);
    for (j = 0; j < len; j++) {
	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, textinfo.xtitle[j]);
    } 
    glPopMatrix();

    glRasterPos2f(textinfo.xlabelloc[0], textinfo.xlabelloc[1]);
    len = strlen(textinfo.xlabel);
    for (j = 0; j < len; j++) {
	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, textinfo.xlabel[j]);
    }
    for (i = 0; i < textinfo.nxxtick; i++) {
	glRasterPos2f(textinfo.xxtickloc[i][0], textinfo.xxtickloc[i][1]);
	len = strlen(textinfo.xxtick[i]);
	for (j = 0; j < len; j++) {
	    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, textinfo.xxtick[i][j]);
	}
    } 
    for (i = 0; i < textinfo.nxytick; i++) {
	glRasterPos2f(textinfo.xytickloc[i][0], textinfo.xytickloc[i][1]);
	len = strlen(textinfo.xytick[i]);
	for (j = 0; j < len; j++) {
	    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, textinfo.xytick[i][j]);
	}
    } 
    
    glRasterPos2f(textinfo.ttitleloc[0], textinfo.ttitleloc[1]);
    len = strlen(textinfo.ttitle);
    for (j = 0; j < len; j++) {
	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, textinfo.ttitle[j]);
    }
    glRasterPos2f(textinfo.tlabelloc[0], textinfo.tlabelloc[1]);
    len = strlen(textinfo.tlabel);
    for (j = 0; j < len; j++) {
	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, textinfo.tlabel[j]);
    }
    for (i = 0; i < textinfo.ntxtick; i++) {
	glRasterPos2f(textinfo.txtickloc[i][0], textinfo.txtickloc[i][1]);
	len = strlen(textinfo.txtick[i]);
	for (j = 0; j < len; j++) {
	    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, textinfo.txtick[i][j]);
	}
    } 
    for (i = 0; i < textinfo.ntytick; i++) {
	glRasterPos2f(textinfo.tytickloc[i][0], textinfo.tytickloc[i][1]);
	len = strlen(textinfo.tytick[i]);
	for (j = 0; j < len; j++) {
	    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, textinfo.tytick[i][j]);
	}
    } 
    glEndList();
}

void drawtime(void) {

    char tmpstring[50];
    float ave = 0;
    int len, i;

//    if ((currenttrial > 2) && (currenttrial < ntrials - 2)) {
//	for (i = currenttrial - 2; i < currenttrial + 3; i++) {
//	    ave += correct[i];
 //       }
//	sprintf(tmpstring, "%.1f %% correct", ave/5 %* 100);
 //   }
  //  else {
//	sprintf(tmpstring, "? %% correct");
 //   }
    if (trialID[currenttrial]) {
	if (correct[currenttrial]) {
	    glColor3f(0.0,1.0,0.0);
	}
	else {
	    glColor3f(1.0,0.0,0.0);

	}
	    sprintf(tmpstring, "trial %d", currenttrial+1);

	glRasterPos2f(textinfo.timeloc[0], textinfo.timeloc[1]);
	len = (int) strlen(tmpstring);
	for (i = 0; i < len; i++) {
	    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, tmpstring[i]);
	}
    }

}

void drawtimeloc(void) {

    glEnable(GL_POINT_SMOOTH);
    glPointSize(75.0f);
    glColor3d(1.0,1.0,1.0);

    glBegin(GL_POINTS);
    glVertex3f(currenttrialtime, winloc.xbottom, 0.0f);
    glEnd();

}

void drawspikes(void) 
{
    int i;
    float tmp2, tmp3;
    int timetmp;

    glColor3d(1.0,1.0,1.0);

    //nspikestmp = 0;

    timetmp = currenttrialtime ;
    for (i = 0; i < timeinc; i++) {
        while ((nspikesused < nspikes) && (*spiketimeptr >= timetmp) && 
		(*spiketimeptr < timetmp+1)) { 
            /* add the spike to the spikelist */
            spikeisi[nspikestmp] = timetmp - lastspiketime;
            lastspiketime = timetmp;
            spiketimeloc[nspikestmp++] = timetmp;
	    fprintf(stderr, "Spike at %d\n", timetmp);
	    spiketimeptr++;
	    nspikesused++;
        }
	timetmp++;
    }
    tmp3 =  peakratex * -0.1;
    tmp2 = winloc.xbottom + tmp3;
    glLineWidth(0.075f);
    glBegin(GL_LINES);
    for(i = 0; i < nspikestmp; i++) {
        /* draw the position spikes */
        glVertex3f(spiketimeloc[i], winloc.xbottom, 0.0f);
        glVertex3f(spiketimeloc[i], tmp2, 0.0f);
        /* draw the temporal spike */
        glVertex3f(leftwall+log10(spikeisi[i]*1e3)*tscale, 0, 0.0f);
        glVertex3f(leftwall+log10(spikeisi[i]*1e3)*tscale, tmp3, 0.0f);
    }

    glEnd();
}

void define_grid() 
{
    int i;
    float tmp, loctmp;

    tmp = winloc.ttop / peakratet; 

    glNewList(GRID, GL_COMPILE);
    glColor3f(0,0,1.0);
    glBegin(GL_LINES);
    /* put in one horizontal line at 0 and one at 1 for the temporal spline */
    glVertex3d(leftwall, 0, 0);
    glVertex3d(rightwall, 0, 0);
    /* line at 1 */
    loctmp = tmp;
    while (loctmp < winloc.ttop) {
	glVertex3d(leftwall, loctmp, 0);
	glVertex3d(rightwall, loctmp, 0);
	loctmp+=tmp;
    }
    /* line at max of t */
    //glVertex3d(leftwall, winloc.ttop, 0);
    //glVertex3d(rightwall, winloc.ttop, 0);

    /* lines every 10 Hz for X spline */
    glVertex3d(leftwall, winloc.xbottom, 0);
    glVertex3d(rightwall, winloc.xbottom, 0); 
    i = 0;
    while ((tmp = winloc.xbottom + 10 * i++) <= winloc.xtop) {
        glVertex3d(leftwall, tmp, 0);
        glVertex3d(rightwall, tmp, 0);
    } 
    /* do the vertical lines */
    /* the bottom plot is on a log scale so put grid lines at zero, .001, .01, 
     * .1, 1, and 10. Scale the leftwall to rightwall distances to be from 0.001 to 10
     * seconds on the log scale or 0 to 4 on a linear scale. To do so, multiply each time
     * by 1000, take the log10, and scale from leftwall to right wall */
    tscale = (rightwall - leftwall) / 3; /* 3 segments: .001->.01, .01->.1, .1->1  */
    glVertex3d(leftwall + log10(.001*1e3)*tscale,  0, 0);
    glVertex3d(leftwall + log10(.001*1e3)*tscale, peakratex , 0);
    glVertex3d(leftwall + log10(.01*1e3)*tscale,  0, 0);
    glVertex3d(leftwall + log10(.01*1e3)*tscale, peakratex , 0);
    glVertex3d(leftwall + log10(.1*1e3)*tscale,  0, 0);
    glVertex3d(leftwall + log10(.1*1e3)*tscale, peakratex , 0);
    glVertex3d(leftwall + log10(1.0*1e3)*tscale,  0, 0);
    glVertex3d(leftwall + log10(1.0*1e3)*tscale, peakratex , 0); 
    glVertex3d(leftwall + log10(10*1e3)*tscale,  0, 0);
    glVertex3d(leftwall + log10(10*1e3)*tscale, peakratex , 0); 

    /* the vertical lines on the top are at the boundries of different periods */
    glVertex3d(300, winloc.xbottom, 0);
    glVertex3d(300, winloc.xtop, 0);
    /* delay */
    glVertex3d(800, winloc.xbottom, 0);
    glVertex3d(800, winloc.xtop, 0);
    /* response */
    glVertex3d(1500, winloc.xbottom, 0);
    glVertex3d(1500, winloc.xtop, 0);
        
    glEnd();
    glEndList();
}

void draw_grid(void) {

    if(grid)
        glCallList(GRID);
    else {
        glColor3d(0,1.0,0);
        glBegin(GL_LINES);
        glVertex3d(-20, 0, 0);
        glVertex3d(370, 0, 0);
        glEnd();
    }
}

void definelabels() 
{
    int    times[] = {1.0, 10, 100, 1000}; 
    float    middle;
    int     i, len; 
    int     charwidth, stringlen;

    charwidth = glutBitmapWidth(GLUT_BITMAP_HELVETICA_18, 'S'); 

    middle = (leftwall + rightwall) / 2;

    textinfo.timeloc[0] = leftwall - (rightwall - leftwall) / 10;
    textinfo.timeloc[1] = winloc.xtop + ABOVE_X_HEIGHT / 2;

    
    sprintf(textinfo.xtitle, "Trial Time");
    stringlen = glutBitmapLength(GLUT_BITMAP_HELVETICA_18, (const char *)textinfo.xtitle);
    len = (int) strlen(textinfo.xtitle) + 2;
    textinfo.xtitleloc[0] = middle-len*2;
    textinfo.xtitleloc[1] = winloc.xtop + ABOVE_X_HEIGHT / 1.3;

    sprintf(textinfo.xlabel, "(ms)");
    len = (int) strlen(textinfo.xlabel);
    textinfo.xlabelloc[0] = middle-len*2;
    textinfo.xlabelloc[1] = winloc.xtop + ABOVE_X_HEIGHT / 3;

    /* do the xxticks */
    textinfo.nxxtick = 3;
    textinfo.xxtickloc = (float **) malloc(sizeof(int *) * textinfo.nxxtick);
    textinfo.xxtick = (char **) malloc(sizeof(char *) * textinfo.nxxtick);
    for (i = 0; i < textinfo.nxxtick; i++) {
	textinfo.xxtickloc[i] = (float *) malloc(sizeof(int) * 2);
	textinfo.xxtickloc[i][1] = winloc.xtop + ABOVE_X_HEIGHT / 6;
	textinfo.xxtick[i] = (char *) malloc(sizeof(char) * 10);
    }
    textinfo.xxtickloc[0][0] = leftwall + 300 - 30; 
    sprintf(textinfo.xxtick[0], "%3d", 300);
    textinfo.xxtickloc[1][0] = leftwall + 800 - 30; 
    sprintf(textinfo.xxtick[1], "%3d", 800);
    textinfo.xxtickloc[2][0] = leftwall + 1500 - 30; 
    sprintf(textinfo.xxtick[2], "%3d", 1500);

    /* do the xyticks */
    textinfo.nxytick = (int) (peakratex / 10) + 1;
    textinfo.xytickloc = (float **) malloc(sizeof(int *) * textinfo.nxytick);
    textinfo.xytick = (char **) malloc(sizeof(char *) * textinfo.nxytick);
    for (i = 0; i < textinfo.nxytick; i++) {
	textinfo.xytickloc[i] = (float *) malloc(sizeof(int) * 2);
	textinfo.xytickloc[i][0] = leftwall - (rightwall - leftwall) / 25;
	textinfo.xytickloc[i][1] = winloc.xbottom + i * 10 - 1;
	textinfo.xytick[i] = (char *) malloc(sizeof(char) * 10);
	sprintf(textinfo.xytick[i], "%d", i * 10);
    }

    sprintf(textinfo.ttitle, "Temporal Component");
    len = (int) strlen(textinfo.ttitle) + 2;
    textinfo.ttitleloc[0] = middle - len * 2 ;
    textinfo.ttitleloc[1] = winloc.ttop;

    sprintf(textinfo.tlabel, "ISI (ms)");
    len = (int) strlen(textinfo.tlabel);
    textinfo.tlabelloc[0] = middle - len * 2;
    textinfo.tlabelloc[1] = winloc.tbottom - BELOW_T_HEIGHT / 1.25;

    /* there are four tick labels for the xaxis of the temporal spline  */
    textinfo.ntxtick = 4;
    textinfo.txtickloc = (float **) malloc(sizeof(int *) * textinfo.ntxtick);
    textinfo.txtick = (char **) malloc(sizeof(char *) * textinfo.ntxtick);
    for (i = 0; i < textinfo.ntxtick; i++) {
	textinfo.txtickloc[i] = (float *) malloc(sizeof(int) * 2);
	textinfo.txtickloc[i][0] = leftwall+log10(times[i])*tscale - 5;
	textinfo.txtickloc[i][1] = winloc.tbottom - BELOW_T_HEIGHT / 3;
	textinfo.txtick[i] = (char *) malloc(sizeof(char) * 10);
	sprintf(textinfo.txtick[i], "%d", times[i]);
    }
    /* move the first one back to the right */
    textinfo.txtickloc[0][0] += 5;


    textinfo.ntytick = round(peakratet);
    textinfo.tytickloc = (float **) malloc(sizeof(int *) * textinfo.ntytick);
    textinfo.tytick = (char **) malloc(sizeof(char *) * textinfo.ntytick);
    for (i = 0; i < textinfo.ntytick; i++) {
	textinfo.tytickloc[i] = (float *) malloc(sizeof(int) * 2);
	textinfo.tytickloc[i][0] = leftwall - (rightwall - leftwall) / 35;
	textinfo.tytickloc[i][1] = winloc.tbottom + (i * winloc.ttop/peakratet) - 1;
	textinfo.tytick[i] = (char *) malloc(sizeof(char) * 10);
	sprintf(textinfo.tytick[i], "%d", i);
    }

}

/* menu callback */
void main_menu(int item) {

    switch(item) {

    case 'r':
        lastspiketime = 0;
	spiketimeptr = spiketime; 
	nspikesused = 0;
	currenttime = 0;
	currenttrial = 0;
        rtxi = 0;
        rtti = 0;
	lastpass = 0;
        estthetaptr = esttheta;
        posi = 0;
        dumpppm = 0;
        presentation = 0;
        break;
    case 'R':
        presentation = 1;
	currentspeed = 24;
	rtxinc = realncpt * currentspeed;
	rttinc = realncpisi * currentspeed;
	timeinc = currentspeed;
        lastspiketime = 0;
	spiketimeptr = spiketime; 
	nspikesused = 0;
        rtxi = 0;
        rtti = 0;
	lastpass = 0;
        estthetaptr = esttheta;
        posi = 0;
        //dumpppm = 0;
        break;
    case 'p':
        //lastspiketime = 0;
	//spiketimeptr = spiketime; 
	//nspikesused = 0;
	//timeptr = time;
        //rtxi = 0;
        //rtti = 0;
	//lastpass = 0;
        //estthetaptr = esttheta;
        //posi = 0;
        dumpppm = 1;
	ppmfilenum = 10000;
        break;
    case 'b':
        /* backup 10 sec */
        break;
        
    case 's':
        animate ^= 1;
        break;
    case 'm':
        grid ^= 1;
        break;
    case 'u':
        currentspeed *= 2;
        rtxinc = realncpt * currentspeed;
        rttinc = realncpisi * currentspeed;
        timeinc = currentspeed;
        break;
    case 'd':
        currentspeed /= 2;
        if (currentspeed < 1) {
            currentspeed = 1;
        }
        rtxinc = realncpt * currentspeed;
        rttinc = realncpisi * currentspeed;
        timeinc = currentspeed;
        break;
    default:
        break;
    }
}

/* keyboard callback */
void Keyboard(unsigned char key, int x, int y) {


}

void calctheta(int ncpt, int ncpisi, float *thetat, float *thetaisi, int nthetat,
               int nthetaisi, float prop, float *thetacurrx, float *thetacurrt)
    /* prop is equal to the current timestep index / ntimestamps */
{

    int i;
    int tstart, isistart;
    int ncp;
    float tmp;

    /* find the pair of thetat and thetaisi variables to interpolate between */
    tstart = (int) floor((nthetat-1) * prop);
    isistart = (int) floor((nthetaisi-1) * prop);

    if (tstart == ((nthetat-1) * prop)) {
        /* this is right at tstart, so just set thetacurrx to the tstart element */
        ncp = ncpt;
        for (i = 0; i < ncp; i++) {
            thetacurrx[i] = thetat[tstart*ncp+i];
        }
    }
    else {
        /* figure out the new proportion between tstart and tstart + 1 */
        tmp = prop * (float) (nthetat-1) - (float) tstart;

        /* interpolate between the two columns of thetat */
        ncp = ncpt;
        for (i = 0; i < ncp; i++) {
            thetacurrx[i] = thetat[tstart*ncp+i] * (1 - tmp) +
                            thetat[(tstart+1)*ncp+i] * tmp;
        }
    }

    if (isistart == ((nthetaisi-1) * prop)) {
        /* this is right at isistart, so just set thetacurrt to the isistart element */
        ncp = ncpisi;
        for (i = 0; i < ncp; i++) {
            thetacurrt[i] = thetaisi[isistart*ncp+i];
        }
    }
    else {
        /* figure out the new proportion between isistart and isistart + 1 */
        tmp = prop * (float) (nthetaisi-1) - (float) isistart;

        /* interpolate between the two columns of thetaisi */
        ncp = ncpisi;
        for (i = 0; i < ncp; i++) {
            thetacurrt[i] = thetaisi[isistart*ncp+i] * (1 - tmp) +
                            thetaisi[(isistart+1)*ncp+i] * tmp;
        }
    }

    return;
}

void writeppm(int *ppmfilenum, unsigned char *pixels)
{
        FILE *ppmfile;
        char tmpstring[100];
        int i;

        sprintf(tmpstring, "a%d.ppm", *ppmfilenum);
        if ((ppmfile = fopen((const char *) tmpstring, "w")) == NULL) {
                fprintf(stderr, "Error opening %s for writing\n", tmpstring);
                exit(1);
        }

        /* write out the header information */
        fprintf(ppmfile, "P6\n%d %d\n255\n", width, height);

        /* write out the bytes from upper left to lower right */
        for (i = height - 1; i >= 0; i--) {
                fwrite(pixels+(i*width*3), sizeof(char), width*3, ppmfile);
        }
        fclose(ppmfile);
        (*ppmfilenum)++;
        return;
}

