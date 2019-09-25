#include <stdio.h>

/* definitions */
#define round(x) (((x) - .5 >= (int) (x)) ? ((int) ceil((x))) : (int) (x))

/***************** CONSTANTS ******/

#define INIT_WINDOW_WIDTH    800
#define INIT_WINDOW_HEIGHT    600

#define PI                  3.14159265359f

#define DEBUG             0 
#define RAT                  1
#define GRID              2
#define GRIDLABELS          3
#define SPIKES              4

/* heights of various parts of the window */
#define BELOW_T_HEIGHT        	(0.2*peakratex)
#define T_HEIGHT        	(peakratex)
#define ABOVE_T_HEIGHT        	(0.2*peakratex)
#define X_HEIGHT        	(peakratex)
#define ABOVE_X_HEIGHT        	(0.3*peakratex)

/* numbers of points */
#define NLAMBDAXP        200
#define NLAMBDATP        200


/***************** STRUCTURES ******/

typedef struct _ViewInfo {
    int simdata;
    char *trialinfofile;
    char *timefile;
    char *spikefile;
    char *estthetafile;
    char *realthetaxfile;
    char *realthetatfile;
} ViewInfo;

ViewInfo info;

/**************** FUNCTIONS ******/

/* adaptive simulation */
void Init(void);
void Iterate(float time);
void Display(void);
void Reshape(int w, int h);
void Keyboard(unsigned char key, int x, int y);
void main_menu(int item);
void drawfieldt(float *cpx, float *cpy, int ncp, float r1, float g1, float b1, float r2, 
        float g2, float b2);
void drawfieldisi(float *cpt, float *cpy, int ncp, float r, float g, float b);
void definesplinepos(void);
void cardinal(float *cpx, float *cpy, int ncp, float *x, float *y, int npoints);
void draw_grid(void);
void define_grid(void);
void definelabels(void);
void drawgridlabels(void);
void drawspikes(void);
void drawtime(void);
void drawtimeloc(void);
void calctheta(int ncpx, int ncpt, float *thetax, float *thetat, int nthetax,
                int nthetat, float prop, float *thetacurrx, float *thetacurrt);


 
/* produce spikes */
void Init(void);
void Iterate(float time);
void Display(void);
void Reshape(int w, int h);
void x_label(float x, float y, char *text);
void y_label(float x, float y, char *text);

/* read input */
int readdata(char *filename, void *data, int nsamples, int size);
int getnelements(char *filename, int esize);

void writeppm(int *ppmfilenum, unsigned char *pixels);


