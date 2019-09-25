// $Id: MovieControl.cc,v 1.1 2008/08/24 19:47:09 chengs Exp $

#include "../model/IsiModel.h"
#include "ModelVis.h"
#include "GUIMain.h"

#include "MovieControl.h"


extern GUIMain *gui;
static const double PLAYRATE_DEFAULT= 500; // [frames/ sec]
//static const unsigned short mpegFrameRate= 12; // [frames/sec]
static const unsigned short mpegFrameRate= 25; // [frames/sec]

MovieControl::MovieControl(int x, int y, int w, int h, const char *l)
	: Display(x,y,w,h,l)
{
	moviePlaying = false;
	
	playBeginPTime.tv_sec = 0;
	playBeginPTime.tv_usec = 0;
	playBeginMTime = 0;
	mTime = 0;
	mTimeMin = 0;
	mTimeMax = 0; 
    playRateNormal= PLAYRATE_DEFAULT;
	playRate = playRateNormal;

    dumpFrames = false;
    pixels = 0;
    writeFrameFct= 0;
    outputFormat= 0;
    row_ptr= 0;
}

MovieControl::~MovieControl() { 
//    hippo_Print((double*) pixels);
    if(model) events->removeMovieEventListener(this); 
    if(pixels) delete[] pixels;
    if(row_ptr) delete[] row_ptr;
//    hippo_Mark;
}

void MovieControl::glInit()
{
	glClearColor (0.0f, 0.0f, 0.0f, 0.0f);	
	glClearDepth (1.0f);

	glDisable(GL_LIGHTING);

	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

//    gl_font(FL_COURIER,10);
	gl_font(FL_HELVETICA_BOLD,16);
	
	model->glInit();
}

// Main methods for graphing the adaptive filter results.
void MovieControl::glDraw()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0,1,0,6,-100,100);
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/* draw geometry plot */
	{
		glPushMatrix();
        glTranslatef(0.6f,1.7f,0.0f);
        glScalef(0.3f,-1.7f,1.0f);

        model->plotGeometry(mTime);

		glPopMatrix();
	}

	/* draw spatial functions */
    glPushMatrix();
    glTranslatef(0.0f, 4.8f, 0.0f);
    glScalef(1.0f,1.1f,1.0f);

//    hippo_Print(mTime);

    model->plotTraj(mTime);
    model->plotTrajActivity(mTime, getPlayDir() );
    glPopMatrix();

	/* draw isi plot */
    glPushMatrix();
    glTranslatef(0.03f,0.4f,0);
    glScalef(0.5f,1.0f,1.0f);
    model->plotIsi(mTime);
    model->plotIsiActivity(mTime, getPlayDir());
    glPopMatrix();
	
    // print real time of movie
    char temp[80];
    double time= data->getTime(mTime)-data->getMinTime();
//    glTranslatef(0.01f, 0.0f, 0.0f);
    glTranslatef(0.85f, 5.85f, 0.0f);
    snprintf(temp,80,"time: %5.0f sec",time);
    glRasterPos3f(0.0f, 0.0f, 0.0f);
    gl_draw(temp);
    
    if (dumpFrames) (this->*writeFrameFct)();
}

void MovieControl::setMTime(int t)
{
    mTime= t;
    setPlayBegin(t);
    updateFrame();
}

void MovieControl::setModel(ModelVis *d) { 
    model= d;
    events= model->getEventHandler();
    events->addMovieEventListener(this); 

    mTimeMin = model->getModel()->getStartindex();  
    mTimeMax = model->getModel()->getEndindex();  
    mTime= mTimeMin;
    data= model->getData();
    nTraj= data->getNTraj();
    
}

void MovieControl::play() {
    setPlayRate(playRateNormal);
    _play();
}

void MovieControl::pause() {
    if (moviePlaying) Fl::remove_idle(MovieControl::cb_movie_idle);
    moviePlaying = false;
}

void MovieControl::stop() {
    pause();
    setMTime(mTimeMin);
}

void MovieControl::ff() 
{
    double r= getPlayRate();
    if (r < 0) r= playRateNormal;
    else r*= 2;
    setPlayRate(r);
    _play();
}

void MovieControl::rew()
{
    double r= getPlayRate();
    if (r > 0) r= -playRateNormal;
    else r*= 2;
    setPlayRate(r);
    _play();
}

void MovieControl::_play()
{
    if (!moviePlaying) {
        // TODO: this probably needs to be changed!!
        Fl::add_idle(MovieControl::cb_movie_idle);
        moviePlaying = true;
    }
}

void MovieControl::updateFrame()
{
    if(!dumpFrames) { // regular movie playback
        if(moviePlaying) {
            timeval t;
            gettimeofday(&t, NULL);
            _dt = 1e6*(t.tv_sec - playBeginPTime.tv_sec) + (t.tv_usec - playBeginPTime.tv_usec);
            mTime = playBeginMTime + (int)(playRate * 1e-6*_dt);
        }
    } else { // frame dump mode, don't link playback to real time
        mTime+= (int) (playRate/ mpegFrameRate);
    }

    // handle out-of-bounds mTime
    if (mTime < mTimeMin || mTime > mTimeMax) {
        // Stop the movie before the time gets out of bounds.
        // (in terms of the number of timesteps that the adaptive model can eval)
        if (moviePlaying) stop(); // stop movie 
        if (mTime >= mTimeMax)
            mTime= mTimeMax;
        else
            mTime= mTimeMin;
    }

    events->notifyMTimeEventListeners();
    for(int n= 0; n < nTraj; n++)
        events->notifyTrajMeshEventListeners(n);
    events->notifyMovieEventListeners();
}

void MovieControl::setDumpFrames(bool b)
{
    dumpFrames= b;
    const int width= w();
    const int height= h();
    if(pixels) { delete[] pixels; pixels= 0; }
    if(row_ptr) { delete[] row_ptr; row_ptr= 0; }
    if(dumpFrames) {
        pixels = new unsigned char[3*width*height];

        switch(outputFormat) {
            case 2:
                hippo_Msg("Writing frames to ppm files");
                writeFrameFct= &MovieControl::writeppm;
            case 1:
                hippo_Msg("Writing frames to png files");
                writeFrameFct= &MovieControl::writepng;
                row_ptr= new png_bytep[height];
                for(int i= 0; i<height; i++) row_ptr[height-1-i]= pixels+ 3*width*i;
        }
    }
}

void MovieControl::writeppm() {
    //    hippo_Print((double*) pixels);
    GLint viewport[4];                  // Where The Viewport Values Will Be Stored
    glGetIntegerv(GL_VIEWPORT, viewport);           // Retrieves The Viewport Values (X, Y, Width, Height)
    glReadPixels(viewport[0],viewport[1],viewport[2],viewport[3], GL_RGB, GL_UNSIGNED_BYTE, (GLvoid *) pixels);
    //glReadPixels(1,1,viewport[2],viewport[3], GL_RGB, GL_UNSIGNED_BYTE, (GLvoid *) pixels);

    FILE *ppmfile;
    char fname[100];

    int width = w();
    int height = h();

    snprintf(fname, 100, "frameDump/m%.7d.ppm", mTime);
    if ((ppmfile = fopen((const char *) fname, "w")) == NULL) {
        fprintf(stderr, "Error opening %s for writing\n", fname);
        exit(1);
    }

    /* write out the header information */
    fprintf(ppmfile, "P6\n%d %d\n255\n", width, height);

    /* write out the bytes from upper left to lower right */
    for (int i = height - 1; i >= 0; i--) {
        fwrite(pixels+(i*width*3), sizeof(char), width*3, ppmfile);
    }
    fclose(ppmfile);
    return;
}

void MovieControl::writepng() {

    // open png file for write
    char fname[100];
    snprintf(fname, 100, "frameDump/m%.7d.png", mTime);
    FILE *fp = fopen(fname, "wb");
    if (!fp) {
        fprintf(stderr, "Error opening %s for writing\n", fname);
        exit(1);
    }

    png_structp png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, png_voidp_NULL,
            png_error_ptr_NULL, png_error_ptr_NULL);
    hippo_Assert(png_ptr, "could not allocate png_ptr");

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
        hippo_ErrorMsg("could not allocate info_ptr");
    }

    //    png_set_IHDR(png_ptr, info_ptr, width, height, bit_depth, color_type, 
    //            interlace_type, compression_type, filter_method)
    png_set_IHDR(png_ptr, info_ptr, w(), h(), 8,  
            PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, 
            PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_set_rows(png_ptr, info_ptr, row_ptr);
    png_init_io(png_ptr, fp);

    //    hippo_Print((double*) pixels);
    GLint viewport[4];                  // Where The Viewport Values Will Be Stored
    glGetIntegerv(GL_VIEWPORT, viewport);           // Retrieves The Viewport Values (X, Y, Width, Height)
    glReadPixels(viewport[0],viewport[1],viewport[2],viewport[3], GL_RGB, GL_UNSIGNED_BYTE, (GLvoid *) pixels);
    //glReadPixels(1,1,viewport[2],viewport[3], GL_RGB, GL_UNSIGNED_BYTE, (GLvoid *) pixels);

    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    png_write_end(png_ptr, info_ptr);

    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
}

void MovieControl::quit()
{
    while( Fl::first_window() ) {
        cerr << Fl::first_window() << "\n";
        Fl::first_window()->hide();
    }
}

void MovieControl::setMexInput(const mxArray *opts)
{
    MX_FieldScalarDefault(opts, "playrate", playRateNormal, PLAYRATE_DEFAULT, double);

    mxArray *mxtmp= mxGetField(opts, 0, "outputFormat");
    char *tmpstring;
    if(mxtmp) {
        MX_StringAllocAssign(mxtmp, tmpstring);
        hippo_Print(tmpstring);
        if(strcmp(tmpstring, "ppm"))
            outputFormat= 2;
        else if(strcmp(tmpstring, "png"))
            outputFormat= 1;
        else
            hippo_Error("Unknown outputFormat", tmpstring);
        mxFree(tmpstring);
    } else outputFormat= 1;

}

void MovieControl::cb_movie_idle(void *v) {
    gui->mcontrol->updateFrame();
}


