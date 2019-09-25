/* $Id: afTrajMovie.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
 *
 * draw 2-d surface plot ?
 */

#ifndef __AFTRAJTEXTURE_H__
#define __AFTRAJTEXTYRE_H__

// OpenGL includes
#include <FL/gl.h>
#include <GL/glu.h>
#include <assert.h>
#include <vector>

// vertex_tmesh, vertex3f, ManageMatrix etc. are defined in here
#include "glUtils.h"
#include "trajInfo.h"

using namespace std;

class AFTrajMovie {
	public:
		AFTrajMovie(int traj, TrajInfo *trajInfo); 
		~AFTrajMovie();
	
		void setColormap(int cLen, color3f *c);
		
		/* 
		*  The following require a current OpenGL context 
		*  */
		/* Get a texture ID, set the texture interpolation modes, etc. */
		void init();
		
		/* Set this as the texture for OpenGL to render. */
		void bind() const { glBindTexture(GL_TEXTURE_2D,texId); }

		/* Call bind and then update the texture */	
		void subImage() const;

		void compute(int time);
		
		/* draws the texture, axes, and label into the rectangle with corners
		*  (0.0, 0.0) and (1.0, 1.0)
		*/
		void draw();
		
	private:
		int traj;
		TrajInfo *trajInfo;

		color3f *colormap;
		int colormapLength;
		
		/* these should be powers of 2 */
		int width;
		int height;
	
		bool initialized;
		
		GLuint texId; // the openGL assigned texture identifier
		color3f *pixels; // the (row-major) image data

		bool active; // is this the active trajectory? 
		float ratPos;

		// # timesteps a spike should remain on screen
		static const int SPIKE_EXIST; 
        // spikes to display
		vector<vertex2f> spikes; 
		// determines spikes on this traj between t1 and t2
		void computeSpikes(int t1, int t2);
        // draw spikes' in position-phase grid
        void drawSpikes();
        // draw only spikes' position on the position axis
        void drawSpikesPos();
		
		// axis drawing
		void drawPosAxis();
		void drawPhaseAxis();
		
		// the coordinates of the texture rect's corners
		static const float TX_L;
		static const float TX_R;
		static const float TX_B;
		static const float TX_T;

		// the offset of the axes from the texture
		static const float AO_L;
		static const float AO_B;
};

#endif
