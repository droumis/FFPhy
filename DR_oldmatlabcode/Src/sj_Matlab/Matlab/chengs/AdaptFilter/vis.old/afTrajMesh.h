/* $Id: afTrajMesh.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
 * 
 * Draw 3-d surface plot.
 */
#ifndef __AFTRAJMESH_H__
#define __AFTRAJMESH_H__

// OpenGL includes
#include <FL/gl.h>
#include <GL/glu.h>

// vertex_tmesh, vertex3f, ManageMatrix etc. are define in here
#include "glUtils.h"
#include "trajInfo.h"

using namespace std;

class AFTrajMesh {
	public:
		AFTrajMesh();
		AFTrajMesh(int traj, TrajInfo *ti); 
		~AFTrajMesh();

		// This must be called before the mesh is drawn 
		// or computed.
		void setTrajInfo(int traj, TrajInfo *ti);
		void setColormap(int colormapLength, color3f *colormap);
		
		/*
		*  The following require a current OpenGL context 
		*  */
		void init();
		void draw();
		
		void compute(int time) { computeMesh(time); }
		
	private:
		void computeMesh(int time);
		void computeSpikes(int t1, int t2); // determines spikes on this traj btw t1 and t2
		
	private:
		void drawMesh();
		// void drawGrid();
		void drawAxes();
		void drawPosAxisGrid();
		void drawPhaseAxisGrid();
		void drawSpikesOnMesh();
		void drawSpikesProjected();
		
	private:
		void createVertexArray();
		void createMeshIndexArray();
		void createGridIndexArray();
		void createSpikeArray();		

	private:

		static const float AO; // axes offset from mesh
		static const float LO; // label offset from mesh
		static const float VE; // extent of value (y) axis

		int traj;
		TrajInfo *trajInfo;

		color3f *colormap;
		int colormapLength;
		
		bool initialized;
		int nVertsX, nVertsZ;
		vertex_tmesh *vertices;
		unsigned int *indices_mesh;
		unsigned int *indices_grid;

#define MAX_SPIKES 256
		vertex2f spikes[MAX_SPIKES];
		int nUsedSpikes; // number of spikes currently used

		vertex3f trans;
		vertex3f rot;
		vertex3f scale; 
};

#endif
