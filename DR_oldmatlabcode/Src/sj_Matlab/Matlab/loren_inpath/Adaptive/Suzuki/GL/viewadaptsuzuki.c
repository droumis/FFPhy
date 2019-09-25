#include <math.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>
#include "adapt.h"

void init_menus(void);

int no_display;
int showrealfield = 1;

int main(int argc, char **argv) {

	int nxtarg = 0;	

	info.simdata = 0;
   	while (++nxtarg < argc) {
		if (strcmp(argv[nxtarg], "-trial") == 0) {
			info.trialinfofile = (char *) malloc(sizeof(char) * 100);
			strcpy(info.trialinfofile, argv[++nxtarg]);
		}
		else if (strcmp(argv[nxtarg], "-time") == 0) {
			info.timefile = (char *) malloc(sizeof(char) * 100);
			strcpy(info.timefile, argv[++nxtarg]);
		}
		else if (strcmp(argv[nxtarg], "-spike") == 0) {
			info.spikefile = (char *) malloc(sizeof(char) * 100);
			strcpy(info.spikefile, argv[++nxtarg]);
		}
		else if (strcmp(argv[nxtarg], "-esttheta") == 0) {
			info.estthetafile = (char *) malloc(sizeof(char) * 100);
			strcpy(info.estthetafile, argv[++nxtarg]);
		}
		else if (strcmp(argv[nxtarg], "-realtheta") == 0) {
			info.realthetaxfile = (char *) malloc(sizeof(char) * 100);
			strcpy(info.realthetaxfile, argv[++nxtarg]);
			info.realthetatfile = (char *) malloc(sizeof(char) * 100);
			strcpy(info.realthetatfile, argv[++nxtarg]);
			info.simdata = 1;
		}
		else {
		fprintf(stderr, "%s\n", argv[nxtarg]);
			fprintf(stderr,"Usage: viewadaptest -trial trialinfofile -time timefile -spike spikefile\n\t\t-esttheta estthetafile [-realtheta realthetax realthetat]\n");
			exit(0);
		}	
	} 

	/* check to make sure that all necessary arguments have been specified */
	if (info.trialinfofile == NULL) {
		fprintf(stderr, "Error: must specify trialinfofile\n");
		exit(1);
	}
	if (info.timefile == NULL) {
		fprintf(stderr, "Error: must specify timefile\n");
		exit(1);
	}
	if (info.spikefile == NULL) {
		fprintf(stderr, "Error: must specify spikefile\n");
		exit(1);
	}
	if (info.estthetafile == NULL) {
		fprintf(stderr, "Error: must specify estthetafile\n");
		exit(1);
	}
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB); 
	glutInitWindowPosition(0,0);
	glutInitWindowSize(INIT_WINDOW_WIDTH, INIT_WINDOW_HEIGHT);
	glutCreateWindow("Adaptive Estimation of Skewed Gaussian Place Field");
	
	Init();

	/* function callbacks */
	glutDisplayFunc(Display);
	glutIdleFunc(Display);
	glutReshapeFunc(Reshape);
	glutKeyboardFunc(Keyboard);
	init_menus();

	glutSwapBuffers();
	glutMainLoop();

	return 1;
}


void init_menus(void) {
	
	/* Menu ID's */
	int Main_Menu	=	1;
	
	Main_Menu = glutCreateMenu(main_menu);
	glutAddMenuEntry("Backup", 'b');
	glutAddMenuEntry("Speed Up", 'u');
	glutAddMenuEntry("Slow Down", 'd');
	glutAddMenuEntry("Start/Stop", 's');
	glutAddMenuEntry("Restart", 'r');
	glutAddMenuEntry("Restart / Presentation", 'R');
	glutAddMenuEntry("Grid On/Off", 'm');
	glutAddMenuEntry("Restart / Dump PPM", 'p');
	
	glutAttachMenu(GLUT_LEFT_BUTTON);
}
