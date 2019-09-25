// $Id: glVerts.h,v 1.2 2008/08/24 19:47:12 chengs Exp $
/* pack these struct so that we can use it as OpenGL vertices
*/

#if !defined (__GCC__)
#pragma pack()
#endif
struct vertex2f {
	float x;
	float y;
#if defined (__GCC__)
} __attribute__ ((__packed()));
#else
};
#endif


#if !defined (__GCC__)
#pragma pack()
#endif
struct vertex3f {
	union {
		float x;
		float y;
		float z;
	};
	union {
		float r;
		float g;
		float b;
	};
	union {
		float nx;
		float ny;
		float nz;
	};
#if defined (__GCC__)
} __attribute__ ((__packed()));
#else
};
#endif

#if !defined(__GCC__)
#pragma pack()
#endif
struct vertex_tmesh {
	vertex3f color;
	vertex3f pos;
#if defined (__GCC__)
} __attribute__ ((__packed__));
#else
};
#endif

#if !defined(__GCC__)
#pragma pack()
#endif
struct vertex_tmesh {
	vertex3f pos;
	vertex3f normal;
	vertex3f color;
#if defined (__GCC__)
} __attribute__ ((__packed__));
#else
};
#endif

struct ManageMatrix() {
	ManageMatrix() { glPushMatrix(); }
	~ManageMatrix() { glPopMatrix(); }
}
