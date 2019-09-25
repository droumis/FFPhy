// $Id: eventListeners.h,v 1.1.1.1 2006/04/24 14:21:48 chengs Exp $
#ifndef __EVENTLISTENERS_H__
#define __EVENTLISTENERS_H__

struct TrajMeshEventListener
{
	virtual void trajMeshEventHandle()=0;
};

struct MovieEventListener
{
	virtual void movieEventHandle()=0;
};

struct MTimeEventListener
{
	virtual void mTimeEventHandle()=0;
};

#endif
