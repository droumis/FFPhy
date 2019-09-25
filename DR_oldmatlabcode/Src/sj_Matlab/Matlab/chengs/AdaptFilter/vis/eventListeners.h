// $Id: eventListeners.h,v 1.2 2008/08/24 19:47:12 chengs Exp $
#ifndef __EVENTLISTENERS_H__
#define __EVENTLISTENERS_H__

struct TrajMeshEventListener
{
    virtual ~TrajMeshEventListener() {};
	virtual void trajMeshEventHandle()=0;
};

struct MovieEventListener
{
    virtual ~MovieEventListener() {};
	virtual void movieEventHandle()=0;
};

struct MTimeEventListener
{
    virtual ~MTimeEventListener() {};
	virtual void mTimeEventHandle()=0;
};

#endif
