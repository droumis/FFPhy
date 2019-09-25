// $Id: VPlot.h,v 1.1 2008/08/24 19:47:09 chengs Exp $
#ifndef AF_VPlot_h
#define AF_VPlot_h


namespace AFilter {

class VPlot {

public:
    VPlot() {};
    virtual ~VPlot() {};

    /*
     * The following require a current OpenGL context.
     */
    virtual void init() =0;
    virtual void draw() =0;
    virtual void compute(int time) =0;

};

}; // namespace AFilter 

#endif // #ifdef AF_VPlot_h
