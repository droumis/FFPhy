/* $Id: Operators.cc,v 1.7 2008/06/19 22:49:26 chengs Exp $
   
   Sen Cheng, Thu Jun  8 21:15:17 PDT 2006
   
   Collection of operators.
*/

#include "../aux/hippoIO.h" 
#include "../aux/mexAux.h" 
#include "../aux/numerics.h" 
#include "AdaptModel.h" 
#include "Operators.h" 

/* ************************************************************
                          class VOp
   ************************************************************ */

static const double maxLog= log(DBL_MAX);

using namespace AFilter;

VOp* VOp::newOperator(const char *name)
{
    VOp *o;
    if (strcmp(name,"Id")==0) {
        o= new OpId;
    } else if (strcmp(name,"Rectify")==0) {
        o= new OpRectify;
    } else if (strcmp(name,"Square")==0) {
        o= new OpSquare;
    } else if (strcmp(name,"Exp")==0) {
        o= new OpExp;
    } else if (strcmp(name,"Exp2")==0) {
        o= new OpExp2;
    } else {
        hippo_Print(name);
        MX_Error("unknown operator: '%s' ", name);
        o= 0;
    }
//    hippo_Print(o->getOpName());
    return o;
}

VOp* VOp::newOperator(const mxArray *in)
{
    VOp* o= 0;
    mxArray *mxtmp= mxGetField(in, 0, "name");
    hippo_Assert(mxtmp, "field 'name' not defined for operator");
    char *tmpstring;
    MX_StringAllocAssign(mxtmp, tmpstring);
//        hippo_Print(tmpstring);
    o= VOp::newOperator(tmpstring);
    mxFree(tmpstring);
    o->setMexInput(in);
    return o;
}

mxArray* VOp::getMexOutput() { 
    mxArray *out= mxCreateStructMatrix(1,1, 0,  0); 
    mxAddField(out, "name");
    mxSetField(out, 0, "name", mxCreateString(getOpName()));
//    mxSetField(out, 0, "name", mxCreateString("Id"));//@@
    return out;
};

/* ************************************************************
                          class OpId
   ************************************************************ */

void OpId::oplog(double &x) const 
{ 
    if(x < tiny) x= tiny;
    x= log(x);
}

void OpId::oplog(double &x, int n, double gradx[]) const 
{ 
    if(x < tiny) x= tiny;
    for(int i=0; i<n; i++) gradx[i]/= x;
    x= log(x);
}

void OpId::op(int n, Grad &g) const
{
    double x= g.f;
    hippo_Assert(x>0, "argument of log must be positive");
    g.logf= log(x);
    for(int i=0; i<n; i++) {
        g.ddf[i]= 0;
        g.dlogf[i]= g.df[i]/x;
        g.ddlogf[i]= - g.dlogf[i]*g.dlogf[i];
    }
}

/* ************************************************************
                          class OpRectify
   ************************************************************ */

const double small= 1e-3;

void OpRectify::oplog(double &x) const
{
    if(x < small) x= small;
    x= log(x);
}

void OpRectify::oplog(double &x, int n, double gradx[]) const
{
    if(x < small) x= small;
    for(int i=0; i<n; i++) gradx[i]/= x;
    x= log(x);
}

void OpRectify::op(int n, Grad &g) const
{
    double x= g.f;
    if(x < small) { x= small; g.f= 0; }
    g.logf= log(x);
    for(int i=0; i<n; i++) {
        g.ddf[i]= 0;
        g.dlogf[i]= g.df[i]/x;
        g.ddlogf[i]= - g.dlogf[i]*g.dlogf[i];
    }
}

/* ************************************************************
                          class OpSquare
   ************************************************************ */

void OpSquare::op(double &x, int n, double gradx[]) const
{
    for(int i=0; i<n; i++) gradx[i]*= 2*x;
    x*= x;
}

void OpSquare::oplog(double &x) const
{
    if(x < tiny) x= tiny;
    x= 2*log(x);
}

void OpSquare::oplog(double &x, int n, double gradx[]) const
{
    if(x < tiny) x= tiny;
    for(int i=0; i<n; i++) gradx[i]/= .5*x;
    x= 2*log(x);
}

void OpSquare::op(int n, Grad &g) const
{
    double x= g.f;
    if(0 < x && x < tiny) x= tiny; 
    if(0 > x && x > -tiny) x= -tiny; 
    for(int i=0; i<n; i++) {
        g.dlogf[i]= 2*g.df[i]/x;
        g.ddlogf[i]= - .5*g.dlogf[i]*g.dlogf[i];
        g.ddf[i]= 2*g.df[i]*g.df[i];
        g.df[i]*= 2*g.f;
    }
    g.f*= g.f;
    g.logf= log(x);
}

/* ************************************************************
                          class OpExp
   ************************************************************ */

void OpExp::op(double &x, int n, double gradx[]) const 
{ 
    x= exp(x);
    for(int i=0; i<n; i++) gradx[i]*= x;
}

void OpExp::op(int n, Grad &g) const
{

//    if(!(g.f < maxLog)) { cerr << g.f << endl; }
    hippo_Assert(g.f < maxLog, "x value too large, exp(x) will be inf.");

//    hippo_Assert((g.f>=0) || (g.f<0), "function value not valid");
    double y= exp(g.f);
    g.logf= g.f;
    for(int i=0; i<n; i++) {
//        if(!(g.df[i]>=0) && !(g.df[i]<0)) {
//            cerr << "g.df[i]=" << g.df[i] << endl;
//        }
//        hippo_Assert((g.df[i]>=0) || (g.df[i]<0), "derivative value not valid");
        g.dlogf[i]= g.df[i];
        g.ddlogf[i]= 0;
        g.ddf[i]= g.df[i]*g.df[i]*y;
        g.df[i]*= y;
    }
    g.f= y;
}

/* ************************************************************
                          class OpExp2
   ************************************************************ */
//  f(x)= exp(a+bx)

void OpExp2::op(double &x, int n, double gradx[]) const 
{ 
    double y= exp(a+b*x);
    for(int i=0; i<n; i++) gradx[i]*= b*y;
    x= y;
}

void OpExp2::op(int n, Grad &g) const
{
    hippo_Assert((g.f>=0) || (g.f<0), "function value not valid");
    g.logf= a+b*g.f;
    double y= exp(g.logf);
    for(int i=0; i<n; i++) {
        hippo_Assert((g.df[i]>=0) || (g.df[i]<0), "derivative value not valid");
        g.dlogf[i]= b*g.df[i];
        g.ddlogf[i]= 0;
        g.ddf[i]= sq(b)*sq(g.df[i])*y;
        g.df[i]*= b*y;
    }
    g.f= y;
}

double OpExp2::lik(AdaptModel *m, AdaptModel *var) const
{
    int start= m->getStartindex();
    int end= m->getEndindex();
    const TData *data= m->getData();
    double dt= data->getTimestep();
    double atmp,  expa, l, x, V, lambda;
    atmp= l= 0;
    for(int t= start; t < end; ++t) {
        x= m->evalCore(t);
        V= var->evalCoreVar(t);
        atmp+= exp(b*x +0.5*b*b*V);
    }
    expa= data->getNValidSpikes()/(dt*atmp);
    atmp= log(expa);
    for(int t= start; t < end; ++t) {
        x= m->evalCore(t);
        V= var->evalCoreVar(t);
        lambda= m->eval(t);
        l+= data->dN[t]*(atmp+ b*x) -dt*expa*exp(b*x +0.5*b*b*V);
    }
    return l;
}

void OpExp2::reestimate(AdaptModel *m, AdaptModel *var)
{
//    int start= m->getStartindex();
//    int end= m->getEndindex();
//    const TData *data= m->getData();
//    double dt= data->getTimestep();
    double aD, bN, bD;
//    double X, V;
    aD= bN= bD= 0;


    // find b by iteration (instable)
    /*
   double lambda;
    for(int i=0; i<20; ++i) {
        for(int t= start; t < end; ++t) {
            x= m->evalCore(t);
            V= var->evalCoreVar(t);
            lambda= m->eval(t);
            bN+= (data->dN[t] -lambda*dt)*x;
            bD+= lambda*V;
        }
        b= bN/(dt*bD);
        cerr << "b= " << b;
        cerr << ", num= " << bN;
        cerr << ", den= " << bD;
        cerr << "\n";
    }
    */

//    for(int q=0; q < 100; ++q) {
//        b= q;
//        printf("%f,\t %f\n", b, lik(m, var));
//    }
//    exit(0);

    // Golden section search
    /*
    const double C= 0.3819660; // golden section ratio
    const double R= 1-C; // golden section ratio
    const double tol= 1e-1;
    const int maxit= 30;
    int it= 0;
    double w= 0, x, y, z= 30;
    double fb, fw, fx, fy, fz;
    fb= lik(m, var);

    b=w; fw= lik(m, var);
    b=z; fz= lik(m, var);
    if(fb < fw) b= w;
    // find bracket

    if( fabs(z-b) > fabs(b-w)) {
        x= b; fx= fb;
        y= b+C*(z-b); 
        b=y; fy= lik(m, var);
    } else {
        x= b-C*(b-w);
        b=x; fx= lik(m, var);
        y= b; fy= fb;
    }
    cerr << "w= " << w << ", x= " << x << ", y= " << y << ", z= " << z << "\n";
    hippo_Print(fb);
    hippo_Print(fw);
    hippo_Print(fz);
    while( it < maxit && fabs(z-w) > tol*(fabs(x)+fabs(y)) ) {
        ++it;
        if(fy>fx) {
            w=x; x= y; y= R*x+C*z;
            fx= fy; b=y; fy= lik(m,var);

            cerr << "y= " << y << ", fy= " << fy << "\n";
        } else {
            z=y; y=x; x= R*y+C*w;
            fy= fx; b=x; fx= lik(m, var);
            cerr << "x= " << x << ", fx= " << fx << "\n";
        }
    }
    b= (fx < fy)? x : y;
    */
    hippo_Print(b);
    

    /*
    for(int t= start; t < end; ++t) {
        X= m->evalCore(t);
        V= var->evalCoreVar(t);
        aD+= exp(b*X +0.5*b*b*V);
    }
    a= log(data->getNValidSpikes()/ (dt* aD));
    */
    hippo_Print(a);
}

void OpExp2::setMexInput(const mxArray *in) 
{
    MX_FieldScalarDefault(in, "a", a, -1, double);
    MX_FieldScalarDefault(in, "b", b, -1, double);
//    hippo_Print(a);
//    hippo_Print(b);
}

mxArray* OpExp2::getMexOutput() 
{
    mxArray *out= VOp::getMexOutput();
    mxAddField(out, "a");
    mxSetField(out, 0, "a", mxCreateScalarDouble(a));
    mxAddField(out, "b");
    mxSetField(out, 0, "b", mxCreateScalarDouble(b));
    return out;
}
