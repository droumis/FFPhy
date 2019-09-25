/* $Id: Gauss.cc,v 1.8 2008/09/01 21:07:44 chengs Exp $
   
   Sen Cheng, 2004/../..
   
   program description
*/

#include "../aux/hippoIO.h"
#include "../aux/mexAux.h"
//#include "AscentFilter.h"
#include "Gauss.h"

/* ************************************************************
                          class Gauss
   ************************************************************ */

using namespace AFilter;

static const double _2_pi= 0.63661977236758;

Gauss::Gauss()
{
//    nFct= 1; //@@ Gauss cannot handle multiple trajectories, yet!
    alpha= mx= my= Sx= Sy= r= 0;
    alloc= false;
}

Gauss::Gauss(Gauss &)
{
    hippo_Message("function empty");
}

Gauss::~Gauss()
{
    dealloc();
}

void Gauss::dealloc()
{
    if( alloc ) {
        if(alpha) free(alpha);
        if(mx) free(mx);
        if(my) free(my);
        if(Sx) free(Sx);
        if(Sy) free(Sy);
        if(r) free(r);
    }
    startindex= endindex= T= 0;
    alpha= mx= my= Sx= Sy= r= 0;
    alloc= false;
}

void Gauss::init(int _startindex, int _endindex, int _T, double *_x, double *_y, int _nid, traj_type *_id)
{
    VDynamics2d::init(_startindex, _endindex, _T, _x, _y, _nid, _id);
    hippo_Message("Gauss cannot handle multiple trajectories, yet");
}


void Gauss::setInitialEst(int t) const
{
    hippo_Empty;
    // initialize parameter values
//     for (int t= 0; t < data->ntimesteps; t++) {
// 	alpha[t]= 3;
// //	alpha[t]= 0.39894;  // 1/sqrt(2 pi)
// 	mx[t]= my[t]= 10;
// 	Sx[t]= Sy[t]= 20;
// 	r[t]= 0;
//     }

//     if(!K) hippo_CAlloc(K, 4*T, double);
    
//     // calculate sigma^-1 for each timestep
//     double det;
//     int index;
//    for(int t=startindex; t <= endindex; t++) {
// 	index= 4*t;
// 	det= calcMatrixDet(sigma+index);
// 	hippo_Assert(det, "sigma[t] is singular");
	
// 	calcMatrixInv(sigma+index, K+index);
//     }
}

double Gauss::eval(int t, double x, double y, traj_type traj) const
{
    hippo_Assert(t >= startindex && t <= endindex, "t out of bounds");
//    double d[2]= {x-mx[t], 0};
    double d[2]= {x-mx[t], y-my[t]};
    double R= _2_pi*atan(r[t]);
    double quad= d[0]*d[0]/(Sx[t]*Sx[t])- 2*R*d[0]*d[1]/(Sx[t]*Sy[t])
	+ d[1]*d[1]/(Sy[t]*Sy[t]);
    hippo_Assert(quad >= 0, "quad < 0 is impossible");
//    hippo_Print(d[0]);
//    hippo_Print(d[1]);
//    hippo_Print(Sx[t]);
//    hippo_Print(Sy[t]);
//    hippo_Print(mx[t]);
//    hippo_Print(my[t]);
//    hippo_Print(exp(alpha[t]- 0.5*quad/(1-R*R)));
    return exp(alpha[t]- 0.5*quad/(1-R*R));
}

// double Gauss::evalPartial(int t, double x, double y, 
// 			       traj_type traj, double partial[]) const
// {
//    hippo_Assert(t >= startindex && t <= endindex, "t out of bounds");
//     double d[2]= {x-mx[t], y-my[t]};
//     double R= _2_pi*atan(r[t]);
//     double _Sx= 1/Sx[t], _Sy= 1/Sy[t], _R2= 1/(1-R*R);
//     double _Sx2= _Sx*_Sx, _Sy2= _Sy*_Sy;
    
//     double quad= d[0]*d[0]*_Sx2 -2*R*d[0]*d[1]*_Sx*_Sy + d[1]*d[1]*_Sy2;
//     hippo_Assert(quad >= 0, "quad < 0 is impossible");
//     double lambda= exp(alpha[t]- 0.5*quad*_R2);
    
//     partial[0]= lambda;
//     partial[1]= (d[0]*_Sx2 -R*d[1]*_Sx*_Sy)*_R2*lambda;
//     partial[2]= (d[1]*_Sy2 -R*d[0]*_Sx*_Sy)*_R2*lambda;
    
//     partial[3]= (d[0]*d[0]*_Sx2*_Sx - R*d[0]*d[1]*_Sx2*_Sy)*_R2*lambda;
//     partial[4]= (d[1]*d[1]*_Sy2*_Sy - R*d[0]*d[1]*_Sx*_Sy2)*_R2*lambda;
//     partial[5]= ((1+R*R)*d[0]*d[1]*_Sx*_Sy -R*d[0]*d[0]*_Sx2 
// 		 -R*d[1]*d[1]*_Sy2)*_R2*_R2*lambda*_2_pi/(1+r[t]*r[t]);
//     return lambda;
// }

double Gauss::evalGrad(double x, double y, 
				  traj_type traj, double partial[]) const
{
    double d[2]= {x-mx[tCurr], y-my[tCurr]};
    double R= _2_pi*atan(r[tCurr]);
    double Sx1= Sx[tCurr], Sy1= Sy[tCurr], R2= 1-R*R;
    double Sx2= Sx1*Sx1, Sy2= Sy1*Sy1, Sxy= Sx1*Sy1;
    
    double quad= d[0]*d[0]/Sx2 -2*R*d[0]*d[1]/Sxy + d[1]*d[1]/Sy2;
    hippo_Assert(quad >= 0, "quad < 0 is impossible");
    double lambda= exp(alpha[tCurr]- 0.5*quad/R2);
    
    partial[0]= 1;
    partial[1]= (d[0]/Sx2 -R*d[1]/Sxy)/R2;
    partial[2]= (d[1]/Sy2 -R*d[0]/Sxy)/R2;
    
    partial[3]= (d[0]*d[0]/Sx2/Sx1 - R*d[0]*d[1]/Sx2/Sy1)/R2;
    partial[4]= (d[1]*d[1]/Sy2/Sy1 - R*d[0]*d[1]/Sx1/Sy2)/R2;
    partial[5]= ((1+R*R)*d[0]*d[1]/Sx1/Sy1 -R*d[0]*d[0]/Sx2 
		 -R*d[1]*d[1]/Sy2)/R2/R2*_2_pi/(1+r[tCurr]*r[tCurr]);

//     if( t== 3055) {
// 	cout << "evalLogPartial, t= " << t << "\t";
// 	cout << alpha[t] << "\t";
// 	cout << mx[t] << "\t";
// 	cout << my[t] << "\t";
// 	cout << Sx[t] << "\t";
// 	cout << Sy[t] << "\t";
// 	cout << R  << "\n";
// 	cout << d[0]  << "\t";
// 	cout << d[1]  << "\t";
// 	cout << partial[2]  << "\t";
// 	cout << partial[4]  << "\t";
// 	cout << lambda  << "\n";
//     }
    return lambda;
}

// double Gauss::evalLogPartial(int t, double x, double y, 
// 				  traj_type traj, double partial[]) const
// {
//    hippo_Assert(t >= startindex && t <= endindex, "t out of bounds");
//     double d[2]= {x-mx[t], y-my[t]};
//     double R= _2_pi*atan(r[t]);
//     double _Sx= 1/Sx[t], _Sy= 1/Sy[t], _R2= 1/(1-R*R);
//     double _Sx2= _Sx*_Sx, _Sy2= _Sy*_Sy;
    
//     double quad= d[0]*d[0]*_Sx2 -2*R*d[0]*d[1]*_Sx*_Sy + d[1]*d[1]*_Sy2;
//     hippo_Assert(quad >= 0, "quad < 0 is impossible");
//     double lambda= exp(alpha[t]- 0.5*quad*_R2);

//     partial[0]= 1;
//     partial[1]= (d[0]*_Sx2 -R*d[1]*_Sx*_Sy)*_R2;
//     partial[2]= (d[1]*_Sy2 -R*d[0]*_Sx*_Sy)*_R2;
    
//     partial[3]= (d[0]*d[0]*_Sx2*_Sx - R*d[0]*d[1]*_Sx2*_Sy)*_R2;
//     partial[4]= (d[1]*d[1]*_Sy2*_Sy - R*d[0]*d[1]*_Sx*_Sy2)*_R2;
//     partial[5]= ((1+R*R)*d[0]*d[1]*_Sx*_Sy -R*d[0]*d[0]*_Sx2 
// 		 -R*d[1]*d[1]*_Sy2)*_R2*_R2*_2_pi/(1+r[t]*r[t]);
//     return lambda;
// }

void Gauss::updateDynamics(int tStart, int tEnd, double diff[])
{
//    const double tiny= 1e-5;


    alpha[tEnd]= alpha[tStart] + diff[0];
    mx[tEnd]= mx[tStart] +  diff[1];
    my[tEnd]= my[tStart] +  diff[2];
    Sx[tEnd]= Sx[tStart] +  diff[3];
    Sy[tEnd]= Sy[tStart] +  diff[4];
    r[tEnd] =  r[tStart] +  diff[5];

//     if(r[tEnd] <= -1) {
// 	cout << "-* r= " << r[tEnd] << "\t";
// 	r[tEnd]= -1 + tiny;
//     } else if (r[tEnd] >= 1) {
// 	cout << "+* r= " << r[tEnd] << "\t";
// 	r[tEnd]= 1 - tiny;
//     }

//     if(tEnd==0) {

      for(int i=0; i < 6; i++)
	  if(fabs(diff[i]) > 5) {
	      cout << "t= " << tStart << "\t";
	      cout << alpha[tStart] << "\t";
	      cout << mx[tStart] << "\t";
	      cout << my[tStart] << "\t";
	      cout << Sx[tStart] << "\t";
	      cout << Sy[tStart] << "\t";
	      cout <<  _2_pi*atan(r[tStart]) << "\n";
	      
	      cout << "t= " << tEnd << "\t";
	      cout << alpha[tEnd] << "\t";
	      cout << mx[tEnd] << "\t";
	      cout << my[tEnd] << "\t";
	      cout << Sx[tEnd] << "\t";
	      cout << Sy[tEnd] << "\t";
	      cout <<  _2_pi*atan(r[tEnd]) << "\n";
	  }
}

// TFilter*  Gauss::allocFilter(char algo[]) 
// {
    
//     if (strcmp(algo,"SteepestAscent")==0) {
// 	return new AFilter::AscentFilter(data, this);
//     } else {
// 	MX_Error("Algorithm \'%s\' has not been implemented for AFilter::Gauss", algo);
// 	return 0;
//     }
// };

void  Gauss::setMexInput(const mxArray *in) 
{
    dealloc();
    VDynamics2d::setMexInput(in);
    int tmp;
    MX_FieldAssign(in, "alpha", alpha, T);
    
//    hippo_Print(alpha);

    MX_FieldAssign(in, "mx", mx, tmp);
    hippo_Assert(tmp==T, "alpha and mx should have same number of steps");
    
    MX_FieldAssign(in, "my", my, tmp);
    hippo_Assert(tmp==T, "alpha and my should have same number of steps");
    
    MX_FieldAssign(in, "Sx", Sx, tmp);
    hippo_Assert(tmp== T, "alpha and Sx should have same length");

    MX_FieldAssign(in, "Sy", Sy, tmp);
    hippo_Assert(tmp== T, "alpha and Sy should have same length");

    MX_FieldAssign(in, "r", r, tmp);
    hippo_Assert(tmp== T, "alpha and r should have same length");

    nFct= 1;
//    hippo_Print(this);
};

mxArray*  Gauss::getMexOutput() 
{
    mxArray *out= VDynamics2d::getMexOutput();
    mxAddField(out, "name");
    mxAddField(out, "alpha");
    mxAddField(out, "mx");
    mxAddField(out, "my");
    mxAddField(out, "Sx");
    mxAddField(out, "Sy");
    mxAddField(out, "r");
    
    mxArray *mxtmp;
    double *ptr;

    mxtmp= mxCreateString("Gauss");
    hippo_Assert(mxtmp, "couldn't allocate space for string");
    mxSetField(out, 0, "name", mxtmp);
    
    MX_CreateDoubleAssign(mxtmp, 1, T, ptr);
    memcpy(ptr, alpha, T* sizeof(double));
    mxSetField(out, 0, "alpha", mxtmp);
    
    MX_CreateDoubleAssign(mxtmp, 1, T, ptr);
    memcpy(ptr, mx, T* sizeof(double));
    mxSetField(out, 0, "mx", mxtmp);

    MX_CreateDoubleAssign(mxtmp, 1, T, ptr);
    memcpy(ptr, my, T* sizeof(double));
    mxSetField(out, 0, "my", mxtmp);
    
    MX_CreateDoubleAssign(mxtmp, 1, T, ptr);
    memcpy(ptr, Sx, T* sizeof(double));
    mxSetField(out, 0, "Sx", mxtmp);
    
    MX_CreateDoubleAssign(mxtmp, 1, T, ptr);
    memcpy(ptr, Sy, T* sizeof(double));
    mxSetField(out, 0, "Sy", mxtmp);
    
    MX_CreateDoubleAssign(mxtmp, 1, T, ptr);
    memcpy(ptr, r, T* sizeof(double));
    mxSetField(out, 0, "r", mxtmp);
    
    return out;
};

