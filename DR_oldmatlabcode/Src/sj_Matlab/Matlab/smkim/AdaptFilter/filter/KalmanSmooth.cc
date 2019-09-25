/* $Id: KalmanSmooth.cc,v 1.8 2008/08/24 19:47:01 chengs Exp $
   
   Sen Cheng, Thu Jun  1 15:44:49 PDT 2006

   program description
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../aux/hippoIO.h" 
#include "../aux/mexAux.h" 
#include "../aux/TData.h" 
#include "../model/AdaptModel.h"
#include "KalmanSmooth.h"

/* ************************************************************
                          class KalmanSmooth
   ************************************************************ */


using namespace AFilter;


KalmanSmooth::KalmanSmooth(TData *_data, AdaptModel *_model)
: data(_data), nState(0), nTotal(0), Q(0), model(_model), var(0)
{
//    hippo_Mark;
    scov= 0;
    nOcc= 0;
    x1= Qt= 0;
    Q1= 0;
}

KalmanSmooth::~KalmanSmooth()
{
    if(var) delete var;
    if(scov) delete[] scov;
    if(nOcc) delete[] nOcc;
//    if(Q) delete[] Q;
}

void KalmanSmooth::init()
{
    hippo_Assert(data, "data object not set");
    hippo_Assert(data->getNValidSpikes(), "data contains no spikes");
    hippo_Assert(model, "model object not set");
    hippo_Assert(nTotal== model->getNAllParam(), "covariance should have same number of elements as model parameters");

    nState= model->getNParam();
    nTotal= model->getNAllParam();
//    for(int i=0; i < nTotal; i++)  {
//        hippo_Assert(Q[i] > 0, "state transition covariance should be positive");
//    }
//    hippo_Mark;
    var= model->link_copy();
    var->setUnreal();
//    hippo_Print(nTotal);
//    hippo_Print(nState);
    scov= new double[nTotal];
    nOcc= new int[nTotal];
//    Q= new double[nTotal];
}

void KalmanSmooth::runFilter()
{
    int tStart= model->getStartindex(),
        tEnd= model->getEndindex();
    model->setInitialEst(tStart);
    var->setAllParam(tStart,Q1); 
//    hippo_Print(Qt);
//    hippo_Print(Q1);
//    hippo_Print(x1);
    
//    for(int i=0; i < nTotal; i++) { nOcc[i]= 0; Q[i]= Qt; }
    for(int i=0; i < nTotal; i++) { nOcc[i]= 0; }
    int *a= new int[nState];
    double *x= new double[nState];
    for(int t= tStart; t < tEnd; ++t) {
        if(data->fieldID[t]>= 0) {
            model->getParam(t, x, a); 
            for(int i=0; i < nState; i++) if(a[i]>= 0) ++nOcc[a[i]];
        }
    }
//    cerr << endl;
//    double maxocc= -1e6;
//    for(int i=0; i < nTotal; ++i) {
//        cerr << nOcc[i] << "\t";
//        if(nOcc[i] > maxocc) maxocc= nOcc[i];
//    }
//    cerr << endl;
//    cerr << "maxOcc= " << maxocc << endl;
    delete[] a;
    delete[] x;

//    runPass(1); // temporal parameters fixed 
//    runPass(0); // spatial parameters fixed 
    runPass(-1); // no fixing
}


void KalmanSmooth::runPass(int fixing)
{
    // run the adaptive filter algorithm until convergence or limit of
    // iterations is reached
    const double            dT= data->deltaTime;
    double *V= new double[nState];
    double *Vd= new double[nState];
    double *x= new double[nState];
    int N= model->getNAllParam();
    double *VT= new double[N];
//    double *xtmp= new double[nState];
    int tStart= model->getStartindex(),
        tEnd= model->getEndindex();
    int *a= new int[nState];
    Grad g(nState);
    hippo_Print(dT);
    hippo_Print(tStart);
    hippo_Print(tEnd);
    for(int i=0; i < nTotal; i++) { scov[i]= 0; }

//    cerr << "\n\tKalmanSmooth forward pass\n";
    // forward pass -- filtering
    for(int t= tStart; t < tEnd; t++) {
        model->moveTo(t);
        for(int i=0; i < nState; i++) { x[i]= Vd[i]= 0; }
        if(data->fieldID[t]>= 0 && var->getParam(t, V, a)) {
            int dNt= data->dN[t+1];

            /*
            hippo_Print(N);
    var->getAllParam(tEnd, VT); 
    cerr << "V at t= " << t << "\n";
    for(int i=0; i<nTotal; ++i) cerr << VT[i] << ' ';
    cerr << "\n---\n";
    for(int i=0; i<nState; ++i) cerr << "a[i]= " << a[i] << ", V= " << V[i] << '\n';
    cerr << endl;
    */

            model->evalNoMove(t, g);

            cerr << ", f=" << g.f << "\t";
            for(int i=0; i < nState; i++) {
            
                if(a[i] < 0) continue;

//                if(a[i] >= nTotal) { 
//                    cerr << "i= " << i << ", a[i]= " << a[i]
//                        << ", nTotal= " << nTotal ;
//                    cerr << endl;
//                    exit(0); 
//                }
                double inno= (dNt-dT*g.f);
                double old= V[i];
//                double Ri= dT*g.ddf[i]-dNt*g.ddlogf[i];
                double Ri= sq(g.dlogf[i])*dT*g.f - g.ddlogf[i]*inno;
                V[i]= 1/(1/(old+Q[a[i]]) + Ri);
                if(V[i]<.5) V[i]= .5;
//                V[i]= (Ri > 0) ? 1/(1/(old+Q[a[i]]) + Ri) : old+Q[a[i]];
//
                Vd[i]= V[i]-old;
//                x[i]= V[i]*(dNt*g.dlogf[i]-dT*g.df[i]);
                x[i]= V[i]*g.dlogf[i]*inno;
//                if(Ri>.1) cerr << "t= " << t << ", a[i]= " << a[i] << ", Ri=" << Ri << ", old V= " << old << ", V= " << V[i] << '\n';
//                cerr << V[i] << '\t';
//                cerr << Ri << "\t";
                if(x[i] > 100 || Ri==0) {
//                if(!(Ri > -tiny) || !(V[i] > -tiny) || old < 0) {
//                if(dNt) {
//                if(0) {
                    cerr << "t= " << t << ":" 
                        << ", Q=" << Q[a[i]]
//                        << ", a[i]=" << a[i]
//                        << ", xtmp[i]=" << xtmp[i] 
                        << ", x[i]=" << x[i] 
//                        << ", dNt=" << dNt 
                        << ", V[i]=" << old 
                        << ", nV[i]=" << V[i] 
//                        << ", dV[i]=" << Vd[i] 
                        << ", R= " << 1/Ri 
//                        << ", 1/R= " << Ri 
//                << ", t1= " << sq(g.dlogf[i])*dT*g.f << ", t2= " << - g.ddlogf[i]*inno
                        << ", logf=" << g.logf
                        << ", df[i]=" << g.df[i] 
                        << ", ddf[i]=" << g.ddf[i] 
                        << ", dlogf[i]=" << g.dlogf[i] 
                        << ", ddlogf[i]=" << g.ddlogf[i] 
                        << ", f=" << g.f << "\t"
                        << ", inno=" << inno << "\t"
                        << endl;
                    exit(0);
                } 
                if(!(Ri > -tiny) || !(V[i] > -tiny) || old < 0) throw 1;

                hippo_Assert(Ri>-tiny, "Ri has to be positive");
            }  // loop over components of state vector
        } // valid timestep 
        model->updateModel(t, t+1, x, fixing);
        var->updateModel(t, t+1, Vd, fixing);
    } // for(int t= tStart; t < tEnd; t++) 

//            return; //@@

    // backward pass -- smoothing
    double *xT= new double[N];
    double maxV=0, minV= DBL_MAX;

//    cerr << "\tKalmanSmooth backward pass\n";

    model->getAllParam(tEnd, xT); 
    var->getAllParam(tEnd, VT); 

    for(int t= tEnd-1; t >= tStart; t--) {
        for(int i=0; i < nState; i++) { x[i]= V[i]= 0; }
        if(data->fieldID[t] >= 0 && model->getParam(t, x, a)) {
            var->getParam(t, V, a);

            // calc state estimate, its variance and covariances
            for(int i=0; i < nState; i++) {
                if(a[i]<0) continue; 

                double J= V[i]/(V[i]+Q[a[i]]);
                double K= Q[a[i]]/(V[i]+Q[a[i]]);
                // x= diff between xT(t) - xT(t+1) 
                double nx= K*x[i] + J*xT[a[i]];
                x[i]= nx- xT[a[i]];
                xT[a[i]]= nx;
                double nV= K*V[i] + J*VT[a[i]]*J;
//                    if(nV<-.01) 
//                        cerr << "t= " << t << "\ti= " << i << "\tnV= " << nV << endl;
                hippo_Assert(nV > -0.1, "variance negative");
                if(nV<0) nV= 0;

//                scov[a[i]]+= x[i]*x[i]+ nV+ (1-2*J)*VT[a[i]] ;
                scov[a[i]]+= x[i]*x[i]+ K*V[i] + sq(K)*VT[a[i]];

                if(0) {
                    cerr << "t= " << t 
//                    << ", a[i]= " << a[i] 
//                        << ", x[i]=" << x[i] 
//                        << ", xT[a[i]]=" << xT[a[i]] 
//                        << ", K=" << K 
//                        << ", J=" << J 
                        << ", V=" << V[i] 
                        << ", nV=" << nV
                        << ", VT=" << VT[a[i]]
//                        << ", dV=" << nV-VT[a[i]]
//                        << ", nx1= " << K*x[i] + J*xT[a[i]]
//                        << ", nx2= " << x[i] + J*(xT[a[i]]-x[i])
//                        << ", nV1= " << K*V[i] + J*VT[a[i]]*J
//                        << ", nV2= " << V[i] + sq(J)*(VT[a[i]] - V[i] - Q[a[i]])
                        << ", cov1= " << x[i]*x[i]+ nV+ (1-2*J)*VT[a[i]] 
//                        << ", cov2= " << x[i]*x[i]+ K*V[i] + sq(K)*VT[a[i]]
//                        << ", x^2= " << x[i]*x[i]
//                        << ", KVt= " << K*V[i] 
//                        << ", K^2Vt+1= " << sq(K)*VT[a[i]]
                        << '\n';
                }

                V[i]= nV-VT[a[i]];
                VT[a[i]]= nV;

                if(maxV<nV) maxV= nV;
                if(nV < minV) minV= nV;

            }
        } // valid trajectory update
        model->updateModel(t+1, t, x, fixing);
        var->updateModel(t+1, t, V, fixing);
    }

    hippo_Print(maxV);
    delete[] a;
    delete[] x;
    delete[] V;
    delete[] Vd;
    delete[] xT;
    delete[] VT;
}


void KalmanSmooth::close()
{
}

void KalmanSmooth::setMexInput(const mxArray* MexParam)
{
//    MX_FieldAssign(MexParam, "Q", Q, nState);
//    MX_FieldScalarDefault(MexParam, "x1", x1, 1, double);
    MX_FieldScalarDefault(MexParam, "Q1", Q1, 1, double);
//    MX_FieldScalarDefault(MexParam, "Qt", Qt, 1, double);
    MX_FieldAssign(MexParam, "Qt", Q, nTotal);

//    cerr << "Q= ";
//    for(int i=0; i<nTotal; ++i) cerr << Q[i] << '\t';
//    cerr << endl;
//    hippo_Mark;

//    Q1= new double[nTotal];
//    for(int i=0; i<1144; ++i) Q1[i]= 1e-3;
//    for(int i=1144; i<1160; ++i) Q1[i]= 1e-3;
}

mxArray* KalmanSmooth::getMexOutput()
{
    /*
    const char *fields[]= { "var" };
    const int nFields= sizeof(fields)/sizeof(fields[0]);
    mxArray *out= mxCreateStructMatrix(1,1, nFields,  fields);

    hippo_Mark;
    mxSetField(out, 0, "var", var->getMexOutput());

    hippo_Mark;
    return out;
    */
    return var->getMexOutput();
}
