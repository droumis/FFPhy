/* $Id: Grad.h,v 1.2 2006/08/24 15:52:09 chengs Exp $
   
   Sen Cheng, Tue Jun 13 17:34:22 PDT 2006
   
*/

#ifndef Grad_H
#define Grad_H

#include "../aux/numerics.h" 

namespace AFilter {

/* ************************************************************
                          struct Grad
   ************************************************************ */

struct Grad {
    int n;
    double f, logf, *df, *dlogf, *ddf, *ddlogf;

    int i, p[2];
    double y[2], logy[2];
    bool valid[2];
    struct { bool f, logf, df, dlogf, ddf, ddlogf; } zero[2];

    bool *ok;

    Grad(int _n) {
        n= _n;
        p[1]= p[2]= 0;
        df=     new double[n];
        ddf=    new double[n];
        dlogf=  new double[n];
        ddlogf= new double[n];
        ok= new bool[n];
        reset();
    }

    ~Grad() {
        delete[] df;
        delete[] ddf;
        delete[] dlogf;
        delete[] ddlogf;
        delete[] ok;
    }

    void reset() {
        for(int i=0; i<n; i++) {
            df[i]= ddf[i]= dlogf[i]= ddlogf[i]= 0;
            ok[i]= true;
        }
        for(int j=1; j<2; j++) {
            y[j]= logy[j]= 0;
            valid[j]= 0;
            zero[j].df= zero[j].ddf= zero[j].dlogf= zero[j].ddlogf= true;
        }
    }

    void markInvalid(int m) {
        for(int i=0; i<m; i++) ok[i]= false;
    }

    const Grad& operator+=(int i) {
        df+= i;
        ddf+= i;
        dlogf+= i;
        ddlogf+= i;
        return *this;
    }

    const Grad& operator-=(int i) {
        df-= i;
        ddf-= i;
        dlogf-= i;
        ddlogf-= i;
        return *this;
    }

}; // struct Grad


} // namespace AFilters

#endif   // Grad_H
