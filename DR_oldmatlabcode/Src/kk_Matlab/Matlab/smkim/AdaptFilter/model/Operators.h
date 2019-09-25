/* $Id: Operators.h,v 1.5 2008/06/19 22:49:26 chengs Exp $
   
   Sen Cheng, Thu Jun  8 21:16:27 PDT 2006
   
   program description

   need to call init before use
*/

#ifndef VOp_H
#define VOp_H

#include "../aux/numerics.h" 
#include "../aux/Grad.h" 

namespace AFilter {

class AdaptModel;

/* ************************************************************
                          class VOp
   ************************************************************ */

class VOp {
public:
    VOp() {};
    virtual ~VOp() {};

    virtual const char* getOpName() const =0;

    virtual double op(const double x) const { hippo_Virtual; }
//    virtual void op(double &x) const { hippo_Virtual; }
    virtual void op(double &x, int n, double gradx[]) const { 
        hippo_Virtual;
    }

    virtual void oplog(double &x) const { hippo_Virtual; }
    virtual void oplog(double &x, int n, double gradx[]) const { 
        hippo_Virtual;
    }

    virtual void op(int n, Grad &g) const { hippo_Virtual; }

    virtual double inv(const double x) const { hippo_Virtual; }

    virtual void reestimate(AdaptModel *m, AdaptModel *var) {};

    static VOp* newOperator(const char *name);
    static VOp* newOperator(const mxArray *);

    virtual void setMexInput(const mxArray *in) {};
    virtual mxArray* getMexOutput();
    
}; // class VOp

class OpId : public VOp {
public:
    OpId() {};
    virtual ~OpId() {};

    virtual const char* getOpName() const { return "Id"; }
    virtual double op(const double x) const { return x; }
//    virtual void op(double &x) const { }
    virtual void op(double &x, int n, double gradx[]) const { }

    virtual void oplog(double &x) const;
    virtual void oplog(double &x, int n, double gradx[]) const;

    virtual void op(int n, Grad &g) const;

    virtual double inv(const double x) const { return x; }
}; // class OpId

class OpRectify : public VOp {
public:
    OpRectify() { };
    virtual ~OpRectify() {};

    virtual const char* getOpName() const { return "Rectify"; }
    virtual double op(const double x) const { return (x<0) ? 0:x; }
//    virtual void op(double &x) const {if(x<0) x= 0;};
    virtual void op(double &x, int n, double gradx[]) const{if(x<0) x= 0;};

    virtual void oplog(double &x) const;
    virtual void oplog(double &x, int n, double gradx[]) const;

    virtual void op(int n, Grad &g) const;

    virtual double inv(const double x) const { return (x<0) ? 0:x; }
}; // class OpRectify

class OpSquare : public VOp {
public:
    OpSquare() {};
    virtual ~OpSquare() {};

    virtual const char* getOpName() const { return "Square"; }
    virtual double op(const double x) const { return x*x; }
//    virtual void op(double &x) const { x*= x; }
    virtual void op(double &x, int n, double gradx[]) const;

    virtual void oplog(double &x) const;
    virtual void oplog(double &x, int n, double gradx[]) const;

    virtual void op(int n, Grad &g) const;

    virtual double inv(const double x) const { return sqrt(x); }
}; // class OpSquare

class OpExp : public VOp {
public:
    OpExp() {};
    virtual ~OpExp() {};

    virtual const char* getOpName() const { return "Exp"; }
    virtual double op(const double x) const { return exp(x); }
//    virtual void op(double &x) const { x= exp(x); }
    virtual void op(double &x, int n, double gradx[]) const;

//    virtual void oplog(double &x) const;
//    virtual void oplog(double &x, int n, double gradx[]) const;

    virtual void op(int n, Grad &g) const;

    virtual double inv(const double y) const { return log(y); }
}; // class OpExp

/** f(x)= exp(a+bx)  */
class OpExp2 : public VOp {
private:
    double a, b;

private:
    double lik(AdaptModel *m, AdaptModel *var) const;
public:
    OpExp2() {};
    virtual ~OpExp2() {};

    virtual const char* getOpName() const { return "Exp2"; }
    virtual double op(const double x) const { return exp(a+b*x); }
    virtual void op(double &x, int n, double gradx[]) const;

//    virtual void oplog(double &x) const;
//    virtual void oplog(double &x, int n, double gradx[]) const;

    virtual void op(int n, Grad &g) const;

    virtual double inv(const double y) const { return (log(y)-a)/b; }

    virtual void reestimate(AdaptModel *m, AdaptModel *var);

    virtual void setMexInput(const mxArray *in);
    virtual mxArray* getMexOutput();
}; // class OpExp2

} // namespace AFilters

#endif   // VOp_H
