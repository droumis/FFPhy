#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "ranlib.h"
#include "barsP_funcs.h"
#include "barsP_funcs1.h"
#include "barsP_utils.h"
#include "mex.h"

/* see barslib.h for struct definitions */

Bars_X_Grid *New_X_Grid(int nd){
  Bars_X_Grid *xg;
  xg = (Bars_X_Grid *)malloc(sizeof(Bars_X_Grid));
  (*xg).n = nd;
  (*xg).x = tssdvec(nd);
  return xg;
}

void Free_X_Grid(Bars_X_Grid *xg){
  free((*xg).x);
  free(xg);
}

Bars_BarsParams *New_BarsParams(){
  Bars_BarsParams *bp;
  bp = (Bars_BarsParams *)malloc(sizeof(Bars_BarsParams));
  (*bp).iter_knots_fname = tsscvec(MAX_VAL_LEN);
  (*bp).iter_mu_fname = tsscvec(MAX_VAL_LEN);
  (*bp).iter_mufine_fname = tsscvec(MAX_VAL_LEN);
  (*bp).iter_params_fname = tsscvec(MAX_VAL_LEN);
  (*bp).summ_mu_fname = tsscvec(MAX_VAL_LEN);
  (*bp).summ_mufine_fname = tsscvec(MAX_VAL_LEN);
  (*bp).summ_params_fname = tsscvec(MAX_VAL_LEN);
  return bp;
}

void Free_BarsParams(Bars_BarsParams *bp){
  free((*bp).iter_knots_fname);
  free((*bp).iter_mu_fname);
  free((*bp).iter_mufine_fname);
  free((*bp).iter_params_fname);
  free((*bp).summ_mu_fname);
  free((*bp).summ_mufine_fname);
  free((*bp).summ_params_fname);
  free(bp);
}

Bars_PriorParams *New_PriorParams(){
  Bars_PriorParams *pp;
  pp = (Bars_PriorParams *)calloc(1L,sizeof(Bars_PriorParams));
  (*pp).prior_id = UNSPECIFIED;
  (*pp).num_dparams = 0;
  (*pp).num_iparams = 0;
  return pp;
}

void Set_PriorParams(Bars_PriorParams *pp, Bars_BarsParams *bp, enum PriorForm pid, double *dpar, int *ipar){
    /*  Sets the form and parameters for the prior. 
        Is responsible for ensuring that practical, as opposed to theoretical, 
            parameter constraints are satisfied. 
	Is responsible for ensuring constraints on non-prior parameters that
	    depend on prior parameters.
	Note that constraints on non-prior parameters that are *not* 
	    dependent on the prior are enforced in imposeParamConstraints().
            
        This function is generally called twice - once for the default
            parameter values, and again if a parameter file is specified. 
	    Caution is urged before adjusting non-prior parameters
	    based on prior parameter values. Currently,
            the only non-prior parameter subject to change is the
	    initial number of knots, and only in the case that the
	    initial number of knots has prior probability of zero.
    */
  enum PriorForm pid2 = pid;
  int nd,ni;
  nd = (*pp).num_dparams;
  ni = (*pp).num_iparams;
  if (nd > 0) free((*pp).dparams);
  if (ni > 0) free((*pp).iparams);
  (*pp).prior_id = pid2;
  switch(pid2){
  case POISSON: 
    nd = 1;
    ni = 0;
    break;
  case UNIFORM:
    nd = 0;
    ni = 2;
    break;
  case USER:
    nd = MAXKNOTS + 1;
    ni = 2;
    break;
  default:
    pid2 = UNSPECIFIED;
    (*pp).prior_id = UNIFORM;
    nd = 0;
    ni = 2;
    break;
  }
  (*pp).num_dparams = nd;
  (*pp).num_iparams = ni;
  if (pid2 == UNSPECIFIED){
    (*pp).iparams = tssivec(ni);
    (*pp).iparams[0] = MAXKNOTS;
  } else {
    if (nd > 0){
        (*pp).dparams = tssdvec(nd);
        memcpy((*pp).dparams,dpar,sizeof(double)*nd);
    }
    if (ni > 0){
        (*pp).iparams = tssivec(ni);
        memcpy((*pp).iparams,ipar,sizeof(int)*ni);
    }
  }
  /* 
  impose prior parameter constraints
  also ensures that initial number of knots has positive prior probability
  parameter constraints that are not dependent on the prior are enforced in
  imposeParamConstraints().
  */
  if ((*pp).prior_id == POISSON){
    asSafePositive((*pp).dparams);
  }
  if ((*pp).prior_id == UNIFORM){
    if ((*pp).iparams[0] < 1) (*pp).iparams[0] = 1;
    if ((*pp).iparams[0] > MAXKNOTS) (*pp).iparams[0] = MAXKNOTS;
    if ((*bp).k < (*pp).iparams[0]) (*bp).k = (*pp).iparams[0];
    if ((*pp).iparams[1] < (*pp).iparams[0]) (*pp).iparams[1] = (*pp).iparams[0];
    if ((*pp).iparams[1] > MAXKNOTS) (*pp).iparams[1] = MAXKNOTS;
    if ((*bp).k > (*pp).iparams[1]) (*bp).k = (*pp).iparams[1];
  }
  if ((*pp).prior_id == USER){
    if ((*pp).iparams[0] < 1) (*pp).iparams[0] = 1;
    if ((*pp).iparams[0] > MAXKNOTS) (*pp).iparams[0] = MAXKNOTS;
    if ((*bp).k < (*pp).iparams[0]) (*bp).k = (*pp).iparams[0];
    if ((*pp).iparams[1] < (*pp).iparams[0]) (*pp).iparams[1] = (*pp).iparams[0];
    if ((*pp).iparams[1] > MAXKNOTS) (*pp).iparams[1] = MAXKNOTS;
    if ((*bp).k > (*pp).iparams[1]) (*bp).k = (*pp).iparams[1];
  }
   
    
}    


void Free_PriorParams(Bars_PriorParams *pp){
  if ((*pp).num_dparams > 0) free((*pp).dparams);
  if ((*pp).num_iparams > 0) free((*pp).iparams);
  free(pp);
}
  
Bars_BinnedData *New_BinnedData(int nd, int nf){
  Bars_BinnedData *bd;
  bd = (Bars_BinnedData *)calloc(1L,sizeof(Bars_BinnedData));
  (*bd).n = nd;
  (*bd).nf = nf;
  (*bd).x_raw = tssdvec(nd);
  (*bd).xg = New_X_Grid(nd);
  (*bd).fg = New_X_Grid(nf);
  (*bd).y = tssdvec(nd);
  (*bd).x_rawmax = 0.0;
  (*bd).x_rawmin = 0.0;
  (*bd).bin_width = 0.0;
  (*bd).trials = 0;
  (*bd).scale_factor = 0.0;
  return bd;
}


void Free_BinnedData(Bars_BinnedData *bd){
  free((*bd).x_raw);
  Free_X_Grid((*bd).xg);
  Free_X_Grid((*bd).fg);
  free((*bd).y);
  free(bd);
}

Bars_Model *New_Model(int nd, int nf){
  Bars_Model *m;
  m = (Bars_Model *)calloc(1L,sizeof(Bars_Model));
  (*m).n = nd;
  (*m).nf = nf;
  (*m).k = 0;
  (*m).fit_info = 0;
  (*m).fit_fine = 0;
  (*m).accept_prob = 0.0;
  (*m).knots_interior = tssdvec(MAXKNOTS);
  (*m).knots_all = tssdvec(MAXKNOTS + 8);
  (*m).basis = tssdmat2(nd,MAXKNOTS + 4);
  (*m).fitted_values = tssdvec(nd);
  (*m).loglikelihood = 0.0;
  (*m).J = tssdmat2(MAXKNOTS + 4,MAXKNOTS + 4);
  (*m).beta = tssdvec(MAXKNOTS + 4);
  (*m).sigma = 1.0;
  (*m).full_bic = 0.0;
  (*m).full_loglikelihood = 0.0;
  (*m).basis_fine = tssdmat2(nf,MAXKNOTS + 4);
  (*m).random_mu_fine = tssdvec(nf);
  (*m).random_mu = tssdvec(nd);
  return m;
}

void Free_Model(Bars_Model *m){
  free((*m).knots_interior);
  free((*m).knots_all);
  free((*m).basis);
  free((*m).fitted_values);
  free((*m).J);
  free((*m).beta);
  free((*m).basis_fine);
  free((*m).random_mu_fine);
  free((*m).random_mu);
  free(m);
}


Bars_WorkSpace *New_WorkSpace(int n,int nf){
  int p = MAXKNOTS + 4,nm;
  Bars_WorkSpace *ws;
  ws = (Bars_WorkSpace *)calloc(1L,sizeof(Bars_WorkSpace));
  nm = imax(n,nf);
  (*ws).n_max = nm;
  (*ws).conwork = tssdmat2(p,3);
  (*ws).coniwork = tssivec(p);
  (*ws).con_basis = tssdmat2(p,p);
  (*ws).t_con_basis = tssdmat2(p,p);
  (*ws).tmp_basis = tssdmat2(nm,4);
  (*ws).tmp_knots = tssdvec(p);
  (*ws).x_con = tssdvec(2);
  (*ws).work = tssdvec(2 * p * n);
  (*ws).tau = tssdvec(p);
  (*ws).deriv = tssivec(nm);
  (*ws).off = tssivec(nm);
  (*ws).index_m = tssivec(4*nm);
  (*ws).deriv_con = tssivec(2);
  (*ws).x_con[0] = 0.0;
  (*ws).x_con[1] = 1.0;
  (*ws).deriv_con[0] = 2;
  (*ws).deriv_con[1] = 2;
  return ws;
}

void Free_WorkSpace(Bars_WorkSpace *ws){
  free((*ws).conwork);
  free((*ws).coniwork);
  free((*ws).con_basis);
  free((*ws).t_con_basis);
  free((*ws).tmp_basis);
  free((*ws).tmp_knots);
  free((*ws).x_con);
  free((*ws).work);
  free((*ws).tau);
  free((*ws).deriv);
  free((*ws).off);
  free((*ws).index_m);
  free((*ws).deriv_con);
  free(ws);
}

Bars_OutputStat *New_OutputStat(int np){
  Bars_OutputStat *os;
  os = (Bars_OutputStat *)calloc(1L,sizeof(Bars_OutputStat));
  (*os).np = np;
  (*os).params = tssdmat(np,4);
  (*os).accept = tssdvec(4);
/* accept is proportion of accepted births,deaths,relocates, and moves */
  return os;
}

void Free_OutputStat(Bars_OutputStat *os){
  int np,i;
  np = (*os).np;
  for(i=0;i<np;i++){
    free((*os).params[i]);
  }
  free((*os).params);
  free((*os).accept);
  free(os);
}

void AddIntercept(int n, double *x){
  while(n--){
    *x++ = 1.0;
  }
}


void SetModelBasis(Bars_Model *m, 
    Bars_BinnedData *bd, 
    Bars_WorkSpace *ws, int fit_version, int *info){
    
/*  
    Using the knots in Model m, and data contained in BinnedData bd,
    forms the natural spline basis in Model m.
    
    If fit_version == 1, the fine grid is used. 
    Otherwise, the design grid is used.
    
    On output, info == 1 if there has been an error, info == 0 otherwise.
    
    IMPORTANT: The basis is formed such that the first two columns, which
    represent the natural spline constraints, are not to be used.
    
    The basis is formed using code by Bates and Venables, 1998. 
    The code falls under the GNU General Public License. See
    bars_utils.c for complete licensing information.  
*/
  
  int i,temp,j,tt,maxp,numx,lwork;
  Bars_X_Grid *xg;
  double *basis;

  if (fit_version == 1){
    xg = (*bd).fg;
    basis = (*m).basis_fine;
  } else {
    xg = (*bd).xg;
    basis = (*m).basis;
  }

  
  numx = (*xg).n;
  temp = (*m).k + 4;
  maxp = MAXKNOTS + 4;
  tt = numx * temp;
  lwork = 2 * maxp * numx;
 
  for(i=0;i<tt;i++) basis[i] = 0.0;
  tt = 2 * temp;
  for(i=0;i<tt;i++) (*ws).con_basis[i] = 0.0;
  
  spline_basis((*m).knots_all,&temp,&iFour,(*xg).x,
           (*ws).deriv,&numx,(*ws).tmp_basis,(*ws).off);
  outer_plus(iFour,(*ws).off,numx,(*ws).index_m);
  basis_cw_trans((*ws).tmp_basis,basis,(*ws).index_m,iFour,numx);

  spline_basis((*m).knots_all,&temp,&iFour,(*ws).x_con,(*ws).deriv_con,&iTwo,
           (*ws).tmp_basis, (*ws).off);
  outer_plus(iFour, (*ws).off, iTwo, (*ws).index_m);
  basis_cw_trans((*ws).tmp_basis, (*ws).con_basis, (*ws).index_m, iFour, iTwo);

  
  /* The resulting basis will have k + 4 columns, the first
     two of which should not be used. The third column is a 
     column of ones for the intercept.

     To isolate the intercept, the following changes were made:
       a) reduced temp by 1 below,
       b) inserted the "+ 2" in the call to transpose_matrix, and
       c) call AddIntercept to put ones in the third column.
       */

  temp -= 1;

  transpose_matrix(iTwo,temp,(*ws).con_basis + 2,(*ws).t_con_basis);

  
  i = 0;
  *info = 0;
  dgeqrf(&temp,&iTwo,(*ws).t_con_basis,&temp,(*ws).tau,
     (*ws).work,&lwork,&i);
  if (i != 0){
    *info = 1;
  } else {
    dormqr(&Right,&NoTrans,&numx,&temp,&iTwo,(*ws).t_con_basis,
        &temp,(*ws).tau,basis + numx,&numx,(*ws).work,&lwork,&i);
    if (i != 0) *info = 1;
  }

  AddIntercept(numx,basis + 2 * numx);
}



void FitPoissonModel(Bars_Model *m, 
         Bars_BinnedData *bd, 
         Bars_WorkSpace *ws, double *mu_start,double *condnum){

  double *g,*h,*J,*theta,loglik0=0.0,loglik1 = 0.0,delta = 0.0,
  anorm = 1.0, rcond= -1.0, *x, *y, eps_cond;
  int info=88, p,p1,p2,i,iter,converged,n;

  p = (*m).k + 2;
  p1 = p + 1;
  p2 = p * p;
  n = (*bd).n;

  eps_cond = EPS_FACTOR * sqrt(DBL_EPSILON);
  /* 
     The fit is rejected if X'WX has condition
     number (in p=1-norm) c such that
     c > 1 / eps_cond
   */



  
  g = tssdmat2(n,p1);
  h = tssdmat2(n,p);
  J = tssdmat2(p,p1);
  theta = tssdmat2(n,1);

  x = (*m).basis + 2*n;
  y = (*bd).y;
/* IMPORTANT: the first two basis columns are NOT to be used. */


  memcpy(g,x,(sizeof(double)*n*p));
  memcpy((*m).fitted_values,mu_start,(sizeof(double)*n));
  i = MAXIT_GLM;
  iter = 0;
  converged = 1;
  while(i--){
    /* calculate Z = adjusted dependent variable. Values are stored
       in the last column of g */
    glmp_getz(n,(g + (n*p)),y, (*m).fitted_values);
    /* calculate h = WX, where W = Diag(mu) */
    glmp_geth(n,p,x,(*m).fitted_values,h);
    /* calculate J = h'g = [X'WX|X'WZ] */
    dgemm(&Trans,&NoTrans,&p,&p1,&n,&One,h,&n,g,&n,&Zero,J,&p);
    anorm = getNorm(J,p);
    /* do Cholesky decomposition of X'WX */
    /* X'WX = U'U. The U matrix is stored in the upper
       triangle of J. The lower triangle remains unchanged,
       that is, it is the lower triangle of X'WX */
    dpotrf(&Upper,&p,J,&p,&info);
    if (info != 0) {
      i=0;
      converged = 0;
    } else {
      /* check condition number of X'WX */
      dpocon(&Upper,&p,J,&p,&anorm,&rcond,(*ws).conwork,(*ws).coniwork,&info);
      if ((info != 0) || (rcond < eps_cond)) {
        i=0;
        converged = 0;
      } else {
        iter++;
	/* solve (X'WX)B = X'WZ for B. Solution is in the last column
	   of J */
        dpotrs(&Upper,&p,&iOne,J,&p,(J + p2),&p,&info);
        if (info != 0){
	  i=0;
	  converged = 0;
        } else {
	  /* calculate eta = XB */
	  dgemm(&NoTrans,&NoTrans,&n,&iOne,&p,&One,x,&n,
                (J + p2),&p,&Zero,theta,&n);
	  /* set mu = exp(eta) */
	  glmp_getmu(n,theta,(*m).fitted_values);
	  /* store Cholesky decomposition of X'WX for
	     later use when we want Cov(beta hat) */
	  memcpy((*m).J,J,(sizeof(double)*p*p1));
	  
	  loglik0 = loglik1;
	  loglik1 = glmp_loglik(n,theta,(*m).fitted_values,y);
	  delta = loglik1 - loglik0;
	  delta = fabs(delta);
	  if (delta < 0.01) i = 0;
        }
      }
    }
  }

  *condnum = rcond;

  (*m).fit_info = iter * converged;
  (*m).loglikelihood = loglik1;
  if (converged){
    (*m).sigma = 1.0;
/*    (*m).sigma = glmp_getsigma(n,p,(*bd).y,(*m).fitted_values); */
/* computes sigma for over-dispersion */

  }

  free(g);
  free(h);
  free(J);
  free(theta);
}

void SetBasisAndFitModel(Bars_Model *m, 
         Bars_BinnedData *bd, 
         Bars_WorkSpace *ws, double *mu_start, double *condnum){
    int info;
  SetModelBasis(m,bd,ws,0,&info);
  if (info == 0){
    FitPoissonModel(m,bd,ws,mu_start,condnum);
  } else {
    *condnum = -1.0;
    (*m).fit_info = 0;
  }
  (*m).fit_fine = 0;
}

void SetFineBasis(Bars_Model *m,
          Bars_BinnedData *bd,
          Bars_WorkSpace *ws){
    int info;
    if ((*m).fit_fine == 0){
        SetModelBasis(m,bd,ws,1,&info);
        (*m).fit_fine = 1;
        if (info != 0) (*m).fit_info = 0;
    }
}

double exp_prior2(Bars_Model *m,double *beta,int p){
  /* flat prior for intercept, unit information for the 
     other coefficients */
  double z,z1,*J;
  int i,i1,n;
  n = (*m).n;
  J = (*m).J;
  z = 0.0;
  for(i=1;i<p;i++){
    z1 = 0.0;
    for(i1=i;i1<p;i1++){
      z1 += beta[i1] * J[(p*i1 + i)];
    }
    z += pow(z1,2);
  }
  z /= ((double) n);
  z *= -0.5;
  return z;
}
  
void RandBeta(Bars_Model *m,Bars_BinnedData *bd, Bars_BarsParams *bp,
	      int *laterej, int iternum){
  int info,p,iters=0,MHiters,MHi,n,countMH;
  double *z,*J,zn,*j2,*curbeta,*lastbeta,*tempbeta,curhi;
  double lasthi,curgi,curfi,lastgi,lastfi,r,u;
  double *beta_one,*beta_two,mlefgh;
  double threshhold;

  MHiters = (*bp).beta_iter;
  threshhold = (*bp).threshhold;
  z = (*m).beta;
  J = (*m).J;
  p = (*m).k + 2;
  n = (*m).n;
  
  beta_one = tssdvec(p);
  beta_two = tssdvec(p);
  curbeta = beta_one;
  lastbeta = beta_two;
   
  memcpy(lastbeta,J + p*p, sizeof(double)*p);

  exb(n,p,(*m).basis + 2 * n,lastbeta,
      (*m).random_mu);
  lastfi = glmp_loglik_3(n,(*m).random_mu,(*bd).y);
  lastgi = exp_prior2(m,lastbeta,p);
  lasthi = 0.0;
  mlefgh = lastfi + lastgi;
  countMH = 0;
  for(MHi=0;MHi<MHiters;MHi++){
    iters = 0;
    info = 0;
    do{
      iters++;
      rmvnorm(p,curbeta);
      curhi = squared_norm(curbeta,p);

      /* solve (X'WX)^(1/2) A = z for A. 
	 So A = (X'WX)^(-1/2)z ~ N(0,(X'WX)^(-1)) 
	 Solution is placed in z */
      dtrtrs(&Upper,&NoTrans,&NoUnit,&p,&iOne,J,&p,curbeta,&p,&info);
    } while((info != 0) && (iters < 20));
    if ((iters == 20) && (info != 0)) {
      mexErrMsgTxt("Maximum number of beta iterations reached.\n");
    }
    /* z = z + beta_hat ~ N(beta_hat,(X'WX)^(-1)) */
    /* sigma = 1 (no overdispersion) */
    curhi *= -0.5;
    affine(p,1.0,curbeta,(J + p*p));
    exb(n,p,(*m).basis + 2 * n,curbeta,
	(*m).random_mu);
    curfi = glmp_loglik_3(n,(*m).random_mu,(*bd).y);
    curgi = exp_prior2(m,curbeta,p);
    r = (curfi - lastfi) + (curgi - lastgi) - (curhi - lasthi);

    if (r > 0.0) r = 0.0;
    u = (double) genunf(0.0,1.0);
    u = (u < DBL_EPSILON) ? (r - 1.0) : log(u);
    if ((MHi == 0) && (r > threshhold)){
      /* automatic acceptance */
      MHi = MHiters;
      u = r - 1.0;
    }
    if (u < r) {
      /* accept MH candidate */
      tempbeta = curbeta;
      curbeta = lastbeta;
      lastbeta = tempbeta;
      lastfi = curfi;
      lastgi = curgi;
      lasthi = curhi;
      countMH++;
    }

  }

  memcpy(z,lastbeta,sizeof(double) * p);

  free(beta_one);
  free(beta_two);
  *laterej = (countMH == 0) ? 1 : 0;
}




void RandMu(Bars_Model *m, Bars_BinnedData *bd, Bars_BarsParams *bp, 
	    int *laterej,int iternum){
  /* Generate random beta, and calculate eta = exp(X*beta)
     for both the design basis and the grid basis */
  RandBeta(m,bd,bp,laterej,iternum);
  exb((*m).n,(*m).k + 2,(*m).basis + 2 * (*m).n,(*m).beta,
      (*m).random_mu);
  exb((*m).nf,(*m).k + 2,(*m).basis_fine + 2 * (*m).nf,(*m).beta,
      (*m).random_mu_fine);
}

void exp_prior(Bars_Model *m){
  double z = 0.0;
  int i=0,n;
  n = (*m).n;
  for(i=0;i<n;i++){
    z += pow(log((*m).random_mu[i]),2) * (*m).fitted_values[i];
  }
  z /= ((double) n);
}

void Full_BIC(Bars_Model *m, Bars_BinnedData *bd){
/* 
   BIC for the FULL model (\xi,k,\beta), as opposed to the
   marginal model (\xi,k). k is always one dimensional, so the
   parameter count is 1 for k, k for \xi plus (k+2) for \beta.
   BIC = l - .5(2*k + 3)log(n)
       = l - (k + 1.5)log(n).
*/
  double v = 0.0, w = 0.0,*tt;
  int i,j,p;
  w = (double)((*m).k) + 1.5;
  v = glmp_loglik_3((*m).n,(*m).random_mu,(*bd).y);
  (*m).full_bic = v - (w * log((*m).n));
  (*m).full_loglikelihood = v;
}

void SetEquallySpacedKnots(Bars_Model *m){
  set_default_knots((*m).k,(*m).knots_interior);
}
void SetLogsplineKnots(Bars_Model *m,Bars_BinnedData *bd, Bars_PriorParams *pp){
  Bars_X_Grid *xg;
  int k;
  xg = (*bd).xg;
  get_logspline_knots((*xg).n,(*xg).x,(*bd).y,(*m).knots_interior,&k,(*bd).bin_width/((*bd).x_rawmax - (*bd).x_rawmin));
  if (((*pp).prior_id == UNIFORM) || ((*pp).prior_id == USER)){
    /* enforce constraint L <= k <= U*/
    if ((k < (*pp).iparams[0]) || (k > (*pp).iparams[1])){
      k = (int) (0.5 * ((*pp).iparams[0] + (*pp).iparams[1]));
      printf("The number of logspline knots has zero prior probability. Using %i equally spaced knots to begin chain.\n\n",k);
      (*m).k = k;
      SetEquallySpacedKnots(m);
    } else {
      (*m).k = k;
    }
  } else {
    /* enforce constraint k <= MAXKNOTS */
    (*m).k = imin(k,MAXKNOTS);
  }
}

void SetKnotsAll(Bars_Model *m){
  boundaries((*m).k,iFour,(*m).knots_interior,(*m).knots_all);
}

void FormFirstKnots(Bars_Model *m, Bars_BarsParams *bp, Bars_BinnedData *bd,Bars_PriorParams *pp){
  if ((*bp).use_logspline){
    SetLogsplineKnots(m,bd,pp);
  } else {
    (*m).k = (*bp).k;
    SetEquallySpacedKnots(m);
  }
  SetKnotsAll(m);
}



void SaveMeanAndMode(Bars_SampStat *ss,double *dout, int pnum){
  double *v,*w;
  v = (*ss).means;
  w = (*ss).modes;
  *dout++ = v[pnum];
  *dout = w[pnum];
}


void SaveOutputStat(Bars_OutputStat *os, Bars_SampStat *ss,int n, int sim, int *countup,int *moves, double conf_level){
  int i,np = (*os).np;
  for(i=0;i<np;i++){
    conf_int((*ss).params[i],sim,conf_level,(*os).params[i]);
    SaveMeanAndMode(ss,(*os).params[i] + 2,i);
  }
  for(i=0;i<4;i++){
    (*os).accept[i] = ((double)(countup[i]))/((double)(moves[i]));
  }
}


Bars_SampStat *New_SampStat (int n, int nf, int np, int sim){
  Bars_SampStat *ss;
  ss = (Bars_SampStat *)calloc(1L,sizeof(Bars_SampStat));
  (*ss).np = np;
  (*ss).sim = sim;
  (*ss).random_mu_mean = tssdvec0(n);
  (*ss).random_mufine_mean = tssdvec0(nf);
  (*ss).random_mu_mode = tssdvec0(n);
  (*ss).random_mufine_mode = tssdvec0(nf);
  (*ss).params = tssdmat(np,sim);
  (*ss).full_bic = tssdvec0(sim);
  (*ss).means = tssdvec0(np);
  (*ss).modes = tssdvec0(np);
  (*ss).max_full_bic = 0.0;
  return ss;
}

void Free_SampStat(Bars_SampStat *ss){
  int i,np;
  free((*ss).random_mu_mean);
  free((*ss).random_mufine_mean);
  free((*ss).random_mu_mode);
  free((*ss).random_mufine_mode);
  np = (*ss).np;
  for(i=0;i<np;i++){
    free((*ss).params[i]);
  }
  free((*ss).params);
  free((*ss).full_bic);
  free((*ss).means);
  free((*ss).modes);
  free(ss);
}

void setSampStatParam(double val, Bars_SampStat *ss, int pnum, int simnum){
  (*ss).params[pnum][simnum] = val;
}

void setSampStatParams(int vlen, double *vals, Bars_SampStat *ss, 
               int pnum, int simnum){
  int i;
  for(i=0;i<vlen;i++){
    (*ss).params[(pnum + i)][simnum] = vals[i];
  }
}

void saveModelStats(Bars_Model *m, Bars_BinnedData *bd, 
          Bars_SampStat *ss, int simnum, double *z){
  
  int nf;
  double *uf,*df,*wf,*pf;
  nf = (*m).nf;
  uf = tssdvec0(nf);
  df = tssdvec0(nf);
  wf = tssdvec0(nf);
  pf = tssdvec0(nf);
/*
  get_mode_stats(nf,(*((*bd).fg)).x,(*m).random_mu_fine,(*bd).x_rawmin, 
         (*bd).x_rawmax, z);
*/
  get_mode_stats0(nf,(*((*bd).fg)).x,(*m).random_mu_fine,z);
  max_spline(nf,(*((*bd).fg)).x,(*m).random_mu_fine,z,uf,df,wf,pf);
  z[0] *= (*bd).x_rawmax - (*bd).x_rawmin;
  z[0] += (*bd).x_rawmin;
  z[1] *= (*bd).scale_factor;
  z[2] = (double)((*m).k);

  setSampStatParams((*ss).np,z,ss,0,simnum);
  free(uf);
  free(df);
  free(wf);
  free(pf);
}



  
void saveBIC(Bars_Model *m, Bars_BinnedData *bd, Bars_SampStat *ss, int simnum){

  Full_BIC(m,bd);
  (*ss).full_bic[simnum] = (*m).full_bic;
}  


void saveSampStatParams(Bars_Model *m, Bars_BinnedData *bd,
            Bars_SampStat *ss, int simnum, double *z){
  int n,nf;
  saveModelStats(m,bd,ss,simnum,z);
  saveBIC(m,bd,ss,simnum);
  n = (*m).n;
  nf = (*m).nf;
  cumsum(n,(*ss).random_mu_mean,(*m).random_mu);
  cumsum(nf,(*ss).random_mufine_mean,(*m).random_mu_fine);
  cumsum((*ss).np,(*ss).means,z);
  if (((*m).full_bic > (*ss).max_full_bic) || (simnum == 0)){
    (*ss).max_full_bic = (*m).full_bic;
    memcpy((*ss).modes,z,sizeof(double) * (*ss).np);
    memcpy((*ss).random_mu_mode,(*m).random_mu,sizeof(double) * n);
    memcpy((*ss).random_mufine_mode,(*m).random_mu_fine,sizeof(double) * nf);
  }
}


void SetMuStart(int n, double *mu, Bars_BinnedData *bd){
  double *y;

  y = (*bd).y;
  while(n--){
    *mu++ = ((*y) < 0.1) ? 0.1 : (*y);
    y++;
  }
}
void CheckFirstFit(Bars_Model *m, Bars_BinnedData *bd,
		   Bars_PriorParams *pp,
           Bars_WorkSpace *ws, double *mu_start){
  int k,i,remk,cnk,kmin;
  double *knot_set, *knot_sub, like,condnum,maxcn;


  if ((*m).fit_info == 0){
    /* initial knots are bad */
    k = (*m).k;
    knot_set = tssdvec(k);
    knot_sub = (*m).knots_interior;
    if (((*pp).prior_id == UNIFORM) || ((*pp).prior_id == USER)){
      kmin = (*pp).iparams[0];
    } else kmin = 1;
    while (((*m).fit_info == 0) && (k > kmin)){
      /* the current set of k knots are bad
	 backwards elimination procedure
	 forms subsets by removing one knot. 
	 Among the subsets that do not induce
	 ill-conditioned matrices, the routine
	 chooses the one with the greatest
	 log likelihood. If all have ill-conditioned
	 matrices, then it chooses the subset with
	 the smallest condition number. In the event
	 that condition numbers can not even be computed
	 for all of the subsets, then the routine removes
	 the first knot in the set.
	 */
      memcpy(knot_set,knot_sub,sizeof(double) * k);
      maxcn = -1.0;
      cnk = -1;
      remk = -1;
      (*m).k -= 1;
      k -= 1;
      for(i=0;i<=k;i++){
	if (i==0) memcpy(knot_sub,knot_set + 1,sizeof(double)*k);
	else if (i==k) memcpy(knot_sub,knot_set,sizeof(double)*k);
	else {
	  memcpy(knot_sub,knot_set,sizeof(double)*i);
	  memcpy(knot_sub + i,knot_set + (i+1),sizeof(double)*(k-i));
	}
	SetKnotsAll(m);
	SetBasisAndFitModel(m,bd,ws,mu_start,&condnum);
	if (condnum > maxcn) {
	  cnk = i;
	  maxcn = condnum;
	}
	if ((*m).fit_info != 0){
	  /* this model is okay */
	  if (remk == -1) {
	    remk = i;
	    like = (*m).loglikelihood;
	  } else {
	    if ((*m).loglikelihood > like){
	      remk = i;
	      like = (*m).loglikelihood;
	    }
	  }
	}
      }
      if (remk == -1) {
	if (cnk == -1){
	  remk = 0;
	} else {
	  remk = cnk;
	}
      }
      /* remove knot remk from knot set and re-fit. */
      if (remk==0) memcpy(knot_sub,knot_set + 1,sizeof(double)*k);
      else if (remk==k) memcpy(knot_sub,knot_set,sizeof(double)*k);
      else {
	memcpy(knot_sub,knot_set,sizeof(double)*remk);
	memcpy(knot_sub + remk,knot_set + (remk+1),sizeof(double)*(k-remk));
      }
      SetKnotsAll(m);
      SetBasisAndFitModel(m,bd,ws,mu_start,&condnum);
    }
    if ((*m).fit_info == 0){
      /* backwards elimination procedure failed to find
	 ANY feasible models. Last resort: try minimum
	 number of possible knots, equally spaced.
	 */
      (*m).k = kmin;
      SetEquallySpacedKnots(m);
      SetKnotsAll(m);
      SetBasisAndFitModel(m,bd,ws,mu_start,&condnum);
      if ((*m).fit_info == 0){
		/* still bad - give up */
		mexErrMsgTxt("Irreparable initial knots.\nThis may have occurred if the prior on the number of knots has too small of a support.\n");
      }
    }
  }
}





void SetPriorRatios(double *pratio, Bars_PriorParams *pp){
  /* pratio[k] = pi(k+1)/pi(k) for k>0. pratio[0] = 0 */
  enum PriorForm prior_id;
  int i,capel,capu;
  double lambda;
  pratio[0] = 0.0;
  pratio[MAXKNOTS] = 0.0;
  prior_id = (*pp).prior_id;
  switch(prior_id){
  case POISSON:
    /* Poisson(lambda) truncated to {1,..,MAXKNOTS} */
    lambda = (*pp).dparams[0];
    for(i=1;i<MAXKNOTS;i++){
      pratio[i] = lambda / ((double)(i+1));
    }
    break;
  case UNIFORM:
    /* uniform(L,..,U) where U <= MAXKNOTS and L >= 1 */
    capel = (*pp).iparams[0];
    capu = (*pp).iparams[1];
    for(i=1;i<MAXKNOTS;i++){
        pratio[i] = ((i >= capel) && (i < capu)) ? 1.0 : 0.0;
    }
    break;
  case USER:
    /* user defined prior */
    capel = (*pp).iparams[0];
    capu = (*pp).iparams[1];
    for(i=1;i<MAXKNOTS;i++){
      pratio[i] = 0.0;
    }
    for(i=capel;i<capu;i++){
      pratio[i] = (*pp).dparams[i+1] / (*pp).dparams[i];
    }
    break;
  default:
    for(i=1;i<MAXKNOTS;i++){
      pratio[i] = 0.0;
    }
    break;
  }
}

void SetRJProbs(double *birth, double *death, Bars_BarsParams *bp,
        Bars_PriorParams *pp){
  int i;
  double c, *pratio;
  /* pratio[k] = pi(k+1)/pi(k) for k>0. pratio[0] = 0 */

  pratio = tssdvec(MAXKNOTS + 1);
  SetPriorRatios(pratio,pp);
  c = (*bp).probbd;
  birth[0] = 0.0;
  death[0] = 0.0;
  birth[1] = (pratio[1] < DBL_EPSILON) ? 0.0 : c * min(1.0,pratio[1]);
  death[1] = 0.0;
  birth[MAXKNOTS] = 0.0;
  if (pratio[MAXKNOTS - 1] < DBL_EPSILON){
    death[MAXKNOTS] = c;
  } else {
    death[MAXKNOTS] = c * min(1.0,(1.0/pratio[MAXKNOTS - 1]));
  }
  for(i=2;i<MAXKNOTS;i++){
    birth[i] = (pratio[i] < DBL_EPSILON) ? 0.0 : c * min(1.0,pratio[i]);
    death[i] = (pratio[i-1] < DBL_EPSILON) ? c : 
      c * min(1.0,(1.0/pratio[i-1]));
  }
  free(pratio);
}

int GetDeathKnotPos(int k){
  return((int) ignuin(0L, (long)(k-1)));
}

void GetBirthKnotInfo(int *k_rand,double *cand, 
               double tau, int k, double *knots){
  double alpha,beta;
  *k_rand = GetDeathKnotPos(k);
  alpha = tau * knots[(*k_rand)];
  beta = tau - alpha;
  do {
    *cand = (double)genbet(alpha,beta);
  } while (check_cand(*cand,k,knots));
}

double GetTransitionDensity(double cand, double tau, int k, double *knots){
  int i;
  double dens = 0.0;
  for(i=0;i<k;i++){
    dens += dbeta(cand,tau*knots[i],tau * (1.0 - knots[i]));
  }
  return dens;
}


void SetBirthModelKnots(double cand, int k, double *knots, Bars_Model *m2){
  double *knots2;
  int newknotpos,nkp1;

  (*m2).k = k + 1;
  knots2 = (*m2).knots_interior;

  newknotpos = binary_search(cand,knots,k);
  if (newknotpos==0 && knots[0]>cand){
    knots2[0] = cand;
    memcpy(knots2 + 1, knots, sizeof(double)*k);
  } else if (newknotpos==k-1){
    memcpy(knots2, knots,sizeof(double)*k);
    knots2[k] = cand;
  } else {
    nkp1 = newknotpos+1;
    memcpy(knots2, knots, sizeof(double)*nkp1);
    knots2[nkp1] = cand;
    memcpy(knots2 + nkp1 + 1,knots + nkp1,
       sizeof(double)*(k - nkp1));
  }
  boundaries(k+1,4,knots2,(*m2).knots_all);
}

void SetDeathModelKnots(int dk_pos, int k, double *knots, Bars_Model *m2){
  double *knots2;

  knots2 = (*m2).knots_interior;
  (*m2).k = k - 1;
  if (dk_pos > 0)
    memcpy(knots2,knots,sizeof(double)*dk_pos);
  if (dk_pos < k - 1)
    memcpy(knots2 + dk_pos,knots + dk_pos + 1,sizeof(double)*(k - dk_pos - 1));
  boundaries(k-1,4,knots2,(*m2).knots_all);
}

void SetRelocationModelKnots(int dk_pos, double birth_cand, int k,
                 double *knots, Bars_Model *m2, double *knots3){
  double *knots2;

  knots2 = (*m2).knots_interior;
  (*m2).k = k;
  if(k==1) {
    knots2[0] = birth_cand;
    boundaries(k,4,knots2,(*m2).knots_all);
  } else {
    SetDeathModelKnots(dk_pos,k,knots,m2);
    memcpy(knots3,knots2, sizeof(double)*(k-1));
    SetBirthModelKnots(birth_cand,k-1,knots3,m2);
  }
}

  
double protected_exp(double x){
  if (x > 80.0) x = 80.0;
  else if (x < -80.0) x = -80.0;
  return exp(x);
}

void SetBirthModel(Bars_Model *m1, Bars_Model *m2, Bars_BinnedData *bd, 
         Bars_BarsParams *bp, Bars_WorkSpace *ws, double *mu_start){
  int k,k_rand;
  double tau,*knots,cand,dens,condnum = -1.0;

  k = (*m1).k;
  knots = (*m1).knots_interior;
  tau = (*bp).tau;
  GetBirthKnotInfo(&k_rand,&cand,tau,k,knots);
  dens = GetTransitionDensity(cand,tau,k,knots);
  SetBirthModelKnots(cand,k,knots,m2);
  SetBasisAndFitModel(m2,bd,ws,mu_start,&condnum);
  (*m2).accept_prob = ((*m2).fit_info > 0) ?
    protected_exp(((*m2).loglikelihood - (*m1).loglikelihood + 
           log(k) - log(dens) - 0.5*log((*m1).n))) : 0.0;
}
void SetDeathModel(Bars_Model *m1, Bars_Model *m2, Bars_BinnedData *bd, 
         Bars_BarsParams *bp, Bars_WorkSpace *ws, double *mu_start){
  int k,dk_pos;
  double tau,*knots,cand,dens,condnum = -1.0;

  k = (*m1).k;
  knots = (*m1).knots_interior;
  tau = (*bp).tau;
  dk_pos = GetDeathKnotPos(k);
  cand = knots[dk_pos];
  SetDeathModelKnots(dk_pos,k,knots,m2);
  dens = GetTransitionDensity(cand,tau,k-1,(*m2).knots_interior);
  SetBasisAndFitModel(m2,bd,ws,mu_start,&condnum);
  (*m2).accept_prob = ((*m2).fit_info > 0) ? 
    protected_exp(((*m2).loglikelihood - (*m1).loglikelihood - 
               log(k-1) + log(dens) + 0.5*log((*m1).n))) : 0.0;
}
void SetRelocateModel(Bars_Model *m1, Bars_Model *m2, 
              Bars_BinnedData *bd, Bars_BarsParams *bp, 
              Bars_WorkSpace *ws, double *mu_start){
  int k,dk_pos;
  double tau,*knots,birth_cand,death_cand,dens1,dens2,condnum = -1.0;

  k = (*m1).k;
  knots = (*m1).knots_interior;
  tau = (*bp).tau;
  GetBirthKnotInfo(&dk_pos,&birth_cand,tau,k,knots);
  death_cand = knots[dk_pos];
  SetRelocationModelKnots(dk_pos,birth_cand,k,knots,m2,(*ws).tmp_knots);
  dens1 = dbeta(birth_cand,tau*death_cand,tau*(1.0 - death_cand));
  dens2 = dbeta(death_cand,tau*birth_cand,tau*(1.0 - birth_cand));
  SetBasisAndFitModel(m2,bd,ws,mu_start,&condnum);
  (*m2).accept_prob = ((*m2).fit_info > 0) ?
    protected_exp(((*m2).loglikelihood - (*m1).loglikelihood +
           log(dens2) - log(dens1))) : 0.0;

}

void FinishSampStats(Bars_SampStat *ss, int samp_iter, int n, int nf){
  double *mu,*wt,*means,div_by,z;
  int np;

  mu = (*ss).random_mu_mean;
  div_by = (double)samp_iter;
  while(n--) *mu++ /= div_by;
  mu = (*ss).random_mufine_mean;
  while(nf--) *mu++ /= div_by;
  np = (*ss).np;
  means = (*ss).means;
  while (np--) *means++ /= div_by;
}

void write_iter_params(Bars_BarsParams *bp, Bars_Model *m, Bars_BinnedData *bd, int np, 
            double *z, int iter_num){
  FILE *the_file;
  double kpos,a,b,*ki;
  int i,k;
  if ((*bp).use_iter_knots){
    the_file = (*bp).iter_knots_file;
    fprintf(the_file,"%i %i ",iter_num,(*m).k);
    a = (*bd).x_rawmax - (*bd).x_rawmin;
    b = (*bd).x_rawmin;
    k = (*m).k;
    ki = (*m).knots_interior;
    for(i=0;i<k;i++){
      kpos = a * ki[i] + b;
      fprintf(the_file,"%le ",kpos);
    }
    fprintf(the_file,"\n");
  }
  if ((*bp).use_iter_params){
    the_file = (*bp).iter_params_file;
    fprintf(the_file,"%i %lf %lf ",iter_num, 
        (*m).full_bic, (*m).full_loglikelihood);
    write_vector(the_file,z,np);
  }
  if ((*bp).use_iter_mu){
    write_vector((*bp).iter_mu_file,(*m).random_mu,(*m).n);
  }
  if ((*bp).use_iter_mufine){
    write_vector((*bp).iter_mufine_file,(*m).random_mu_fine,(*m).nf);
  } 
}

void write_summ_params (Bars_BarsParams *bp, Bars_BinnedData *bd,
            Bars_SampStat *ss,
            Bars_OutputStat *os, int n, int nf){
  double a,b,*x1,*x2;
  int i,np;
  Bars_X_Grid *xg;
  if ((*bp).use_summ_mu){
    write_vector((*bp).summ_mu_file,(*bd).x_raw,n);
    write_vector((*bp).summ_mu_file,(*ss).random_mu_mean,n);
    write_vector((*bp).summ_mu_file,(*ss).random_mu_mode,n);
  }
  if ((*bp).use_summ_mufine){
    x2 = tssdvec(nf);
    xg = (*bd).fg;
    x1 = (*xg).x;
    a = (*bd).x_rawmax - (*bd).x_rawmin;
    b = (*bd).x_rawmin;
    for(i=0;i<nf;i++){
      x2[i] = a * x1[i] + b;
    }
    write_vector((*bp).summ_mufine_file,x2,nf);
    write_vector((*bp).summ_mufine_file,(*ss).random_mufine_mean,nf);
    write_vector((*bp).summ_mufine_file,(*ss).random_mufine_mode,nf);
    free(x2);
  }
  if ((*bp).use_summ_params){
    np = (*os).np;
    for(i=0;i<np;i++){
      write_vector((*bp).summ_params_file,(*os).params[i],4);
    }
  }
  
}

void OpenFiles (Bars_BarsParams *bp, int iters, int summs){
  if (iters){
    if ((*bp).use_iter_knots)
      (*bp).iter_knots_file = fopen((*bp).iter_knots_fname,"w");
    if ((*bp).use_iter_mu) 
      (*bp).iter_mu_file = fopen((*bp).iter_mu_fname,"w");
    if ((*bp).use_iter_mufine) 
      (*bp).iter_mufine_file = fopen((*bp).iter_mufine_fname,"w");
    if ((*bp).use_iter_params) 
      (*bp).iter_params_file = fopen((*bp).iter_params_fname,"w");
  }
  if (summs){
    if ((*bp).use_summ_mu) 
      (*bp).summ_mu_file = fopen((*bp).summ_mu_fname,"w");
    if ((*bp).use_summ_mufine) 
      (*bp).summ_mufine_file = fopen((*bp).summ_mufine_fname,"w");
    if ((*bp).use_summ_params) 
      (*bp).summ_params_file = fopen((*bp).summ_params_fname,"w");
  }
}

void CloseFiles (Bars_BarsParams *bp, int iters, int summs){
  if (iters){
    if ((*bp).use_iter_knots) fclose((*bp).iter_knots_file);
    if ((*bp).use_iter_mu) fclose((*bp).iter_mu_file);
    if ((*bp).use_iter_mufine) fclose((*bp).iter_mufine_file);
    if ((*bp).use_iter_params) fclose((*bp).iter_params_file);
  }
  if (summs){
    if ((*bp).use_summ_mu) fclose((*bp).summ_mu_file);
    if ((*bp).use_summ_mufine) fclose((*bp).summ_mufine_file);
    if ((*bp).use_summ_params) fclose((*bp).summ_params_file);
  }
}

void Bars_MCMC(Bars_Model *m1, Bars_Model *m2, 
	       Bars_BinnedData *bd, 
	       Bars_WorkSpace *ws, 
	       Bars_BarsParams *bp, 
	       Bars_PriorParams *pp,
	       Bars_SampStat *ss,
	       Bars_OutputStat *os){
  
  Bars_Model *m_old, *m_new, *m_temp;
  double *mu_start,*birth_probs,*death_probs,u,*z;
  int i,n,nf,samp_iter,burn_iter,tot_iter,iter_info,k,np;
  int vb=0,*countup,*moves,bdr,late_rejection;
  
  OpenFiles(bp,1,0);

  n = (*bd).n;
  iter_info = 0;
  countup = tssivec(4);
  moves = tssivec(4);
  for(i=0;i<4;i++){
    countup[i] = 0;
    moves[i] = 0;
  }
  samp_iter = (*bp).samp_iter;
  burn_iter = (*bp).burn_iter;
  tot_iter = samp_iter + burn_iter;
  nf = (*bd).nf;

  m_old = m1;
  m_new = m2;


  /* set starting knots using logspline or evenly spaced
     knots, depending on a user-defined parameter.
     Default is to use logspline */
  FormFirstKnots(m_old,bp,bd,pp);

  birth_probs = tssdvec(MAXKNOTS + 1);
  death_probs = tssdvec(MAXKNOTS + 1);
  SetRJProbs(birth_probs,death_probs, bp,pp);

  /* calcuate the initial value of mu to be used
     for every glm fit. It is possible that one can
     do better by using, for example, the solution from
     the previous glm as the starting value. We encountered
     numerical difficulties with this approach. */
  mu_start = tssdvec(n);
  SetMuStart(n,mu_start,bd);
  np = (*ss).np;
  z = tssdvec(np);

  u = -1.0;

  SetBasisAndFitModel(m_old,bd,ws,mu_start,&u);
  CheckFirstFit(m_old,bd,pp,ws,mu_start);

  for(i=0;i<tot_iter;i++){
/*    printf("[[%i]] ",i);  */
    k = (*m_old).k;
    u = (double) genunf(0.0,1.0);
    if (u <= birth_probs[k]){
      /* Birth Step */
      bdr = 0;
      SetBirthModel(m_old,m_new,bd,bp,ws,mu_start);
    } else if ((1.0 - u) < death_probs[k]){
      /* Death Step */
      bdr = 1;
      SetDeathModel(m_old,m_new,bd,bp,ws,mu_start);
    } else {
      /* Relocation Step */
      bdr = 2;
      SetRelocateModel(m_old,m_new,bd,bp,ws,mu_start);
    }
    moves[bdr]++;
    moves[3]++;
    u = (double) genunf(0.0,1.0);
/*    u = (*m_new).accept_prob + 1.0; */
    if (u < (*m_new).accept_prob){
      /* New knot set accepted. Swap Models. */
      m_temp = m_old;
      m_old = m_new;
      m_new = m_temp;
      countup[bdr]++;
      countup[3]++;
    }
    if (i >= burn_iter){
      SetFineBasis(m_old,bd,ws);
      late_rejection = 0;
      RandMu(m_old,bd,bp,&late_rejection,i);
      if (late_rejection){
	i--;
      } else {
	saveSampStatParams(m_old,bd,ss,i-burn_iter,z);
	write_iter_params(bp,m_old,bd,np,z,i);
      }
    }
  }

  FinishSampStats(ss,samp_iter,n,nf);
  SaveOutputStat(os, ss, n,samp_iter,countup,moves,(*bp).conf_level);

  CloseFiles(bp,1,0);
  OpenFiles(bp,0,1);
  write_summ_params(bp,bd,ss,os,n,nf);
  CloseFiles(bp,0,1);
  
  free(z);
  free(mu_start);
}

void imposeParamConstraints(Bars_BarsParams *bp, Bars_PriorParams *pp){
    if ((*bp).k < 1) (*bp).k = 1;
    else if ((*bp).k > MAXKNOTS) (*bp).k = MAXKNOTS;
    if ((*bp).probbd < 0.0001) (*bp).probbd = 0.0001;
    else if ((*bp).probbd > 0.4999) (*bp).probbd = 0.4999;
    if ((*bp).burn_iter < 0) (*bp).burn_iter = 0;
    if ((*bp).samp_iter < 1) (*bp).samp_iter = 1;
    if ((*bp).tau < 0.0001) (*bp).tau = 0.0001;
    if ((*bp).nf < 1) (*bp).nf = 1;
    if ((*bp).beta_iter < 1) (*bp).beta_iter = 1;
    if ((*bp).conf_level < 0.000001) (*bp).conf_level = 0.000001;
    if ((*bp).conf_level > 0.999999) (*bp).conf_level = 0.999999;
}

void setDefaultParamValues(Bars_BarsParams *bp, Bars_PriorParams *pp){
  double *dpar;
  int *ipar;
  enum PriorForm pform;
  (*bp).burn_iter = 0;
  (*bp).samp_iter = 2000;
  (*bp).k = 3;
  (*bp).probbd = 0.40;
  (*bp).tau = 50.0;
  (*bp).conf_level = 0.95;
  (*bp).nf = N_FINE;
  (*bp).beta_iter = 3;
  (*bp).threshhold = -10.0;
  (*bp).use_iter_knots = 0;
  (*bp).use_iter_mu = 1;
  (*bp).use_iter_mufine = 0;
  (*bp).use_iter_params = 1;
  (*bp).use_summ_mu = 1;
  (*bp).use_summ_mufine = 0;
  (*bp).use_summ_params = 1;
  (*bp).use_logspline = 1;
  strcpy((*bp).iter_knots_fname,"none");
  strcpy((*bp).iter_mu_fname,"samp_mu");
  strcpy((*bp).iter_mufine_fname,"none");
  strcpy((*bp).iter_params_fname,"samp_params");
  strcpy((*bp).summ_mu_fname,"summ_mu");
  strcpy((*bp).summ_mufine_fname,"none");
  strcpy((*bp).summ_params_fname,"summ_params");
/*
  dpar = tssdvec(1);
  dpar[0] = 6.0;
*/
  ipar = tssivec(2);
  ipar[0] = 1;
  ipar[1] = MAXKNOTS;
  imposeParamConstraints(bp,pp);
  pform = UNIFORM;
  Set_PriorParams(pp,bp,pform,dpar,ipar);
  free(ipar);
}












