static char mc_version[] = "MATLAB Compiler 1.0 infun";
/*
 *  MATLAB Compiler: 1.0
 *  Date: Oct 20, 1995
 *  Arguments: starts 
 */
#include <math.h>
#include "mex.h"
#include "mcc.h"


void
mexFunction(
    int nlhs_,
    Matrix *plhs_[],
    int nrhs_,
    Matrix *prhs_[]
)
{
   int ci_, i_, j_;
   unsigned flags_;
   Matrix *Mplhs_[32], *Mprhs_[32];
   for (ci_=i_=0; i_<nrhs_; ++i_)
   {
      if (prhs_[i_]->pi)
      {
         ci_ = 1;
         break;
      }
      if (prhs_[i_]->pr)
      {
         break;
      }
   }
   if (ci_)
   {
/***************** Compiler Assumptions ****************
 *
 *       C0_         	complex scalar temporary
 *       C1_         	complex scalar temporary
 *       I0_         	integer scalar temporary
 *       RM0_        	real vector/matrix temporary
 *       c           	integer scalar
 *       del         	complex vector/matrix
 *       i           	integer scalar
 *       len         	complex vector/matrix
 *       length      	<function>
 *       n           	integer scalar
 *       nargin      	<function>
 *       starts      	<function being defined>
 *       x           	complex vector/matrix
 *       y           	complex vector/matrix
 *       yi          	real vector/matrix
 *******************************************************/
      Matrix y;
      Matrix yi;
      Matrix x;
      Matrix len;
      Matrix del;
      int c;
      int n;
      int i;
      int I0_;
      double C0__r, C0__i;
      double C1__r, C1__i;
      Matrix RM0_;
      
      mccComplexInit(x);
      mccImport(&x, ((nrhs_>0) ? prhs_[0] : 0), 0, 0);
      mccComplexInit(len);
      mccImportCopy(&len, ((nrhs_>1) ? prhs_[1] : 0), 0, 0);
      mccComplexInit(del);
      mccImportCopy(&del, ((nrhs_>2) ? prhs_[2] : 0), 0, 0);
      mccComplexInit(y);
      mccRealInit(yi);
      mccRealInit(RM0_);
      
      /* % STARTS -- find sequences of entries in a vector */
      
      /* % [y, yi] = starts(x) returns a vector of entries and indices in x.   */
      /* % Each index marks the beginning of a sequence of numbers, each one 1  */
      /* % greater than the previous.  */
      
      /* % [y, yi] = starts(x, min, del) returns a vector of entries and */
      /* % indices of sequences in x, that have at least min adjacent elements */
      /* % separated by delta. (x, x+delta, x+2*delta ... x + (min-1)*delta ...) */
      
      /* % y = starts(find(x)) finds the beginnings of sequences of non zero  */
      /* % elements in x.  */
      
      /* if nargin < 2  len = 1; end */
      if ((mccNargin() < 2))
      {
         {
            double tr_ = 1;
            mccAllocateMatrix(&len, 1, 1);
            *len.pr = tr_;
         }
         *len.pi = 0.;
         len.dmode = mxNUMBER;
      }
      /* if nargin < 3  del = 1; end */
      if ((mccNargin() < 3))
      {
         {
            double tr_ = 1;
            mccAllocateMatrix(&del, 1, 1);
            *del.pr = tr_;
         }
         *del.pi = 0.;
         del.dmode = mxNUMBER;
      }
      
      /* y = []; */
      mccCreateEmpty(&y);
      /* yi = []; */
      mccCreateEmpty(&yi);
      
      /* c = 1;			% candidate index */
      c = 1;
      /* n = 1;			% number in sequence */
      n = 1;
      
      /* for i = 2:length(x) */
      if( x.flags & mccNOTSET )
      {
         mexErrMsgTxt( "variable x undefined, line 24" );
      }
      I0_ = mccGetLength(&x);
      for (i = 2; i <= I0_; i = i + 1)
      {
         /* if x(i) == (x(i-1) + del) */
         C0__r = (mccGetRealVectorElement(&x, (int)i));
         C0__i = mccGetImagVectorElement(&x, (int)i);
         C1__r = (mccGetRealVectorElement(&x, (int)(i - 1)));
         C1__i = mccGetImagVectorElement(&x, (int)(i - 1));
         if( del.flags & mccNOTSET )
         {
            mexErrMsgTxt( "variable del undefined, line 25" );
         }
         RM0_.dmode = mxNUMBER;
         {
            int m_=1, n_=1, cx_ = 0;
            double t_;
            double *p_RM0_;
            int I_RM0_=1;
            double *p_del;
            int I_del=1;
            double *q_del;
            m_ = mcmCalcResultSize(m_, &n_, del.m, del.n);
            mccAllocateMatrix(&RM0_, m_, n_);
            I_RM0_ = (RM0_.m != 1 || RM0_.n != 1);
            p_RM0_ = RM0_.pr;
            I_del = (del.m != 1 || del.n != 1);
            p_del = del.pr;
            q_del = del.pi;
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_del+=I_del, q_del+=I_del)
               {
                  *p_RM0_ = ((C0__r == (C1__r + *p_del)) && (C0__i == (C1__i + *q_del)));
                  ;
               }
            }
         }
         RM0_.dmode = mxNUMBER;
         if (mccIfCondition(&RM0_))
         {
            /* n = n+1; */
            n = (n + 1);
            /* else */
         }
         else
         {
            /* if n >= len */
            RM0_.dmode = mxNUMBER;
            {
               int m_=1, n_=1, cx_ = 0;
               double t_;
               double *p_RM0_;
               int I_RM0_=1;
               double *p_len;
               int I_len=1;
               double *q_len;
               m_ = mcmCalcResultSize(m_, &n_, len.m, len.n);
               mccAllocateMatrix(&RM0_, m_, n_);
               I_RM0_ = (RM0_.m != 1 || RM0_.n != 1);
               p_RM0_ = RM0_.pr;
               I_len = (len.m != 1 || len.n != 1);
               p_len = len.pr;
               q_len = len.pi;
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_len+=I_len, q_len+=I_len)
                  {
                     *p_RM0_ = (n >= *p_len);
                     ;
                  }
               }
            }
            RM0_.dmode = mxNUMBER;
            if (mccIfCondition(&RM0_))
            {
               /* y  = [y, x(c)]; */
               mccCatenateColumns(&y, &y, mccTempVectorElement(&x, c));
               /* yi = [yi, c]; */
               mccCatenateColumns(&yi, &yi, mccTempMatrix((double)(c), 0., mxNUMBER));
               /* end */
            }
            /* c = i; */
            c = i;
            /* n = 1; */
            n = 1;
            /* end */
         }
         /* end */
      }
      
      /* if n >= len */
      if( len.flags & mccNOTSET )
      {
         mexErrMsgTxt( "variable len undefined, line 37" );
      }
      RM0_.dmode = mxNUMBER;
      {
         int m_=1, n_=1, cx_ = 0;
         double t_;
         double *p_RM0_;
         int I_RM0_=1;
         double *p_len;
         int I_len=1;
         double *q_len;
         m_ = mcmCalcResultSize(m_, &n_, len.m, len.n);
         mccAllocateMatrix(&RM0_, m_, n_);
         I_RM0_ = (RM0_.m != 1 || RM0_.n != 1);
         p_RM0_ = RM0_.pr;
         I_len = (len.m != 1 || len.n != 1);
         p_len = len.pr;
         q_len = len.pi;
         for (j_=0; j_<n_; ++j_)
         {
            for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_len+=I_len, q_len+=I_len)
            {
               *p_RM0_ = (n >= *p_len);
               ;
            }
         }
      }
      RM0_.dmode = mxNUMBER;
      if (mccIfCondition(&RM0_))
      {
         /* y = [y, x(c)]; */
         mccCatenateColumns(&y, &y, mccTempVectorElement(&x, c));
         /* yi = [yi, c]; */
         mccCatenateColumns(&yi, &yi, mccTempMatrix((double)(c), 0., mxNUMBER));
         /* end */
      }
      
      
      mccReturnFirstValue(&plhs_[0], &y);
      mccReturnValue(&plhs_[1], &yi);
   }
   else
   {
/***************** Compiler Assumptions ****************
 *
 *       I0_         	integer scalar temporary
 *       R0_         	real scalar temporary
 *       R1_         	real scalar temporary
 *       RM0_        	real vector/matrix temporary
 *       c           	integer scalar
 *       del         	real vector/matrix
 *       i           	integer scalar
 *       len         	real vector/matrix
 *       length      	<function>
 *       n           	integer scalar
 *       nargin      	<function>
 *       starts      	<function being defined>
 *       x           	real vector/matrix
 *       y           	real vector/matrix
 *       yi          	real vector/matrix
 *******************************************************/
      Matrix y;
      Matrix yi;
      Matrix x;
      Matrix len;
      Matrix del;
      int c;
      int n;
      int i;
      int I0_;
      double R0_;
      double R1_;
      Matrix RM0_;
      
      mccRealInit(x);
      mccImport(&x, ((nrhs_>0) ? prhs_[0] : 0), 0, 0);
      mccRealInit(len);
      mccImportCopy(&len, ((nrhs_>1) ? prhs_[1] : 0), 0, 0);
      mccRealInit(del);
      mccImportCopy(&del, ((nrhs_>2) ? prhs_[2] : 0), 0, 0);
      mccRealInit(y);
      mccRealInit(yi);
      mccRealInit(RM0_);
      
      /* % STARTS -- find sequences of entries in a vector */
      
      /* % [y, yi] = starts(x) returns a vector of entries and indices in x.   */
      /* % Each index marks the beginning of a sequence of numbers, each one 1  */
      /* % greater than the previous.  */
      
      /* % [y, yi] = starts(x, min, del) returns a vector of entries and */
      /* % indices of sequences in x, that have at least min adjacent elements */
      /* % separated by delta. (x, x+delta, x+2*delta ... x + (min-1)*delta ...) */
      
      /* % y = starts(find(x)) finds the beginnings of sequences of non zero  */
      /* % elements in x.  */
      
      /* if nargin < 2  len = 1; end */
      if ((mccNargin() < 2))
      {
         {
            double tr_ = 1;
            mccAllocateMatrix(&len, 1, 1);
            *len.pr = tr_;
         }
         len.dmode = mxNUMBER;
      }
      /* if nargin < 3  del = 1; end */
      if ((mccNargin() < 3))
      {
         {
            double tr_ = 1;
            mccAllocateMatrix(&del, 1, 1);
            *del.pr = tr_;
         }
         del.dmode = mxNUMBER;
      }
      
      /* y = []; */
      mccCreateEmpty(&y);
      /* yi = []; */
      mccCreateEmpty(&yi);
      
      /* c = 1;			% candidate index */
      c = 1;
      /* n = 1;			% number in sequence */
      n = 1;
      
      /* for i = 2:length(x) */
      if( x.flags & mccNOTSET )
      {
         mexErrMsgTxt( "variable x undefined, line 24" );
      }
      I0_ = mccGetLength(&x);
      for (i = 2; i <= I0_; i = i + 1)
      {
         /* if x(i) == (x(i-1) + del) */
         R0_ = (mccGetRealVectorElement(&x, i));
         R1_ = (mccGetRealVectorElement(&x, (i-1)));
         if( del.flags & mccNOTSET )
         {
            mexErrMsgTxt( "variable del undefined, line 25" );
         }
         RM0_.dmode = mxNUMBER;
         {
            int m_=1, n_=1, cx_ = 0;
            double t_;
            double *p_RM0_;
            int I_RM0_=1;
            double *p_del;
            int I_del=1;
            m_ = mcmCalcResultSize(m_, &n_, del.m, del.n);
            mccAllocateMatrix(&RM0_, m_, n_);
            I_RM0_ = (RM0_.m != 1 || RM0_.n != 1);
            p_RM0_ = RM0_.pr;
            I_del = (del.m != 1 || del.n != 1);
            p_del = del.pr;
            for (j_=0; j_<n_; ++j_)
            {
               for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_del+=I_del)
               {
                  *p_RM0_ = (R0_ == (R1_ + *p_del));
                  ;
               }
            }
         }
         RM0_.dmode = mxNUMBER;
         if (mccIfCondition(&RM0_))
         {
            /* n = n+1; */
            n = (n + 1);
            /* else */
         }
         else
         {
            /* if n >= len */
            RM0_.dmode = mxNUMBER;
            {
               int m_=1, n_=1, cx_ = 0;
               double t_;
               double *p_RM0_;
               int I_RM0_=1;
               double *p_len;
               int I_len=1;
               m_ = mcmCalcResultSize(m_, &n_, len.m, len.n);
               mccAllocateMatrix(&RM0_, m_, n_);
               I_RM0_ = (RM0_.m != 1 || RM0_.n != 1);
               p_RM0_ = RM0_.pr;
               I_len = (len.m != 1 || len.n != 1);
               p_len = len.pr;
               for (j_=0; j_<n_; ++j_)
               {
                  for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_len+=I_len)
                  {
                     *p_RM0_ = (n >= *p_len);
                     ;
                  }
               }
            }
            RM0_.dmode = mxNUMBER;
            if (mccIfCondition(&RM0_))
            {
               /* y  = [y, x(c)]; */
               mccCatenateColumns(&y, &y, mccTempVectorElement(&x, c));
               /* yi = [yi, c]; */
               mccCatenateColumns(&yi, &yi, mccTempMatrix((double)(c), 0., mxNUMBER));
               /* end */
            }
            /* c = i; */
            c = i;
            /* n = 1; */
            n = 1;
            /* end */
         }
         /* end */
      }
      
      /* if n >= len */
      if( len.flags & mccNOTSET )
      {
         mexErrMsgTxt( "variable len undefined, line 37" );
      }
      RM0_.dmode = mxNUMBER;
      {
         int m_=1, n_=1, cx_ = 0;
         double t_;
         double *p_RM0_;
         int I_RM0_=1;
         double *p_len;
         int I_len=1;
         m_ = mcmCalcResultSize(m_, &n_, len.m, len.n);
         mccAllocateMatrix(&RM0_, m_, n_);
         I_RM0_ = (RM0_.m != 1 || RM0_.n != 1);
         p_RM0_ = RM0_.pr;
         I_len = (len.m != 1 || len.n != 1);
         p_len = len.pr;
         for (j_=0; j_<n_; ++j_)
         {
            for (i_=0; i_<m_; ++i_, p_RM0_+=I_RM0_, p_len+=I_len)
            {
               *p_RM0_ = (n >= *p_len);
               ;
            }
         }
      }
      RM0_.dmode = mxNUMBER;
      if (mccIfCondition(&RM0_))
      {
         /* y = [y, x(c)]; */
         mccCatenateColumns(&y, &y, mccTempVectorElement(&x, c));
         /* yi = [yi, c]; */
         mccCatenateColumns(&yi, &yi, mccTempMatrix((double)(c), 0., mxNUMBER));
         /* end */
      }
      
      
      mccReturnFirstValue(&plhs_[0], &y);
      mccReturnValue(&plhs_[1], &yi);
   }
   return;
}
