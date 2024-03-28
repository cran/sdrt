//----------------------------------------------------------------//
// file name: FMTS-clibs.c                                        //
//                                                                //
// C functions needed by FMTS.R                                   //
// To complie this file,                                          //
//       use "Rcmd SHLIB FMTS-clibs.c -lRblas" (windows)          //
//      or "R CMD SHLIB FMTS-clibs.c" (unix/linux)                //
//                                                                //
// Reference:                                                     //
//   Samadi and Priyan (2020).                                    //
//    Fourier Methods for Estimating the Central Subspace        //
//     in Time series                                             //
//                                                                //
// 09/21/2019                                                     //
//----------------------------------------------------------------//
#define USE_FC_LEN_T
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef FCONE
# define FCONE
#endif
/*=========================== Interface ===========================*/

void vecMfmtscms(const int* den_est,const double* x, const double* y, const double* yful,const int* nn,
              const int* pp, const double* ssww2,const double* dlogi,
              double* vecM,double* var);
void vecMfmtscs(const int* den_est,const double* x, const double* y, const double* yful, const int* nn,
                 const int* pp, const double* ssww2,const double* sstt2,const double* dlogi,
                 double* vecM,double* var);
/*====================== supported function =======================*/

void Afm(const double* x, const double* y, const int* nn,
         const int* pp, const double* ssww2, double* A,double* var);

void Afc(const double* x, const double* y, const int* nn,
         const int* pp, const double* ssww2, const double* sstt2,
         double* A);

void Mnorm(const double* x, const double* A, const int* nn,
           const int* pp, const double* ssww2,
           const double* scale,const double* dlogi,
           double* vecM,double* var);

void Mnormal(const double* x, const double* A, const int* nn,
             const int* pp, const double* ssww2, const double* yful, const double* scale,
             const double* dlogi, double* vecM,double* var);

void IplusMnorm(const int* den_est,const double* x, const double* A, const int* nn,
                const int* pp, const double* ssww2, const double* yful,
                const double* scale,const double* dlogi, double* vecM,double* var);

/*======================== CODE BEGIN HERE ========================*/
/*=================================================================*/
/*      calculate matrix M_{FMTSCS},  output a p-by-p matrix       */
/*=================================================================*/

void vecMfmtscs(const int* den_est,const double* x, const double* y,const double* yful, const int* nn,
               const int* pp, const double* ssww2,const double* sstt2,const double* dlogi,
               double* vecM,double* var)
{
  const double one = 1.0;
  double* A;
  A = (double*)calloc((*nn)*(*nn), sizeof(double));

  Afc(x, y, nn, pp, ssww2, sstt2, A);

  IplusMnorm(den_est,x, A, nn, pp, ssww2,yful ,&one,dlogi,vecM,var);
  free(A);
}

/*=================================================================*/
/*      calculate matrix M_{FMTSCMS},  output a p-by-p matrix      */
/*=================================================================*/


void vecMfmtscms(const int* den_est,const double* x, const double* y, const double* yful, const int* nn,
              const int* pp, const double* ssww2,const double* dlogi,
              double* vecM,double* var)
{
  const double one = 1.0;
  double* A;
  A = (double*)calloc((*nn)*(*nn), sizeof(double));
  Afm(x, y, nn, pp, ssww2, A,var);
  IplusMnorm(den_est,x, A, nn, pp, ssww2, yful,&one,dlogi,vecM,var);
  free(A);
}

/********************************************************************/
/********************************************************************/
/**********************************************************************
 * Estimation derivative of log joint density
 * using Gaussian kernel function.
 *
 * inputs:
 *      x    ---- n-by-p matrix(data)
 *      out  ---- 1 or 0 (1: exclude i==j from summation, 0: otherwise)
 *
 * Output:
 *    dlogfjoint  --  1-by-p matrix  (estimated joint d-log-den at y_{s-1}\y_{t-1})
 **********************************************************************/
void dlogfjoint(const double* x, const int* n, const int* p, const int* t,
            const int* s, const int* out,
            double den, double* dlogfjoint,double* var)
{

  const int nn=(*n);
  const int k=abs((*s)-(*t));
  const int ss=(*s);
  const int tt=(*t);
  const int pp=(*p);
  const int d=pp+k;
  double a, b, sum=0,sumtt=0,sum2=0.000001,sum3=0.000001,sumt=0,sums=0;
  double xnew[d];
  int i,j;
  double an1[*p], sum1[d],sum22[pp],dlogk[pp],den1=0;
  double d4=d+4;
  double d2=d+2;
  for(j=0;j<pp;j++){
    an1[j] =var[j]*pow(nn,(-1/(d4)))*pow((4/d2),(1/(d2)));
  }
  //printf(" %f ",an1[1]);
  for(i=0;i<d;i++){
    sum1[i]=0;
    if(i<pp){
      sum22[i]=0;
      xnew[i]=x[ss+i*nn];
    }else{
      xnew[i]=x[tt+(i-k)*nn];
    }
  }

  for(j=0;j<nn;j++){
    sumtt=0;
    sums=0;
    sumt=0;
    sum=0;
    for(i=0;i<d;i++){
      if(i<pp){
        a=-0.5/an1[i]/an1[i];
        sums=sums+(a*(xnew[i]-x[j+i*nn])*(xnew[i]-x[j+i*nn]));
        sumtt+=a*(x[tt+i*nn]-x[j+i*nn])*(x[tt+i*nn]-x[j+i*nn]);
      }else{
        a=-0.5/an1[i-k]/an1[i-k];
        sumt=sumt+(a*(xnew[i]-x[j+(i-k)*nn])*(xnew[i]-x[j+(i-k)*nn]));

      }
    }

    sum3+=exp(sumtt);
    sum=sums+sumt;
    sum2+=exp(sum);
      for(i=0;i<d;i++){
        if(i<pp){
           sum1[i]+=exp(sum)*((x[j+i*nn]));
          sum22[i]+=exp(sumtt)*((x[j+i*nn]));
        }else{
           sum1[i]+=exp(sum)*((x[j+(i-k)*nn]));
        }
      }
  }// End of sample loop

  den=(sum2);

  den1=(sum3);
  if((*out) == 0)
  {
    for(i = 0; i < pp; i++){
      dlogfjoint[i] = xnew[i];
      dlogk[i] = x[tt+i*nn];
      //b = -1.0 / an1[i] /an1[i];
      b = -1.0 / an1[i];
      a = b / (den);
      dlogfjoint[i] =dlogfjoint[i]*b-a*sum1[i];
      a = b / (den1);
      dlogk[i] =dlogk[i]*b-a*sum22[i];
    }
  }
  else
  {
    for(i = 0; i < pp; i++){
      dlogfjoint[i] = 0.0;
    //b = -1.0 / an1[i] /an1[i];
      b = -1.0 / an1[i];
      a = b / (den);
      dlogfjoint[i] =dlogfjoint[i]*b-a*sum1[i];
      a = b / (den1);
      dlogk[i] =dlogk[i]*b-a*sum22[i];
    }
  }
 for(i=0;i<pp;i++){
    if(i<=k){
      dlogfjoint[i]=dlogfjoint[i];
    }else{
      dlogfjoint[i]=dlogfjoint[i]-dlogk[i];
    }
  }
}

/**********************************************************************
 * Estimation derivative of log density
 * using Gaussian kernel function.
 *
 * inputs:
 *      x    ---- n-by-p matrix(data)
 *      out  ---- 1 or 0 (1: exclude i==j from summation, 0: otherwise)
 *
 * Output:
 *    dlogi  --  n-by-p matrix  (estimated d-log-den at y_{t-1})
 **********************************************************************/

void dlogfmarg(const double* x, const int* n, const int* p,double* var,
               const int* out,double* den, double* dlogi)
{
  const int nn=(*n);
  const int pp=(*p);
  const int np=nn*pp;
  //const int tt=(*ni);
  double sc;
  double a, b,sum=0,sum2=0.000001;
  double prod=1;
  //double wio;
  double sum1[np];
  int i, j,k;
  double an1[*p];//,xixi[nn];
  double p4=(*p)+4;//,p2=pp+2;

  for(j=0;j<(*p);j++){
    an1[j] =var[j]*pow((*n),(-1/(p4)));//*pow((4/p2),(1/(p2)));
    //sum1[j]=0;
    prod=prod*an1[j];
  }
  sc= R_pow_di(M_1_SQRT_2PI, *p) / (double)((*n) - (*out))/prod;
 //*************************New codes**************************
 for(i = 0; i < np; i++){
      sum1[i]=0.0;
 }

  for(i = 0; i < (*n); i++){
    sum2=0;
    for(j=0;j<(*n);j++){
      sum=0;
      for(k=0;k<(*p);k++){
        a=-0.5/an1[k]/an1[k];
        sum=sum+a*(x[i+k*nn]-x[j+k*nn])*(x[i+k*nn]-x[j+k*nn]);
      }
      sum2=sum2+exp(sum);

      for(k=0;k<pp;k++){
        sum1[i+k*nn]+=exp(sum)*((x[j+k*nn]));
      }
    }//end of sample

    den[i]=sum2;


  for(k=0;k<pp;k++)dlogi[i+k*nn] =x[i+k*nn];
  for(k = 0; k < (*p); k++)
  {
    //b = -1.0 / an1[k] /an1[k];
    b = -1.0 / an1[k];
    a = b / (den[i]);
    dlogi[i+k*nn] =dlogi[i+k*nn]*b -a*sum1[i+k*nn];
  }
  den[i]=den[i]*sc;
  }//end of loop for i
}
/*========================== SUBROUTINES ==========================*/

/*=================================================================*/
/* Calculate matrix of type                                        */
/* $$ n^{-2} \sum_{i, j} A_{i, j}                                  */
/*   ((sw2 + 1) x_i - sw2 x_j) ((sw2 + 1) x_j - sw2 x_i)^{\tau} $$ */
/*                                                                 */
/* x is a n-by-p matrix (or n*p vector)                            */
/* A is a n-by-n matrix (or n*n vector)                            */
/* n is the sample size and p is the dimension                     */
/* sw2 is the tuning parameter                                     */
/* scale, times the resulted matrix by scale                       */
/* result are stored in vecM, which is a p-by-p matrix             */
/*                                                                 */
/* for (i, j)th element of vecM                                    */
/* first calculate                                                 */
/*  sum_{ni < nj}                                                  */
/*  A[ni, nj] * (sc1 * (x[ni,i] * x[nj,j] + x[nj,i] * x[ni,j]))    */
/*                                                                 */
/* then plus                                                       */
/*  (sc1 * A[ni, ni] - sc2 * (ni rowsum)) * x[ni,i] * x[ni,j]      */
/*=================================================================*/

void Mnorm(const double* x, const double* A, const int* nn,
           const int* pp, const double* ssww2, const double* scale,
           const double* dlogi,double* vecM,double* var)
{
  const int n = (*nn);
  const int p = (*pp);
  const int out=0;
  const double sw2 = (*ssww2);
  const double ndiv = (n * n) / (*scale);
  int i, j, ni, nj;

  double temp;
  double* dlogfij;
  double den=0.0;

  dlogfij= (double*)calloc((*pp), sizeof(double));
  /* initialize vecM */
  for(i = 0; i < (p*p); i++) vecM[i] = 0.0;

  /* begin to calculate vecM */

  for(i=0;i<p;i++){
    for(j=0;j<p;j++){
      temp=0.0;
      for(ni=0;ni<n;ni++){ // ni do the roll of t
        for(nj=0;nj<n;nj++){//nj do the roll of s
          if((abs(ni-nj)<p) & (abs(ni-nj)!=0)){
             dlogfjoint(x,&n, &p,&ni,&nj,&out,den,dlogfij,var);
             //temp=temp+A[ni+nj*n]*(dlogi[i]*dlogfij[j]+sw2*dlogi[i]*(x[ni+j*n]-x[nj+j*n])-sw2*dlogfij[j]*(x[ni+i*n]-x[nj+i*n])-sw2*sw2*(x[ni+i*n]*x[ni+j*n]-x[ni+i*n]*x[nj+j*n]-x[nj+i*n]*x[ni+j*n]+x[nj+i*n]*x[nj+j*n]));
             temp=temp+A[ni+nj*n]*(dlogi[ni+i*n]*dlogfij[j]+sw2*dlogi[ni+i*n]*(x[ni+j*n]-x[nj+j*n])-sw2*dlogfij[j]*(x[ni+i*n]-x[nj+i*n])-sw2*sw2*(x[ni+i*n]*x[ni+j*n]-x[ni+i*n]*x[nj+j*n]-x[nj+i*n]*x[ni+j*n]+x[nj+i*n]*x[nj+j*n]));
          }else{
              //temp=temp+A[ni+nj*n]*(dlogi[i]*dlogi[j]+sw2*dlogi[i]*(x[ni+j*n]-x[nj+j*n])-sw2*dlogi[j]*(x[ni+i*n]-x[nj+i*n])-sw2*sw2*(x[ni+i*n]*x[ni+j*n]-x[ni+i*n]*x[nj+j*n]-x[nj+i*n]*x[ni+j*n]+x[nj+i*n]*x[nj+j*n]));
                temp=temp+A[ni+nj*n]*(dlogi[ni+i*n]*dlogi[nj+j*n]+sw2*dlogi[ni+i*n]*(x[ni+j*n]-x[nj+j*n])-sw2*dlogi[nj+j*n]*(x[ni+i*n]-x[nj+i*n])-sw2*sw2*(x[ni+i*n]*x[ni+j*n]-x[ni+i*n]*x[nj+j*n]-x[nj+i*n]*x[ni+j*n]+x[nj+i*n]*x[nj+j*n]));

          }

        }//end of nj loop
      }//end of ni loop
      vecM[j+i*p]=temp/ndiv;

    }//end of j loop
  }//end i loop
}
/*=================================================================*/
/*For calculating Determinant of the Matrix */
int LUPDecompose(double *A, int N, double Tol, int *P) {

  int i, j, k, imax;
  double maxA, ptr, absA;

  for (i = 0; i <= N; i++)
    P[i] = i; //Unit permutation matrix, P[N] initialized with N

  for (i = 0; i < N; i++) {
    maxA = 0.0;
    imax = i;

    for (k = i; k < N; k++){
      if ((absA = fabs(A[k+i*N])) > maxA) {
        maxA = absA;
        imax = k;
      }

      if (maxA < Tol){ return 0;
      }
      //failure, matrix is degenerate

      if (imax != i) {
        //pivoting P
        j = P[i];
        P[i] = P[imax];
        P[imax] = j;

        //pivoting rows of A
        ptr = A[i];
        A[i] = A[imax];
        A[imax] = ptr;

        //counting pivots starting from N (for determinant)
        P[N]++;
      }
    }
      for (j = i + 1; j < N; j++) {
        A[j+i*N] /= A[i+i*N];

        for (k = i + 1; k < N; k++)
          A[j+k*N] -= A[j+i*N] * A[i+k*N];
      }
  }

  return 1;  //decomposition done
}

/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
 * OUTPUT: x - solution vector of A*x=b
 */
void LUPSolve(double **A, int *P, double *b, int N, double *x) {

  for (int i = 0; i < N; i++) {
    x[i] = b[P[i]];

    for (int k = 0; k < i; k++)
      x[i] -= A[i][k] * x[k];
  }

  for (int i = N - 1; i >= 0; i--) {
    for (int k = i + 1; k < N; k++)
      x[i] -= A[i][k] * x[k];

    x[i] /= A[i][i];
  }
}

/* INPUT: A,P filled in LUPDecompose; N - dimension
 * OUTPUT: IA is the inverse of the initial matrix
 */
void LUPInvert(double *A, int *P, int N, double *IA) {

  for (int j = 0; j < N; j++) {
    for (int i = 0; i < N; i++) {
      IA[i+j*N] = P[i] == j? 1.0 : 0.0;

      for (int k = 0; k < i; k++)
        IA[i+j*N] -= A[i+k*N] * IA[k+j*N];
    }

    for (int i = N - 1; i >= 0; i--) {
      for (int k = i + 1; k < N; k++)
        IA[i+j*N] -= A[i+k*N] * IA[k+j*N];

      IA[i+j*N] /= A[i+i*N];
    }
  }
}


/* INPUT: A,P filled in LUPDecompose; N - dimension.
 * OUTPUT: Function returns the determinant of the initial matrix
 */
double LUPDeterminant(double *A, int *P, int N) {

  double det = A[0];

  for (int i = 1; i < N; i++)
    det *= A[i+i*N];
   //if((P[N] - N)%2==0){
     //return(det);
   //}else
  return (P[N] - N)%2 == 0?det:-det;
}
void Gk(const double* x,const int* nn,const int* pp,const double* sigma,double* gki)
{
  int i,j,ni;
  int p=(*pp);
  int n=(*nn);

  double temp;
  for(ni=0;ni<n;ni++){
    for(i=0;i<(*pp);i++){
      gki[ni+i*n]=0;
    }
  }
  for(ni=0;ni<n;ni++){
    for(i=0;i<(*pp);i++){
    temp=0;
        for(j=0;j<(*pp);j++){
        temp+=x[ni+j*n]*sigma[j+i*p];
      }
    gki[ni+i*n]=temp;
    }
  }
}
void dspyy(const double* x, const int* nn,const int* pp,
           double* inv_sigmayy)
{
  int i,j,t;
  double temp1;
  double* sigmayy;
  //double d=0;
  double Tol=0.000001;
  int n=(*nn);
  int p=(*pp);
  int P[n];
  sigmayy= (double*)calloc(p*p, sizeof(double));

  for(i=0;i<p;i++){
    for(j=0;j<p;j++){
      sigmayy[j+i*p]=0;
      inv_sigmayy[j+i*p]=0;

    }
  }
  for(i=0;i<p;i++){
    for(j=0;j<p;j++){
      temp1=0.0;
      for(t=0;t<n;t++){
        temp1+=x[t+i*n]*x[t+j*n];
      }
      sigmayy[j+i*p]=temp1/n;
      //sigmayy[j+i*p]=temp1/n;
      //printf("%f",sigmayy[i+j*p]);
      //inv_sigmayy[i+j*p]=sigmayy[i+j*p];
      //inv_sigmayy[j+i*p]=sigmayy[j+i*p];
    }
  }
  //cofactor(sigmayy, *pp,inv_sigmayy);
  LUPDecompose(sigmayy,p,Tol,P);
  //d=LUPDeterminant(sigmayy,P,p);
  //printf("%f\t",d);
  LUPInvert(sigmayy, P, p, inv_sigmayy);
  free(sigmayy);
}
void dspxy(const double* x,const double* y, const int* nn,const int* pp,
           const int* ni, const int* nj,double* inv_sigmaxy)
{
  int t=(*ni);
  int s=(*nj);
  int i,j,k,r,h=0;
  int n=(*nn);
  int p=(*pp);
  //double d=0;
  double Tol=0;
  int P[n];
  double temp1;
  double* sigmaxx;
  double* inv_sigmaxx;
  double* A;
  sigmaxx= (double*)calloc(p*p, sizeof(double));
  inv_sigmaxx= (double*)calloc(p*p, sizeof(double));
  A= (double*)calloc(p*p, sizeof(double));
  k=abs(t-s);
  for(i=0;i<p;i++){
    for(j=0;j<p;j++){
      sigmaxx[j+i*p]=0;
      inv_sigmaxx[j+i*p]=0;
      inv_sigmaxy[j+i*p]=0;
    }
  }

  for(i=0;i<(*pp);i++){
    for(j=0;j<(*pp);j++){

      if(i>=j){

        //for(l=k;l<(k+p);l++){
           h=k+(i-j);
          temp1=0;
           for(r=0;r<(n-h);r++){
              temp1+=y[r]*y[r+h];
            }
           //printf("%f\t",temp1);
            inv_sigmaxy[j+i*p]=temp1/(n-h);
        //}
      }else{

        //for(l=j;l<p;l++){
          h=k-j;
          temp1=0;
          for(r=0;r<(n-h);r++){
            temp1+=y[r]*y[r+h];
          }
          //printf("%f\t",temp1);
          inv_sigmaxy[j+i*p]=temp1/(n-h);
        //}
      }
      //printf("%f\t",temp1);

      //sigmaxy[j+i*p]=temp1/(n-k);

    }
  }
  dspyy(x,&n,&p,sigmaxx);
  //cofactor(sigmaxx, *pp,inv_sigmaxx);
  LUPDecompose(sigmaxx,p,Tol,P);

  LUPInvert(sigmaxx, P, p, inv_sigmaxx);

  for(i=0;i<p;i++){
    temp1=0;
    for(j=0;j<p;j++){
      temp1+=inv_sigmaxy[i+j*p]*inv_sigmaxx[j+i*p];
    }
    //printf("%f\t",temp1);
    A[i]=temp1;
  }

  for(i=0;i<p;i++){
    temp1=0;
    for(j=0;j<p;j++){
      temp1+=A[i+j*p]*inv_sigmaxy[j+i*p];
    }
    inv_sigmaxy[i]=sigmaxx[i]-temp1;
  }
  //cofactor(inv_sigmaxy, *pp,inv_sigmaxy);
  LUPDecompose(inv_sigmaxy,p,Tol,P);
  //d=LUPDeterminant(inv_sigmaxy,P,p);

  //if(d>0 || d<0){
    //printf("%f\t",d);
    LUPInvert(inv_sigmaxy, P, p, inv_sigmaxy);
  //}else{
    //inv_sigmaxy[0]=inv_sigmaxy[0]+0.1;
    //LUPInvert(sigmaxx, P, p, inv_sigmaxy);
  //}
  free(A);
  free(sigmaxx);

}


/*=================================================================*/
/* Calculate matrix of type                                        */
/* $$ n^{-2} \sum_{i, j} A_{i, j}                                  */
/*   ((sw2 + 1) x_i - sw2 x_j) ((sw2 + 1) x_j - sw2 x_i)^{\tau} $$ */
/*                                                                 */
/* x is a n-by-p matrix (or n*p vector)                            */
/* A is a n-by-n matrix (or n*n vector)                            */
/* n is the sample size and p is the dimension                     */
/* sw2 is the tuning parameter                                     */
/* scale, times the resulted matrix by scale                       */
/* result are stored in vecM, which is a p-by-p matrix             */
/*                                                                 */
/* for (i, j)th element of vecM                                    */
/* first calculate                                                 */
/*  sum_{ni < nj}                                                  */
/*  A[ni, nj] * (sc1 * (x[ni,i] * x[nj,j] + x[nj,i] * x[ni,j]))    */
/*                                                                 */
/* then plus                                                       */
/*  (sc1 * A[ni, ni] - sc2 * (ni rowsum)) * x[ni,i] * x[ni,j]      */
/*=================================================================*/

void Mnormal(const double* x, const double* A, const int* nn,
           const int* pp, const double* ssww2,const double* yful,
           const double* scale,const double* dlogi,
           double* vecM,double* var)
{
  const int n = (*nn);
  const int p = (*pp);
  //const int out=0;
  const double sw2 = (*ssww2);
  const double ndiv = (n * n) / (*scale);
  int i, j, ni, nj;
  //double dlogi[*pp];
  double temp,temp1,temp2;
  double* dlogfij;
  double* dlogfii;
  double* gki;
  double* gkij;
  //double den=0.0;

  dlogfij= (double*)calloc(p*p, sizeof(double));
  dlogfii= (double*)calloc(p*p, sizeof(double));
  gki= (double*)calloc(n*p, sizeof(double));
  gkij= (double*)calloc(n*p, sizeof(double));
  /* initialize vecM */
  for(i = 0; i < (p*p); i++) vecM[i] = 0.0;

  /* begin to calculate vecM */

  for(i=0;i<p;i++){
    for(j=0;j<p;j++){
      temp=0.0;
      dspyy(x,&n,&p,dlogfii);


      for(ni=0;ni<n;ni++){ // ni do the roll of t
        Gk(x,&n,&p,dlogfii,gki);
        for(nj=0;nj<n;nj++){//nj do the roll of s

          //printf("%f",dlogfij[ni+i*n]);
          temp1=-gki[ni+i*n];//-dlogi[ni+i*n];//x[ni+i*p];//
          if((abs(ni-nj)<p) & (abs(ni-nj)!=0)){
            //dlogfjoint(x,&n, &p,&ni,&nj,&out,den,dlogfij,var);
            //if(j<=abs(ni-nj)){
            dspxy(x,yful,&n,&p,&ni,&nj,dlogfij);
            Gk(x,&n,&p,dlogfij,gkij);
               temp2=-gkij[nj+j*n];//-gkij[ni+i*n];
            //printf("%f\t",dlogfij[j]);
            //}else{
              // temp2=-x[nj+j*p];
          //  //}
               temp=temp+A[ni+nj*n]*(temp1*temp2+sw2*temp1*(x[ni+j*n]-x[nj+j*n])-sw2*temp2*(x[ni+i*n]-x[nj+i*n])-sw2*sw2*(x[ni+i*n]*x[ni+j*n]-x[ni+i*n]*x[nj+j*n]-x[nj+i*n]*x[ni+j*n]+x[nj+i*n]*x[nj+j*n]));
          }else{
            //dlogi[i]=-x[ni+j*n];
            temp2=-gki[nj+j*n];//-gkij[j];
            //temp2=-x[nj+j*p];
            temp=temp+A[ni+nj*n]*(temp1*temp2+sw2*temp1*(x[ni+j*n]-x[nj+j*n])-sw2*temp2*(x[ni+i*n]-x[nj+i*n])-sw2*sw2*(x[ni+i*n]*x[ni+j*n]-x[ni+i*n]*x[nj+j*n]-x[nj+i*n]*x[ni+j*n]+x[nj+i*n]*x[nj+j*n]));
          }

        }//end of nj loop
      }//end of ni loop

      vecM[j+i*p]=temp/ndiv;
      //printf("%f\t",temp);
    }//end of j loop
  }//end i loop
}

/*=================================================================*/
void Afc(const double* x, const double* y, const int* nn,
         const int* pp, const double* ssww2, const double* sstt2,
         double* A)
{
  const int n = (*nn);
  const int p = (*pp);
  const double sw2 = (*ssww2);
  const double st2 = (*sstt2);
  int i, j, k;
  double temp1, temp2, temp3;

  for(i = 0; i < (n*n); i+=(n+1)) A[i] = 1.0;

  /* calculate the distance matrix A */
  for(i = 0; i < n; i++)
    for(j = (i+1); j < n; j++)
    {
      temp1 = 0.0;
      for(k = 0; k < p; k++)
      {
        temp2 = (x[i + k*n] - x[j + k*n]);
        temp1 += (temp2 * temp2);
      }

      temp2 = (y[i] - y[j]);
      temp3 = exp(-0.5 * (st2 * temp2 * temp2 + sw2 * temp1));
      A[i + j*n] = temp3;
      A[j + i*n] = temp3;

    }
}
/*=================================================================*/
/* calculate the weight matrix                                     */
/* the (i, j)th element is                                         */
/*    $ y[i] * y[j] * exp(- 0.5 * sw2 * |x[i] - x[j]|^2) $         */
/*=================================================================*/

void Afm(const double* x, const double* y, const int* nn,
         const int* pp, const double* ssww2, double* A,double* var)
{
  const int n = (*nn);
  const int p = (*pp);
  //const int one=1;
  const double sw2 = (*ssww2);
  //const double b=-1;
  int i, j;
  int k;
  double temp1, temp2;
  double temp3;
  //double xij[p];
  //double wij;

  for(i = 0; i < n; i++) A[i + i*n] = y[i] * y[i];

  /* calculate the distance matrix A */
  for(i = 0; i < n; i++)
    for(j = (i+1); j < n; j++)
    {
      temp1 = 0.0;
      for(k = 0; k < p; k++)
      {
        temp2 = (x[i + k*n] - x[j + k*n]);
        temp1 += (temp2 * temp2);
      }
      temp3 = y[i] * y[j] * exp(-0.5 * sw2 * temp1);

      A[i + j*n] = temp3;
      A[j + i*n] = temp3;
    }
}
/*=================================================================*/
/* Calculate matrix of type                                        */
/* $$ n^{-2} \sum_{i, j} A_{i, j} [ sw2 I_p                        */
/*  ((sw2 + 1) x_i - sw2 x_j) ((sw2 + 1) x_j - sw2 x_i)^{\tau} ]$$ */
/*                                                                 */
/* x is a n-by-p matrix (or n*p vector)                            */
/* A is a n-by-n matrix (or n*n vector)                            */
/* n is the sample size and p is the dimension                     */
/* sw2 is the tuning parameter                                     */
/* scale, times the resulted matrix by scale                       */
/* result are stored in vecM, which is a p-by-p matrix             */
/*=================================================================*/

void IplusMnorm(const int* den_est,const double* x, const double* A, const int* nn,
                const int* pp, const double* ssww2, const double* yful,
                const double* scale,const double* dlogi,double* vecM,double* var)
{
  const int n = (*nn);
  const int p = (*pp);
  const int denest=(*den_est);
  int i;
  double temp1, Asum;
  if(denest==1){
    Mnorm(x, A, nn, pp, ssww2, scale,dlogi, vecM,var);
  }else{
    Mnormal(x, A, nn, pp, ssww2,yful, scale,dlogi,vecM,var);
  }
  Asum = 0.0;
  for(i = 0; i < (n*n); i++) Asum += A[i];
  temp1 = (*scale) * Asum * (*ssww2) / (n * n);
  for(i = 0; i < p; i++) vecM[i + i*p] += temp1;
}

//*********************************** END OF FILE **************************

