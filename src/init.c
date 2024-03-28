#define USE_FC_LEN_T
#include "sdrt.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <stdlib.h>
#include "R_ext/BLAS.h"
#include <math.h>
#include <stddef.h> //for NULL

//void FMTS_interface(
//    const int* den_est,
//    const double* x,
//  const double* y,
//    const double* yful,
//    const int* n,
//    const int* p,
//    const double* ssww2,
//    double* dlogi,
//    double* den,
//    const int* out,
//    double* vecM,double* var
//){
//  int status1=dlogfmarg(x,n,p,var,out,den,dlogi);
//  int status2=vecMfmtscms(den_est,x,y,yful,n,
//                          p,ssww2,dlogi,vecM,var);
//  if(status1!=0){
//    error("Dlogfmarg has non zero exit");
//  }
//  else if(status2!=0)
//  {
//    error("vecMfmtscms has non zero exit");
//  }
//}
//static const R_CallMethodDef callMethods[]={
//  {"vecMfmtscms",(DL_FUNC) &vecMfmtscms,10},
//  {"vecMfmtscs",(DL_FUNC) &vecMfmtscs,11},
 // {"dlogfmarg",(DL_FUNC) &dlogfmarg,7},
 // {NULL,NULL,0}
//};

static R_NativePrimitiveArgType vecMfmtscms_t[] = {
  INTSXP,REALSXP,REALSXP,REALSXP,INTSXP,INTSXP,REALSXP,REALSXP,REALSXP,REALSXP
};
static R_NativePrimitiveArgType vecMfmtscs_t[] = {
  INTSXP,REALSXP,REALSXP,REALSXP,INTSXP,INTSXP,REALSXP,REALSXP,REALSXP,REALSXP,REALSXP
};
static R_NativePrimitiveArgType dlogfmarg_t[] = {
  REALSXP,INTSXP,INTSXP,REALSXP,INTSXP,REALSXP,REALSXP
};
static const R_CMethodDef cMethods[] = {
  {"vecMfmtscms", (DL_FUNC) &vecMfmtscms, 10,vecMfmtscms_t},
  {"vecMfmtscs", (DL_FUNC) &vecMfmtscs, 11,vecMfmtscs_t},
  {"dlogfmarg", (DL_FUNC) &dlogfmarg, 7,dlogfmarg_t},
  {NULL, NULL, 0, NULL}
};


//static R_NativePrimitiveArgType vecMfmn_t[] = {
//  REALSXP, REALSXP,INTSXP,INTSXP,REALSXP,REALSXP
//};

//static const R_CMethodDef cMethods[] = {
//  {"vecMfmn", (DL_FUNC) &dr2_interface, 6,vecMfmn_t},
//  {NULL, NULL, 0, NULL}
//};

//void R_init_dr2(DllInfo *info)
//{
//  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
//}


void R_init_sdrt(DllInfo *dll)
{
  R_registerRoutines(dll, cMethods, NULL
                       , NULL,NULL);
  R_useDynamicSymbols(dll, FALSE);
  //R_forceSymbols(dll, TRUE);
}
