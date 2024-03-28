/*=========== Interface =============*/

void vecMfmtscms(
    const int* den_est,
    const double* x,
    const double* y,
    const double* yful,
    const int* nn,
    const int* pp,
    const double* ssww2,
    const double* dlogi,
    double* vecM,double* var);

void vecMfmtscs(
    const int* den_est,
    const double* x,
    const double* y,
    const double* yful,
    const int* nn,
    const int* pp,
    const double* ssww2,
    const double* sstt2,
    const double* dlogi,
    double* vecM,double* var);

void dlogfmarg(
    const double* x,
    const int* n,
    const int* p,
    double* var,
    const int* out,
    double* den, double* dlogi);

