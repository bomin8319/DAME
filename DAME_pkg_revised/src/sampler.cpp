#include <cmath>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


void R_init_markovchain(DllInfo* info) {
	R_registerRoutines(info, NULL, NULL, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);	
}

extern void dsyev_( char *jobz, char *uplo, int *n, double *a, int *lda,
                   double *w, double *work, int *lwork, int *info );

