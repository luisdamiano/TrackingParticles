/**
 * @file main.h
 * @authors Luis Damiano
 * @version 0.1
 * @details
 *
 * General particleawe settings & includes. For settings, see first
 * two blocks of defines.
 */

#ifndef C_MAIN_H_
#define C_MAIN_H_

/* particleawe settings */
/* #define DEBUG */
/* #define CSVOUT */

/* GLS Settings */
#define GSL_RANGE_CHECK_OFF
#define HAVE_INLINE

/* Awesome macros */
#define GSL_VEC_TO_CSV(x, filename) do {\
	FILE *fp = fopen((filename), "w");\
	 	if (fp == NULL)\
		fatal("couldn't create file to store the resulting vector");\
	for (int i = 0; i < ((x))->size; i++) {\
		fprintf(fp, "% 19.17f\n", gsl_vector_get((x), i));\
	}\
	fclose(fp);\
} while(0)

#define GSL_MAT_TO_CSV(x, filename) do {\
	FILE *fp = fopen((filename), "w");\
	 	if (fp == NULL)\
		fatal("couldn't create file to store the resulting matrix");\
	for (int i = 0; i < ((x))->size1; i++) {\
		for (int j = 0; j < ((x))->size2; j++)\
			fprintf(fp, "% 19.17f,", gsl_matrix_get((x), i, j));\
		fprintf(fp, "\n");\
	}\
	fclose(fp);\
} while(0)

#ifdef CSVOUT

FILE *outFile;
#define HEADEROUT "k,i,importanceMu,importanceMu,importanceMu,importanceMu,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importance_r,importance_r,importance_r,importance_r,yk,yk,measurementMu,measurementMu,measurementL,measurementL,measurementL,measurementL,measurement_lpdf,xk,xk,xk,xk,stateMu,stateMu,stateMu,stateMu,stateL,stateL,stateL,stateL,stateL,stateL,stateL,stateL,stateL,stateL,stateL,stateL,stateL,stateL,stateL,stateL,state_lpdf,xk,xk,xk,xk,xkm1,xkm1,xkm1,xkm1,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importanceL,importance_lpdf,lpdf1,lpdf2,lpdf3,wkm1i,logwkm1i,lognewweight,newweight,errno\n"
#define INITOUT() outFile = fopen("filter.csv", "w");fprintf(outFile, HEADEROUT);
#define EXITOUT() fclose(outFile);
#define EOUT() fprintf(outFile, "\n"); // End of line
#define SOUT(x) fprintf(outFile, "%s,", x);
#define IOUT(x) fprintf(outFile, "%i,", x);
#define DOUT(x) fprintf(outFile, "%.17g,", x); // https://stackoverflow.com/a/21162120/2860744
#define VOUT(x) for(int tmpint = 0; tmpint < (x)->size; tmpint++) { DOUT(gsl_vector_get((x), tmpint)) }
#define MOUT(x) for(int tmpint1 = 0; tmpint1 < (x)->size1; tmpint1++) { for(int tmpint2 = 0; tmpint2 < (x)->size2; tmpint2++) { DOUT(gsl_matrix_get((x), tmpint1, tmpint2)) } }

#else

#define EOUT() do {} while(0);
#define INITOUT() do {} while(0);
#define EXITOUT() do {} while(0);
#define SOUT(x) do {} while(0);
#define IOUT(x) do {} while(0);
#define DOUT(x) do {} while(0);
#define VOUT(x) do {} while(0);
#define MOUT(x) do {} while(0);

#endif

/* Includes */

#include <errno.h>
#include <math.h>
#include <stdlib.h> /* malloc, free */
#include <string.h> /* strcpy */
#include <unistd.h> /* getopt */

#include <gsl/gsl_blas.h>
#include <gsl/gsl_blas_types.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h> /* gsl_linalg_cholesky_decomp */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_vector.h>

#include "interface.h"
#include "noiseless.h"
#include "load.h"
#include "model.h"
#include "tracking.h"
#include "filter.h"

#endif /* C_MAIN_H_ */
