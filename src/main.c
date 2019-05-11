/**
 * @file main.c
 * @authors Luis Damiano
 * @version 0.1
 * @details
 *
 * DESCRIPTION
 *
 * Usage:
 *	./??
 *
 * Compile:
 *  ./??
 *
 * Troubleshooting:
 *  ????
 */

#include "main.h"

/* Files */
#define MEASUREMENT_FILE_IN "../R/data/measurements.txt"
#define ESS_FILE_OUT "essOut.txt"
#define WEIGHTS_FILE_OUT "wOut.txt"
#define STATEMEAN_FILE_OUT "xMeanOut.txt"
#define BASELINE_FILE_OUT "baselineOut.txt"

/* Measurement model constants */
#define DT 1.0 /* Keep it double, will you? */
#define LOCATION_1_X -93.2494663765932f
#define LOCATION_1_Y  41.5563518606521f
#define LOCATION_2_X -93.2475338232000f
#define LOCATION_2_Y  41.5576632356000f
#define MEASUREMENT_ERROR_1 0.01;

/* State model */
#define STATE_DIFFUSION_1 0.0005
#define STATE_DIFFUSION_2 0.0005

/* State prior */
#define STATEPRIOR_MU_X -93.24952047
#define STATEPRIOR_MU_Y 41.55575337
#define STATEPRIOR_L_00 5.0E-09
#define STATEPRIOR_L_11 3.5E-08
#define STATEPRIOR_L_22 5.0E-04
#define STATEPRIOR_L_33 5.0E-04

/* Importance distribution */
#define IMPORTANCE_L_00 3 * 5.00E-10
#define IMPORTANCE_L_11 3 * 1.75E-08
#define IMPORTANCE_L_22 3 * 5.00E-05
#define IMPORTANCE_L_33 3 * 5.00E-05

/* Particle filter constants */
#define NPARTICLES 100

int main(int argc, char** argv)
{
	INITOUT()

	/* Read data */
	gsl_matrix* y = (gsl_matrix*)malloc(sizeof(gsl_matrix));
	load_data(MEASUREMENT_FILE_IN, &y);
	int T = y->size1;

	/* Compute noiseless solution */
	gsl_vector* location1 = gsl_vector_alloc(MEASUREMENT_DIM);
	gsl_vector* location2 = gsl_vector_alloc(MEASUREMENT_DIM);

	gsl_vector_set(location1, 0, LOCATION_1_X);
	gsl_vector_set(location1, 1, LOCATION_1_Y);
	gsl_vector_set(location2, 0, LOCATION_2_X);
	gsl_vector_set(location2, 1, LOCATION_2_Y);

	gsl_matrix *baseline = gsl_matrix_alloc(T, MEASUREMENT_DIM);
	noiseless(y, location1, location2, baseline);

	/* Initialize model */
	model_param param;
	param.baseline = baseline;
	param.dt = DT;
	param.l1x = LOCATION_1_X;
	param.l1y = LOCATION_1_Y;
	param.l2x = LOCATION_2_X;
	param.l2y = LOCATION_2_Y;
	param.sr = MEASUREMENT_ERROR_1;

	param.q1 = STATE_DIFFUSION_1;
	param.q2 = STATE_DIFFUSION_2;

	param.statepriorMuX = STATEPRIOR_MU_X;
	param.statepriorMuY = STATEPRIOR_MU_Y;
	param.statepriorL00 = STATEPRIOR_L_00;
	param.statepriorL11 = STATEPRIOR_L_11;
	param.statepriorL22 = STATEPRIOR_L_22;
	param.statepriorL33 = STATEPRIOR_L_33;

	param.importanceL00 = IMPORTANCE_L_00;
	param.importanceL11 = IMPORTANCE_L_11;
	param.importanceL22 = IMPORTANCE_L_22;
	param.importanceL33 = IMPORTANCE_L_33;

	importance_init(&param);
	state_init(&param);
	measurement_init(&param);

	/* Run particle filter */
	gsl_matrix *xMeanOut = gsl_matrix_alloc(T + 1, STATE_DIM);
	gsl_matrix *wOut = gsl_matrix_alloc(T + 1, NPARTICLES);
	gsl_vector* essOut = gsl_vector_alloc(T + 1);

	filter(y, NPARTICLES, &param, &xMeanOut, &wOut, &essOut);

	/* Write results to disk */
	GSL_MAT_TO_CSV(baseline, BASELINE_FILE_OUT);
	GSL_MAT_TO_CSV(xMeanOut, STATEMEAN_FILE_OUT);
	GSL_MAT_TO_CSV(wOut, WEIGHTS_FILE_OUT);
	GSL_VEC_TO_CSV(essOut, ESS_FILE_OUT);

	/* Clean up */
	gsl_matrix_free(xMeanOut);
	gsl_matrix_free(wOut);
	gsl_vector_free(essOut);

	importance_free(&param);
	state_free(&param);
	measurement_free(&param);

	gsl_vector_free(location2);
	gsl_vector_free(location1);
	gsl_matrix_free(baseline);
	free(y);

	/* Say goodbye */
	EXITOUT()
	return EXIT_SUCCESS;
}
