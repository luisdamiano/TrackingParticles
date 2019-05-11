/**
 * @file Rfilter.c
 * @authors Luis Damiano
 * @version 0.1
 * @details R wrapper for the (awesome) Particle Filter.
 */

#include "main.h"

void Rfilter(double *Ry1, double *Ry2, int *RT,
		double *LOCATION_1_X, double *LOCATION_1_Y,
		double *LOCATION_2_X, double *LOCATION_2_Y,
		double *DT,
		double *MEASUREMENT_ERROR_1,
		double *STATE_DIFFUSION_1, double *STATE_DIFFUSION_2,
		double *STATEPRIOR_MU_X, double *STATEPRIOR_MU_Y,
		double *STATEPRIOR_L_00, double *STATEPRIOR_L_11,
		double *STATEPRIOR_L_22, double *STATEPRIOR_L_33,
		double *IMPORTANCE_L_00, double *IMPORTANCE_L_11,
		double *IMPORTANCE_L_22, double *IMPORTANCE_L_33,
		int* NPARTICLES,
		double *noiselessOut,
		double *RxMeanOut, double *RwOut, double *RessOut);

void Rfilter(double *Ry1, double *Ry2, int *RT,
		double *LOCATION_1_X, double *LOCATION_1_Y,
		double *LOCATION_2_X, double *LOCATION_2_Y,
		double *DT,
		double *MEASUREMENT_ERROR_1,
		double *STATE_DIFFUSION_1, double *STATE_DIFFUSION_2,
		double *STATEPRIOR_MU_X, double *STATEPRIOR_MU_Y,
		double *STATEPRIOR_L_00, double *STATEPRIOR_L_11,
		double *STATEPRIOR_L_22, double *STATEPRIOR_L_33,
		double *IMPORTANCE_L_00, double *IMPORTANCE_L_11,
		double *IMPORTANCE_L_22, double *IMPORTANCE_L_33,
		int* NPARTICLES,
		double *RnoiselessOut,
		double *RxMeanOut, double *RwOut, double *RessOut) {

	/* Read data from R*/
	gsl_matrix *y = gsl_matrix_alloc(*RT, MEASUREMENT_DIM);
	int T = *RT;

	/* RECALL: R is col-major order while GSL is row-major order. */
	for (int i = 0; i < T; i++) {
		gsl_matrix_set(y, i, 0, Ry1[i]);
		gsl_matrix_set(y, i, 1, Ry2[i]);
	}

	gsl_vector* location1 = gsl_vector_alloc(MEASUREMENT_DIM);
	gsl_vector* location2 = gsl_vector_alloc(MEASUREMENT_DIM);

	gsl_vector_set(location1, 0, *LOCATION_1_X);
	gsl_vector_set(location1, 1, *LOCATION_1_Y);
	gsl_vector_set(location2, 0, *LOCATION_2_X);
	gsl_vector_set(location2, 1, *LOCATION_2_Y);

	gsl_matrix *baseline = gsl_matrix_alloc(T, MEASUREMENT_DIM);
	noiseless(y, location1, location2, baseline);

	/* Initialize model */
	model_param param;
	param.baseline = baseline;
	param.dt = *DT;
	param.l1x = *LOCATION_1_X;
	param.l1y = *LOCATION_1_Y;
	param.l2x = *LOCATION_2_X;
	param.l2y = *LOCATION_2_Y;
	param.sr = *MEASUREMENT_ERROR_1;

	param.q1 = *STATE_DIFFUSION_1;
	param.q2 = *STATE_DIFFUSION_2;

	param.statepriorMuX = *STATEPRIOR_MU_X;
	param.statepriorMuY = *STATEPRIOR_MU_Y;
	param.statepriorL00 = *STATEPRIOR_L_00;
	param.statepriorL11 = *STATEPRIOR_L_11;
	param.statepriorL22 = *STATEPRIOR_L_22;
	param.statepriorL33 = *STATEPRIOR_L_33;

	param.importanceL00 = *IMPORTANCE_L_00;
	param.importanceL11 = *IMPORTANCE_L_11;
	param.importanceL22 = *IMPORTANCE_L_22;
	param.importanceL33 = *IMPORTANCE_L_33;

	importance_init(&param);
	state_init(&param);
	measurement_init(&param);

	/* Run particle filter */
	gsl_matrix *xMeanOut = gsl_matrix_alloc(T + 1, STATE_DIM);
	gsl_matrix *wOut = gsl_matrix_alloc(T + 1, *NPARTICLES);
	gsl_vector* essOut = gsl_vector_alloc(T + 1);

	filter(y, *NPARTICLES, &param, &xMeanOut, &wOut, &essOut);

	/* Write results to R */
	for (int i = 0; i < T; i++)
		for (int j = 0; j < MEASUREMENT_DIM; j++)
			RnoiselessOut[i + j * T] =
					gsl_matrix_get(baseline, i, j);

	for (int i = 0; i < T + 1; i++)
		for (int j = 0; j < STATE_DIM; j++)
			RxMeanOut[i + j * (T + 1)] =
					gsl_matrix_get(xMeanOut, i, j);

	for (int i = 0; i < T + 1; i++)
		for (int j = 0; j < *NPARTICLES; j++)
			RwOut[i + j * (T + 1)] = gsl_matrix_get(wOut, i, j);

	for (int i = 0; i < T + 1; i++)
		RessOut[i] = gsl_vector_get(essOut, i);

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

	// Say goodbye?
}
