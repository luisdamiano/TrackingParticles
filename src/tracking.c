/**
 * @file tracking.c
 * @authors Luis Damiano
 * @version 0.1
 * @details
 *
 * Density function and random number generation for the densities involved in
 * the particle filter for the bearing-only tracking problem.
 *
 * Naming convention for mathematical variables:
 *   `y` refers to measurements and `x` to states.
 *
 *   C NAME	MATH NOTATION	TYPE	DESCRIPTION
 *   x		x		Matrix 	From start to end.
 *
 *   xk		x_{k}  		Vector 	For the current time step.
 *   xkm1	x_{k-1}		Vector 	For the previous time step.
 *
 *   x1tok	x_{1:k}		Matrix 	From start up to k included.
 *   x1tokm1	x_{1:(k-1)}	Matrix 	From start up to k-1 included.
 *
 * Design convention for functions:
 *
 *   Expect gsl_vector* and gsl_matrix* as arguments -- don't ask for views.
 *   Write results directly -- avoid returns to avoid excessive copying.
 */

#include "main.h"

/** FIRST PART: RANDOM GENERATION AND DENSITY FUNCTIONS --------------------- */

void stateprior_r(const gsl_rng *r, model_param *param, gsl_vector *xOut) {
	gsl_ran_multivariate_gaussian(r, param->statepriorMu,
			param->statepriorL, xOut);
	VOUT(param->statepriorMu)
	MOUT(param->statepriorL)
	VOUT(xOut)
}

void importance_r(const gsl_rng *r, gsl_vector *xkm1, gsl_matrix *y1tok,
		model_param *param, gsl_vector *xOut) {
	gsl_vector_view mu = gsl_matrix_row(param->baseline, y1tok->size1);
	double padded[] = {
			gsl_vector_get(&mu.vector, 0),
			gsl_vector_get(&mu.vector, 1),
			0,
			0
	};

	gsl_vector_view center = gsl_vector_view_array(padded, 4);
	gsl_ran_multivariate_gaussian(r, &center.vector, param->importanceL,
			xOut);
	VOUT(&center.vector)
	MOUT(param->importanceL)
	VOUT(xOut)
}

void importance_lpdf(gsl_vector *xk, gsl_vector *xkm1, gsl_matrix *y1tok,
		model_param *param, double *lpdf) {
	gsl_ran_multivariate_gaussian_log_pdf(xk, xkm1, param->importanceL,
			lpdf, param->stateWork);
	VOUT(xk)
	VOUT(xkm1)
	MOUT(param->importanceL)
	DOUT(*lpdf)
}

void measurement_lpdf(gsl_vector *yk, gsl_vector *xk, model_param *param,
		double *lpdf) {
	gsl_ran_multivariate_gaussian_log_pdf(yk, param->measurementMu,
			param->measurementL, lpdf, param->measurementWork);
	VOUT(yk)
	VOUT(param->measurementMu)
	MOUT(param->measurementL)
	DOUT(*lpdf)
}

void state_lpdf(gsl_vector *xk, gsl_vector *xkm1, model_param *param,
		double *lpdf) {
	gsl_ran_multivariate_gaussian_log_pdf(xk, param->stateMu,
			param->stateL, lpdf, param->stateWork);
	VOUT(xk)
	VOUT(param->stateMu)
	MOUT(param->stateL)
	DOUT(*lpdf)
}

/** SECOND PART: PARAMETER UPDATING FUNCTIONS ------------------------------- */

void importance_init(model_param *param) {
	/* Allocate & populate importance covariance matrix */
	param->importanceL = gsl_matrix_calloc(STATE_DIM, STATE_DIM);

	gsl_matrix_set(param->importanceL, 0, 0, param->importanceL00);
	gsl_matrix_set(param->importanceL, 1, 1, param->importanceL11);
	gsl_matrix_set(param->importanceL, 2, 2, param->importanceL22);
	gsl_matrix_set(param->importanceL, 3, 3, param->importanceL33);

	gsl_linalg_cholesky_decomp(param->importanceL);
}

void importance_free(model_param *param) {
	gsl_matrix_free(param->importanceL);
}

void measurement_init(model_param *param) {
	/* Allocate mean vector */
	param->measurementMu = gsl_vector_calloc(MEASUREMENT_DIM);

	/* Allocate & populate covariance matrix */
	param->measurementL = gsl_matrix_alloc(MEASUREMENT_DIM,
							MEASUREMENT_DIM);

	gsl_matrix_set(param->measurementL, 0, 0, param->sr);
	gsl_matrix_set(param->measurementL, 0, 1, 0);
	gsl_matrix_set(param->measurementL, 1, 0, 0);
	gsl_matrix_set(param->measurementL, 1, 1, param->sr);
	gsl_linalg_cholesky_decomp(param->measurementL);

	/* Allocate work vector */
	param->measurementWork = gsl_vector_alloc(MEASUREMENT_DIM);
}

void measurement_free(model_param *param) {
	gsl_vector_free(param->measurementWork);
	gsl_matrix_free(param->measurementL);
	gsl_vector_free(param->measurementMu);
}

void measurement_update(gsl_vector *yk, gsl_vector *xk, model_param *param) {
	/* RECALL:
	 * Observation vector = (angle1, angle2)
	 * State vector       = (x-coord, y-coord, x-velocity, y-velocity)
	 */

	/* Update mean vector: the accuracy transformation was verified. */
	double xkXCoord = gsl_vector_get(xk, 0); /**< x-coord at k from
								state model */
	double xkYCoord = gsl_vector_get(xk, 1); /**< y-coord at k from
								state model */
	double mean1 = atan2(xkYCoord - param->l1y, xkXCoord - param->l1x);
	double mean2 = atan2(xkYCoord - param->l2y, xkXCoord - param->l2x);

	gsl_vector_set(param->measurementMu, 0, mean1);
	gsl_vector_set(param->measurementMu, 1, mean2);

	/* Update covariance matrix */
	/* We make no update here as the covariance matrix is fixed
	 * in our case. Yet, we leave this comment as a "placeholder"
	 * remainder in case the model is to be extended in the future.
	 */
}

void state_init(model_param *param) {
	/* Allocate mean vector and set to baseline */
	param->stateMu = gsl_vector_alloc(STATE_DIM);
	gsl_vector_view baseline1 = gsl_matrix_row(param->baseline, 0);

	gsl_vector_set(param->stateMu, 0,
			gsl_vector_get(&baseline1.vector, 0));
	gsl_vector_set(param->stateMu, 1,
			gsl_vector_get(&baseline1.vector, 1));
	gsl_vector_set(param->stateMu, 2, 0);
	gsl_vector_set(param->stateMu, 3, 0);

	/* Allocate & populate covariance matrix */
	param->stateL = gsl_matrix_alloc(STATE_DIM, STATE_DIM);

	double dt3 = param->dt * param->dt * param->dt / 3;
	double q1dt3 = param->q1 * dt3;
	double q2dt3 = param->q2 * dt3;

	double dt2 = param->dt * param->dt / 2;
	double q1dt2 = param->q1 * dt2;
	double q2dt2 = param->q2 * dt2;

	/* Careful here -- getting the Q matrix right is super tricky */
	gsl_matrix_set_zero(param->stateL);
	gsl_matrix_set(param->stateL, 0, 0, q1dt3);
	gsl_matrix_set(param->stateL, 2, 0, q1dt2);

	gsl_matrix_set(param->stateL, 1, 1, q2dt3);
	gsl_matrix_set(param->stateL, 3, 1, q2dt2);

	gsl_matrix_set(param->stateL, 0, 2, q1dt2);
	gsl_matrix_set(param->stateL, 2, 2, param->q1 * param->dt);

	gsl_matrix_set(param->stateL, 1, 3, q2dt2);
	gsl_matrix_set(param->stateL, 3, 3, param->q2 * param->dt);

	gsl_linalg_cholesky_decomp(param->stateL);

	/* Allocate & populate transition matrix */
	param->stateTransition = gsl_matrix_alloc(STATE_DIM, STATE_DIM);

	gsl_matrix_set_identity(param->stateTransition);
	gsl_matrix_set(param->stateTransition, 0, 2, param->dt);
	gsl_matrix_set(param->stateTransition, 1, 3, param->dt);

	/* Allocate & populate state prior mean vector */
	param->statepriorMu = gsl_vector_calloc(STATE_DIM);
	gsl_vector_set(param->statepriorMu, 0, param->statepriorMuX);
	gsl_vector_set(param->statepriorMu, 1, param->statepriorMuY);

	/* Allocate & populate state prior covariance matrix */
	param->statepriorL = gsl_matrix_calloc(STATE_DIM, STATE_DIM);
	gsl_matrix_set(param->statepriorL, 0, 0, param->statepriorL00);
	gsl_matrix_set(param->statepriorL, 1, 1, param->statepriorL11);
	gsl_matrix_set(param->statepriorL, 2, 2, param->statepriorL22);
	gsl_matrix_set(param->statepriorL, 3, 3, param->statepriorL33);

	/* Allocate work vector */
	param->stateWork = gsl_vector_alloc(STATE_DIM);
}

void state_free(model_param *param) {
	gsl_vector_free(param->stateWork);
	gsl_matrix_free(param->statepriorL);
	gsl_vector_free(param->statepriorMu);
	gsl_matrix_free(param->stateTransition);
	gsl_matrix_free(param->stateL);
	gsl_vector_free(param->stateMu);
}

#if FALSE
void state_update(gsl_vector *xk, gsl_vector *xkm1, model_param *param) {
	/* Update mean vector
	 *
	 * NOTE: Transition equation (p. 93) is a matrix-vector operation.
	 * Level 2 Matrix-vector operations
	 * D	Double precision
	 * GE	General matrix
	 * MV	Matrix-vector product A*x
	 */
	gsl_blas_dgemv(CblasNoTrans, 1, param->stateTransition, xkm1, 0, xk);

	/* Update covariance matrix */
	/* We make no update here as the covariance matrix is fixed
	 * in our case. Yet, we leave this comment as a "placeholder"
	 * remainder in case the model is to be extended in the future.
	 */
}
#endif
