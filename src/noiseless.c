/**
 * @file noiseless.c
 * @authors Luis Damiano
 * @version 0.1
 * @details
 *
 * Implementation of the noiseless solution of the bearing-only tracking
 * problem.
 */

#include "main.h"
/**
 * Compute the noiseless solution of the bearing-only tracking problem.
 *
 * @param angles A 2-column matrix with the sensor measurements.
 * @param location1 A 2-element vector with the coordinates (x, y)
 * of the first sensor.
 * @param location2 A 2-element vector with the coordinates (x, y)
 * of the second sensor.
 * @param solutionOut Pointer to the 2-column matrix where the solutions
 * will be stored.
 */
void noiseless(gsl_matrix* angles, gsl_vector* location1,
		gsl_vector* location2, gsl_matrix* solutionOut) {
	/* Recover locations */
	double l1x, l1y;
	l1x = gsl_vector_get(location1, 0);
	l1y = gsl_vector_get(location1, 1);

	/* Allocate work matrices */
	double a1, a2;
	gsl_matrix* derivatives = gsl_matrix_alloc(MEASUREMENT_DIM,
			MEASUREMENT_DIM);
	gsl_vector* differences = gsl_vector_alloc(MEASUREMENT_DIM);
	gsl_vector* tau = gsl_vector_alloc(MEASUREMENT_DIM); /**< QR factors */
	gsl_vector* coefficients = gsl_vector_alloc(MEASUREMENT_DIM); /**< Solution x* of
	linear system `derivatives * x = differences` */

	/* Populate differences vector */
	gsl_vector_set(differences, 0, gsl_vector_get(location2, 0) - l1x);
	gsl_vector_set(differences, 1, gsl_vector_get(location2, 1) - l1y);

	for (int k = 0; k < angles->size1; k++) { /* for each time step k */
		/* Populate derivatives matrix
		 * RECALL: gsl_matrix are stored in column order */
		a1 = gsl_matrix_get(angles, k, 0);
		a2 = gsl_matrix_get(angles, k, 1);
		double dx1 = cos(a1);
		double dy1 = sin(a1);

		gsl_matrix_set(derivatives, 0, 0, dx1);
		gsl_matrix_set(derivatives, 1, 0, dy1);
		gsl_matrix_set(derivatives, 0, 1, cos(a2));
		gsl_matrix_set(derivatives, 1, 1, sin(a2));

		/* Solve linear system */
		gsl_linalg_QR_decomp(derivatives, tau);
		gsl_linalg_QR_solve(derivatives, tau, differences,
								coefficients);

		/* Write solutions */
		double term = gsl_vector_get(coefficients, 0);
		gsl_matrix_set(solutionOut, k, 0, l1x + dx1 * term);
		gsl_matrix_set(solutionOut, k, 1, l1y + dy1 * term);
	}

	/* Clean up */
	gsl_vector_free(coefficients);
	gsl_vector_free(tau);
	gsl_vector_free(differences);
	gsl_matrix_free(derivatives);
}
