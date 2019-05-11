/**
 * @file model.h
 * @authors Luis Damiano
 * @version 0.1
 * @details
 *
 * Struct holding the model parameters.
 */

#ifndef C_MODEL_H_
#define C_MODEL_H_

typedef struct model_parameters {
	/* Model constants */
	double l1x, l1y, l2x, l2y, dt;

	/* State prior distributions */
	gsl_vector *statepriorMu; /**< Location for initial state prior */
	gsl_matrix *statepriorL; /**< Cholesky factor for initial state prior */
	double statepriorMuX, statepriorMuY;
	double statepriorL00, statepriorL11, statepriorL22, statepriorL33;

	/* Importance distribution parameters */
	gsl_matrix *importanceL; /**< Cholesky factor for importance pdf */
	gsl_matrix *baseline; /**< Noiseless approximation used by the
							importance pdf */
	double importanceL00, importanceL11, importanceL22, importanceL33;

	/* State model parameters */
	double q1; /**< Diffusion coefficient for transition in x-coords */
	double q2; /**< Diffusion coefficient for transition in y-coords */
	gsl_vector *stateMu; /**< Location for state model */
	gsl_matrix *stateL; /**< Cholesky factor for state model */
	gsl_matrix *stateTransition; /**< Transition matrix for state model */
	gsl_vector *stateWork; /**< Work variable vector of size STATE_DIM */

	/* Measurement model parameters */
	double sr; /**< Standard deviation for measurement error */
	gsl_vector *measurementMu; /**< Location for measurement model */
	gsl_matrix *measurementL; /**< Cholesky factor for measurement model */
	gsl_vector *measurementWork; /**< Work variable vector of size
							MEASUREMENT_DIM */
} model_param;

#endif /* C_MODEL_H_ */
