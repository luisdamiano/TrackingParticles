/**
 * @file tracking.h
 * @authors Luis Damiano
 * @version 0.1
 * @details
 *
 * Specifics to the bearing-only tracking problem.
 */

#ifndef C_TRACKING_H_
#define C_TRACKING_H_

/* Hard constants (can't be changed without modifying the source code first) */
#define MEASUREMENT_DIM 2 /* int */
#define STATE_DIM 4 /* int */

void importance_init(model_param *param);
void importance_free(model_param *param);
void importance_r(const gsl_rng *r, gsl_vector *xkm1, gsl_matrix *y1tok,
		model_param *param, gsl_vector *xOut);
void importance_lpdf(gsl_vector *x, gsl_vector *xkm1, gsl_matrix *y1tok,
		model_param *param, double *lpdf);

void state_init(model_param *param);
void state_update(gsl_vector *xk, gsl_vector *xkm1, model_param *param);
void state_free(model_param *param);
void state_lpdf(gsl_vector *xk, gsl_vector *xkm1, model_param *param,
		double *lpdf);
void stateprior_r(const gsl_rng *r, model_param *param, gsl_vector *xOut);

void measurement_init(model_param *param);
void measurement_free(model_param *param);
void measurement_update(gsl_vector *yk, gsl_vector *xk, model_param *param);
void measurement_lpdf(gsl_vector *yk, gsl_vector *xk, model_param *param,
		double *lpdf);

#endif /* C_TRACKING_H_ */
