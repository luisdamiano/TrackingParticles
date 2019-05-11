/**
 * @file filter.h
 * @authors Luis Damiano
 * @version 0.1
 * @details
 *
 * Sequential Importance Resampling, also known as Particle Filter.
 */

#ifndef C_FILTER_H_
#define C_FILTER_H_

void filter(gsl_matrix *y, int nParticles, model_param *param,
		gsl_matrix **xMeanOut, gsl_matrix **wOut,
		gsl_vector **essOut);
void filter_free(gsl_matrix *xMeanOut, gsl_matrix *wOut, gsl_vector *essOut);

#endif /* C_FILTER_H_ */
