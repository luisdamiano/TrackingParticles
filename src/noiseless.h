/**
 * @file noiseless.h
 * @authors Luis Damiano
 * @version 0.1
 * @details
 *
 * Header for noiseless solution.
 */

#ifndef C_NOISELESS_H_
#define C_NOISELESS_H_

void noiseless(gsl_matrix* angles, gsl_vector* location1,
		gsl_vector* location2, gsl_matrix* solutionOut);

#endif /* C_INTERFACE_H_ */
