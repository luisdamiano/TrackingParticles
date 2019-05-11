/**
 * @file filter.c
 * @authors Luis Damiano
 * @version 0.1
 * @details
 *
 * Sequential Importance Resampling, also known as Particle Filter.
 */

#include "main.h"

/**
 * Compute the posterior mean of the latent matrix via a Particle Filter.
 *
 * @param y The measurement vector.
 * @param nParticles The number of particles (MC samples) to use.
 * @param param The model parameters.
 * @param xMeanOut Pointer to the T x STATE_DIM matrix where the resulting
 * posterior mean matrix will be stored.
 * @param wOut Pointer to the T x nParticles matrix where the weights will be
 * stored.
 * @param essOut Pointer to the T sized vector where the resulting effective
 * sample size will be stored.
 *
 * @note This function allocates several data structures. Don't forget to call
 * `filter_free`.
 */
void filter(gsl_matrix *y, int nParticles, model_param *param,
		gsl_matrix **xMeanOut, gsl_matrix **wOut,
		gsl_vector **essOut) {
	/* Notation and indexing rules
	 *
	 * NAME INDEXING   : DESCRIPTION			(EXAMPLE  )
	 * N i = 1, ..., N : number of particles (MC samples)	(N =  1000)
	 * T k = 1, ..., T : series length                      (T = 11027)
	 * m               : measurement model vector dimension (m =     2)
	 * y[k, m]         : measurement vector
	 * w[i, k]         : weights
	 * n               : system state vector dimension      (n =     4)
	 * x[i, k, n]      : state vector
	 */

	/* Allocate error handlers */
	gsl_error_handler_t *oldHandler;
	gsl_sf_result res;
	int check;

	/* Initialize random number generator */
	const gsl_rng_type *rType;
	rType = gsl_rng_default;

	gsl_rng *r;
	gsl_rng_env_setup();
	r = gsl_rng_alloc(rType);

	/* Preallocate and initialize filtering quantities */
	int T = y->size1;
	gsl_vector_view yk, xk, xkm1, wk;
	gsl_matrix_view y1tok;
	double xkMean1, xkMean2, xkMean3, xkMean4,
		lpdf1, lpdf2, lpdf3, wkm1i, lwkm1i, wki, w0, wSum, wSumSq;

	gsl_matrix **x = (gsl_matrix **)malloc(nParticles *
							sizeof(gsl_matrix *));
	for (int i = 0; i < nParticles; i++)
		x[i] = gsl_matrix_alloc(T + 1, STATE_DIM);

	/* k = 0 (previous-to-first step) */
	/* Draw initial state -- Sarkka Eq. 7.28 */
	w0 = 1.0 / nParticles;
	wk = gsl_matrix_row(*wOut, 0);
	gsl_vector_set_all(&wk.vector, w0);
	for (int i = 0; i < nParticles; i++) {
		xk = gsl_matrix_row(x[i], 0);
		IOUT(0); IOUT(i)
		stateprior_r(r, param, &xk.vector);
		EOUT()
	}

	/* k = 1, 2, ..., T (each time step) */
	for (int k = 1; k < T + 1; k++) {
		yk = gsl_matrix_row(y, k - 1); /* Note: k - 1! */
		y1tok = gsl_matrix_submatrix(y, 0, 0, k - 1, y->size2);

		for (int i = 0; i < nParticles; i++) {
			IOUT(k);IOUT(i)

			/* Draw candidates -- Sarkka Step 1 Eq. 7.29 */
			xk = gsl_matrix_row(x[i], k);
			xkm1 = gsl_matrix_row(x[i], k - 1);

			importance_r(r, &xkm1.vector, &y1tok.matrix, param,
								&xk.vector);

			//			state_update(&xk.vector, &xkm1.vector, param);
			measurement_update(&yk.vector, &xk.vector, param);

			/* Update weights -- Sarkka Step 2 Eq. 7.30 */
			/* (1) Precompute quantities */
			measurement_lpdf(&yk.vector, &xk.vector, param, &lpdf1);
			state_lpdf(&xk.vector, &xkm1.vector, param, &lpdf2);
			importance_lpdf(&xk.vector, &xkm1.vector,
						&y1tok.matrix, param, &lpdf3);
			wkm1i = gsl_matrix_get(*wOut, k - 1, i);
			wki = 0;

			/* (2) Calculate new weight */
			oldHandler = gsl_set_error_handler_off();

			/**
			 * TODO Design a unified strategy to deal with
			 * numerical errors.
			 *
			 * Small weights produce numerical errors with both log
			 * (here) and exp (below). Currently, we deal with
			 * them independently sometimes fixing it twice.
			 */
			check = gsl_sf_log_e(wkm1i, &res);
			if (check) { /* numerical error */
				 /* Assume underflow & replace with the
				  * smallest representation of log(x) */
				lwkm1i = GSL_LOG_DBL_MIN;
#ifdef DEBUG
				printf("Numerical error gsl_sf_log_e: k % 5i, t % 5i, code % 5i, wkm1i %8.2f, lpdf1 %8.2f, lpdf2 %8.2f, lpdf3 %8.2f\n", k, i, check, wkm1i, lpdf1, lpdf2, lpdf3);
#endif
			} else {
				lwkm1i = res.val;
			}

			check = gsl_sf_exp_e(lwkm1i + lpdf1 + lpdf2 - lpdf3,
									&res);
			if (check) { /* numerical error */
				/**
				 * TODO Implement a better strategy to deal
				 * with under/overflows.
				 */
				if (check == GSL_EUNDRFLW) {
					/* Replace with the representation of
					 * the smallest positive number. */
					wki = GSL_DBL_MIN;
				} else { /* gsl_sf_exp_e only returns
						underflows or overflows */
					wki = GSL_DBL_MAX;
				}
#ifdef DEBUG
				printf("Numerical error gsl_sf_exp_e: k % 5i, t % 5i, code % 5i, wkm1i %8.2f, lpdf1 %8.2f, lpdf2 %8.2f, lpdf3 %8.2f\n", k, i, check, wkm1i, lpdf1, lpdf2, lpdf3);
#endif
			} else {
				wki = res.val;
			}

			gsl_set_error_handler(oldHandler);

			/* (3) Update weight matrix */
			gsl_matrix_set(*wOut, k, i, wki);

			DOUT(lpdf1);DOUT(lpdf2);DOUT(lpdf3);
			DOUT(wkm1i);DOUT(lwkm1i);
			DOUT(lwkm1i + lpdf1 + lpdf2 - lpdf3);
			DOUT(wki);
			IOUT(check);
			EOUT()
		} /* for each particle i */

		/* Normalize weights -- Sarkka Step 2 Eq. 7.30 */
		/* NOTE: We keep k (time step) fixed and normalize
		 * over i (particles).
		 */
		wk = gsl_matrix_row(*wOut, k);

		wSum = 0;
		for (int i = 0; i < nParticles; i++)
			wSum += gsl_vector_get(&wk.vector, i);

		gsl_vector_scale(&wk.vector, 1 / wSum);

		/* Adaptive resampling -- Sarkka Step 3 */
		/* (1) Compute effective sample size Sarkka Eq. 7.27 */

		/* TODO Implement weight squaring as a dot product
		 * for increased efficiency.
		 * Reference: https://www.gnu.org/software/gsl/manual/html_node/Level-1-GSL-BLAS-Interface.html
		 */
		wSumSq = 0;
		for (int i = 0; i < nParticles; i++)
			wSumSq += gsl_matrix_get(*wOut, k, i) *
					gsl_matrix_get(*wOut, k, i);

		gsl_vector_set(*essOut, k, 1 / wSumSq);

		/* (2) Resample */
		/* TODO Implement adaptive resampling */

		/* Compute posterior mean -- Sarkka Eq. 7.32 */
		/* TODO Implement posterior mean approximation as a matrix
		 * operation for increased efficiency. */
		xkMean1 = 0; xkMean2 = 0; xkMean3 = 0; xkMean4 = 0;
		for (int i = 0; i < nParticles; i++) {
			wki = gsl_matrix_get(*wOut, k, i);
			xk = gsl_matrix_row(x[i], k);
			xkMean1 += gsl_vector_get(&xk.vector, 0) * wki;
			xkMean2 += gsl_vector_get(&xk.vector, 1) * wki;
			xkMean3 += gsl_vector_get(&xk.vector, 2) * wki;
			xkMean4 += gsl_vector_get(&xk.vector, 3) * wki;
		}

		gsl_matrix_set(*xMeanOut, k, 0, xkMean1);
		gsl_matrix_set(*xMeanOut, k, 1, xkMean2);
		gsl_matrix_set(*xMeanOut, k, 2, xkMean3);
		gsl_matrix_set(*xMeanOut, k, 3, xkMean4);

#ifdef DEBUG
		printf("k = % 5i, total wSum %0.8f, ESS: % 10.6f \t \t %0.8f\t%0.8f\t%0.8f\t%0.8f\n", k, wSum, 1 / wSumSq, xkMean1, xkMean2, xkMean3, xkMean4);
#endif
	} /* for each time step k */

	/* Cleanup */
	/* NOTE: Don't free xMean, w, ess -- pointers to these are returned. */
	for (int i = 0; i < nParticles; i++)
		gsl_matrix_free(x[i]);
	free(x);
	gsl_rng_free (r);
}
