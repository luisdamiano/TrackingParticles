/**
 * @file load.c
 * @authors Luis Damiano
 * @version 0.1
 * @details
 *
 * Data reading functions.
 */

#include "main.h"

/**
 * Read measurements from file.
 *
 * @param filename Path to the file with measurements.
 * @param y Pointer where the measurement matrix will be stored.
 */
void load_data(char *filename, gsl_matrix **y)
{
	/* Read file */
	FILE *pFile = fopen(filename, "r");

	if (pFile == NULL)
		fatal("cannot open file, is it accessible?");

	/* Determine length */
	int nrow = 0;
	char c;

	for (c = getc(pFile); c != EOF; c = getc(pFile))
		if (c == '\n')
			nrow++;

	/* Allocate rows */
	*y = gsl_matrix_alloc(nrow, MEASUREMENT_DIM);

	/* Populate rows */
	fseek(pFile, 0L, SEEK_SET);

	double xTmp, yTmp;
	for (int row = 0; row < nrow; row++) {
			int res = fscanf(pFile, "%lf %lf\n", &xTmp, &yTmp);
			gsl_matrix_set(*y, row, 0, xTmp);
			gsl_matrix_set(*y, row, 1, yTmp);
	}

	fclose(pFile);
}
