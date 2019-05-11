/**
 * @file interface.c
 * @authors Luis Damiano
 * @version 0.1
 * @details
 *
 * Printing functions.
 */

#include "main.h"

/**
 * Print a fatal error message and exit with failure.
 * @param message Error message.
 */
void fatal(char *message)
{
	fprintf(stderr, "FATAL: %s.\n", message);
	fprintf(stderr,
			"\nUsage:\n	./particle\n");

	exit(EXIT_FAILURE);
}

/**
 * Print a warning message.
 * @param message Warning message.
 */
void debug(char *message)
{
	fprintf(stderr, "DEBUG: %s.\n", message);
}

/**
 * Print a warning message.
 * @param message Warning message.
 */
void warning(char *message)
{
	fprintf(stderr, "WARNING: %s.\n", message);
}
