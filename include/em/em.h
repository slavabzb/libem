#ifndef EM_H
#define EM_H

#include <matrix/matrix.h>

/**
 * Solve the optimization problem.
 * 
 * Uses the ellipsoid method algorithm to find out an optimum
 * of the given optimization problem.
 * 
 * @param niters number of iterations spent by algorithm,
 * @param fx minimum value of the target function,
 * @param x point where the target function has it's minimum
 * @param a matrix of constraint coefficients
 * @param b non-negative vector
 * @param c objective function coefficients
 * @param eps precision
 * @return 0 on success, non-zero otherwise
 */
int em_optimize(mpz_t* const niters,
		mpfr_t* const fx,
		struct mtx* const x,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c,
		double eps);

#endif /* EM_H */

