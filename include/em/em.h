#ifndef EM_H
#define EM_H

#include <matrix/matrix.h>

/**
 * Solve an optimization problem.
 * 
 * Uses the ellipsoid method algorithm to find out an optimum
 * of the given optimization problem.
 * 
 * @param fx minimum value of the target function,
 * @param x point where the target function has it's minimum
 * @param prec answer precision
 * @param niters number of iterations spent by algorithm,
 * @param a matrix of constraint coefficients
 * @param b non-negative vector
 * @param c objective function coefficients
 * @param eps accuracy
 * @return 0 on success, non-zero otherwise
 */
int em_optimize(mpfr_t* const fx,
		struct mtx* const x,
		mpfr_prec_t* const prec,
		mpz_t* const niters,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c,
		double eps);

#endif /* EM_H */

