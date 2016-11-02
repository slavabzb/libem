#ifndef EM_H
#define EM_H

#include <matrix/matrix.h>

int em_optimize(mpfr_t* const fx,
		struct mtx* const x,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c,
		double eps);

#endif /* EM_H */

