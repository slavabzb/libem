#ifndef EM_H
#define EM_H

#include <matrix/matrix.h>

int em_oprimize(mpfr_t* const fx,
		struct mtx* const x,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c);

#endif /* EM_H */

