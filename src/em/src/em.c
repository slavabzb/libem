#include <math.h>
#include <stdlib.h>
#include <em/em.h>
#include <errno.h>

int em_set_maxiter(size_t* const maxiter,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c,
		double const R2,
		double const eps);

int em_set_size(size_t* const size,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c);

int em_accum_size(size_t* const size,
		struct mtx const m);

int em_check_dimensions(struct mtx* const x,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c);

int em_optimize(mpfr_t* const fx,
		struct mtx* const x,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c,
		double eps)
{
	if (em_check_dimensions(x, a, b, c))
		return -1;
	
	size_t size;
	em_set_size(&size, a, b, c);
	
	mpfr_prec_t const prec = 2 * log2(a.ncols) + 2 * size;
	if (errno == EDOM || errno == ERANGE)
		return -1;
	
	mpfr_set_default_prec(prec);
	
	struct mtx y;
	if (mtx_init(&y, a.nrows, 1, prec))
		return -1;
	
	if (mtx_fill(y, 0.f, 0.f))
		return -1;
	
	struct mtx H;
	if (mtx_init(&H, a.ncols, a.ncols, prec))
		return -1;
	
	double const R2 = pow(a.ncols, 2) * pow(2, 2 * size);
	
	if (mtx_fill(H, 0, R2))
		return -1;
	
	size_t iter;
	size_t maxiter;
	if (em_set_maxiter(&maxiter, a, b, c, R2, eps))
		return -1;
	
	for (iter = 0; iter < maxiter; ++iter)
	{
		
	}
	
	return 0;
}

int em_set_maxiter(size_t* const maxiter,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c,
		double const R2,
		double const eps)
{
	mpfr_t tmp0;
	mpfr_init_set_ui(tmp0, 0, MPFR_RNDD);
	
	mpfr_t tmp1;
	mpfr_init(tmp1);
	
	int i;
	for (i = 0; i < c.ncols; ++i)
	{
		mpfr_t* const ptr = c.storage + 0 * c.ncols + i;
		mpfr_mul(tmp1, *ptr, *ptr, MPFR_RNDD);
		mpfr_add(tmp0, tmp0, tmp1, MPFR_RNDD);
	}
	
	mpfr_sqrt(tmp0, tmp0, MPFR_RNDD);	
	mpfr_set_ui(tmp1, 0, MPFR_RNDD);
	
	for (i = 0; i < b.nrows; ++i)
	{
		mpfr_t* const ptr = b.storage + i * b.ncols + 0;
		mpfr_add(tmp1, tmp1, *ptr, MPFR_RNDD);
	}
	
	mpfr_sqrt(tmp1, tmp1, MPFR_RNDD);
	mpfr_ui_div(tmp1, 1, tmp1, MPFR_RNDD);
	mpfr_mul(tmp0, tmp0, tmp1, MPFR_RNDD);
	mpfr_mul_d(tmp0, tmp0, R2, MPFR_RNDD);
	mpfr_log(tmp0, tmp0, MPFR_RNDD);
	
	mpfr_set_d(tmp1, eps, MPFR_RNDD);
	mpfr_abs(tmp1, tmp1, MPFR_RNDD);
	mpfr_log(tmp1, tmp1, MPFR_RNDD);
	
	mpfr_add(tmp0, tmp0, tmp1, MPFR_RNDD);
	
	mpfr_set_ui(tmp1, a.ncols + 1, MPFR_RNDD);
	mpfr_sqr(tmp1, tmp1, MPFR_RNDD);
	mpfr_mul_ui(tmp1, tmp1, 2, MPFR_RNDD);
	mpfr_mul(tmp0, tmp0, tmp1, MPFR_RNDD);
	
	*maxiter = mpfr_get_ui(tmp0, MPFR_RNDD);
	
	if (errno == EDOM || errno == ERANGE)
	{
		return -1;
	}
	
	mpfr_clear(tmp0);
	mpfr_clear(tmp1);
	
	return 0;
}

int em_set_size(size_t* const size,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c)
{
	*size = 1;
	if (em_accum_size(size, a)
			|| em_accum_size(size, b)
			|| em_accum_size(size, c))
	{
		return -1;
	}
	
	*size = log10(*size);
	if (errno == EDOM || errno == ERANGE)
	{
		return -1;
	}
	
	*size = *size + a.nrows * a.ncols;
	
	return 0;
}

int em_accum_size(size_t* const size,
		struct mtx const m)
{	
	int i, j;
	
	for (i = 0; i < m.nrows; ++i)
	{
		for (j = 0; j < m.ncols; ++j)
		{
			long x = mpfr_get_si(*(m.storage + i * m.ncols + j), MPFR_RNDD);
			
			if (errno == ERANGE)
			{
				return -1;
			}
			
			if (x != 0)
			{
				*size = *size * labs(x);
			}
		}
	}
	
	return 0;
}

int em_check_dimensions(struct mtx* const x,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c)
{
	if (a.nrows != b.nrows
			|| a.ncols != c.ncols
			|| b.ncols != 1
			|| c.nrows != 1
			|| x->nrows != a.nrows
			|| x->ncols != 1)
	{
		return -1;
	}
	
	return 0;
}
