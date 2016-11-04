#include <math.h>
#include <stdlib.h>
#include <em/em.h>
#include <errno.h>

int em_set_r2(double* const R2,
		size_t n,
		mpz_t size);

int em_set_default_prec(mpfr_prec_t* const prec,
		struct mtx const a,
		mpz_t size);

int em_set_maxiter(size_t* const maxiter,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c,
		double const R2,
		double const eps);

int em_set_size(mpz_t* const size,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c);

int em_accum_size(mpz_t* const size,
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
	
	mpz_t size;
	mpz_init(size);
	
	if (em_set_size(&size, a, b, c))
		return -1;
	
	mpfr_prec_t prec;
	if (em_set_default_prec(&prec, a, size))
		return -1;
	
	struct mtx y;
	if (mtx_init(&y, a.nrows, 1, prec))
		return -1;
	
	if (mtx_fill_d(y, 0.f, 0.f))
		return -1;
	
	struct mtx H;
	if (mtx_init(&H, a.ncols, a.ncols, prec))
		return -1;
	
	double R2;
	if (em_set_r2(&R2, a.ncols, size))
		return -1;
	
	mpz_clear(size);
	
	if (mtx_fill_d(H, 0, R2))
		return -1;
	
	size_t iter = 0;
	size_t maxiter = 0;
	if (em_set_maxiter(&maxiter, a, b, c, R2, eps))
		return -1;
	
	for (iter = 0; iter < maxiter; ++iter)
	{
		printf("%d ", iter);
	}
	
	return 0;
}

int em_set_r2(double* const R2,
		size_t n,
		mpz_t size)
{
	mpz_t tmp0;
	mpz_init_set(tmp0, size);
	mpz_mul_ui(tmp0, tmp0, 2);
	
	mpfr_t tmp1;
	mpfr_init_set_ui(tmp1, n, MPFR_RNDD);
	mpfr_sqr(tmp1, tmp1, MPFR_RNDD);
	
	mpfr_t tmp2;
	mpfr_init_set_ui(tmp2, 2, MPFR_RNDD);
	mpfr_pow_z(tmp2, tmp2, tmp0, MPFR_RNDD);
	
	mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDD);
	
	*R2 = mpfr_get_d(tmp1, MPFR_RNDD);
	
	mpfr_clear(tmp2);
	mpfr_clear(tmp1);
	
	mpz_clear(tmp0);
	
	return 0;
}

int em_set_default_prec(mpfr_prec_t* const prec,
		struct mtx const a,
		mpz_t size)
{
	mpfr_t tmp0;
	mpfr_t tmp1;
	
	mpfr_init_set_z(tmp0, size, MPFR_RNDD);
	mpfr_mul_ui(tmp0, tmp0, 2, MPFR_RNDD);
	
	mpfr_init_set_ui(tmp1, a.ncols, MPFR_RNDD);
	mpfr_mul_ui(tmp1, tmp1, 2, MPFR_RNDD);
	
	mpfr_add(tmp0, tmp0,  tmp1, MPFR_RNDD);
	
	*prec = mpfr_get_si(tmp0, MPFR_RNDD);
	
	if (errno == ERANGE)
		return -1;
	
	mpfr_set_default_prec(*prec);
	
	mpfr_clear(tmp0);
	mpfr_clear(tmp1);
	
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

int em_set_size(mpz_t* const size,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c)
{
	mpz_set_ui(*size, 1);
	
	if (em_accum_size(size, a)
			|| em_accum_size(size, b)
			|| em_accum_size(size, c))
	{
		return -1;
	}
	
	mpfr_t tmp0;
	mpfr_init_set_z(tmp0, *size, MPFR_RNDD);
	mpfr_log10(tmp0, tmp0, MPFR_RNDD);
	mpfr_add_ui(tmp0, tmp0, a.nrows * a.ncols, MPFR_RNDD);
	
	mpfr_get_z(*size, tmp0, MPFR_RNDD);
	
	if (errno == ERANGE)
	{
		return -1;
	}
	
	mpfr_clear(tmp0);
	
	return 0;
}

int em_accum_size(mpz_t* const size,
		struct mtx const m)
{	
	int i, j;
	
	mpfr_t tmp0;
	mpfr_init(tmp0);
	
	mpz_t tmp1;
	mpz_init(tmp1);
	
	for (i = 0; i < m.nrows; ++i)
	{
		for (j = 0; j < m.ncols; ++j)
		{
			mpfr_t* const ptr = m.storage + i * m.ncols + j;
						
			if (0 != mpfr_cmp_ui(*ptr, 0))
			{
				mpfr_abs(tmp0, *ptr, MPFR_RNDD);
				mpfr_get_z(tmp1, tmp0, MPFR_RNDD);
				mpz_add(*size, *size, tmp1);
			}
		}
	}
	
	mpz_clear(tmp1);
	mpfr_clear(tmp0);
	
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
