#include <math.h>
#include <stdlib.h>
#include <em/em.h>
#include <errno.h>

int em_check_constraints(size_t* const nrow,
		struct mtx const y,
		struct mtx const a,
		struct mtx const b);

int em_set_r2(mpfr_t* const R2,
		size_t n,
		mpz_t size);

int em_set_default_prec(mpfr_prec_t* const prec,
		struct mtx const a,
		mpz_t size);

int em_set_maxiter(mpz_t* const maxiter,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c,
		mpfr_t R2,
		double eps);

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
		mpfr_prec_t* const prec,
		mpz_t* const niters,
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
	
	if (em_set_default_prec(prec, a, size))
		return -1;
		
	if (mtx_fill_d(*x, 0.f, 0.f))
		return -1;
	
	struct mtx H;
	if (mtx_init(&H, a.nrows, a.nrows, *prec))
		return -1;
	
	mpfr_t R2;
	mpfr_init2(R2, *prec);
	if (em_set_r2(&R2, a.nrows, size))
		return -1;
	
	mpz_clear(size);
	
	mpfr_t tmp;
	mpfr_init2(tmp, *prec);
	mpfr_set_ui(tmp, 0, MPFR_RNDD);
	if (mtx_fill(H, tmp, R2))
		return -1;
	
	mpz_t iter;
	if (em_set_maxiter(niters, a, b, c, R2, eps))
		return -1;
	
	mpfr_clear(R2);

	struct mtx g;
	g.nrows = c.nrows;
	g.ncols = c.ncols;
	
	struct mtx gt;
	if (mtx_init(&gt, g.ncols, g.nrows, *prec))
		return -1;
	
	struct mtx Hg;
	if (mtx_init(&Hg, H.nrows, gt.ncols, *prec))
		return -1;
	
	struct mtx Hgtmp;
	if (mtx_init(&Hgtmp, Hg.nrows, Hg.ncols, *prec))
		return -1;
		
	struct mtx gtH;
	if (mtx_init(&gtH, gt.nrows, H.ncols, *prec))
		return -1;
	
	mpfr_t Hgg;
	mpfr_init(Hgg);
	
	for (mpz_init_set_ui(iter, 0); mpz_cmp(iter, *niters) < 0; mpz_add_ui(iter, iter, 1))
	{
		g.storage = c.storage;
		
		size_t idx = 0;
		if (em_check_constraints(&idx, *x, a, b))
		{
			g.storage = a.storage + idx * a.ncols;
		}

		if (mtx_tr(gt, g))
			return -1;
		
		if (mtx_mul(Hg, H, gt))
			return -1;
		
		mpfr_set_ui(Hgg, 0, MPFR_RNDD);
		
		size_t i;
		for (i = 0; i < Hg.nrows; ++i)
		{
			mpfr_t* const Hgptr = Hg.storage + i * Hg.ncols;
			mpfr_t* const gptr = g.storage + 0 * g.ncols + i;
			
			mpfr_mul(tmp, *Hgptr, *gptr, MPFR_RNDD);
			mpfr_add(Hgg, Hgg, tmp, MPFR_RNDD);
		}
		
		mpfr_sqrt(tmp, Hgg, MPFR_RNDD);
		mpfr_mul_ui(tmp, tmp, a.ncols + 1, MPFR_RNDD);
		mpfr_si_div(tmp, -1, tmp, MPFR_RNDD);
		
		if (mtx_mulval(Hgtmp, Hg, tmp))
			return -1;
		
		if (mtx_add(*x, *x, Hgtmp))
			return -1;
		
		
		mtx_fprint(stdout, Hg);
		printf("\nHg size: %d x %d\n", Hg.nrows, Hg.ncols);
		printf("\ng size: %d x %d\n", g.nrows, g.ncols);
		
		gmp_printf("%Zd ", iter);
	}
	
	mpfr_clear(Hgg);
	mpz_clear(iter);
	mpfr_clear(tmp);
	
	if (mtx_clear(Hg) ||
			mtx_clear(gt) ||
			mtx_clear(gtH) ||
			mtx_clear(Hgtmp))
	{
		return -1;
	}
	
	return 0;
}

int em_check_constraints(size_t* const nrow,
		struct mtx const y,
		struct mtx const a,
		struct mtx const b)
{
	size_t i, j;
	
	mpfr_t tmp0;
	mpfr_t sum;
	mpfr_t max;
	
	mpfr_init(tmp0);
	mpfr_init(sum);
	mpfr_init_set_ui(max, 0, MPFR_RNDD);
	
	for (i = 0; i < a.nrows; ++i)
	{
		mpfr_t* const yptr = y.storage + i * y.ncols;
		mpfr_t* const bptr = b.storage + i * b.ncols;
				
		mpfr_set_ui(sum, 0, MPFR_RNDD);
		
		for (j = 0; j < a.ncols; ++j)
		{
			mpfr_t* const aptr = a.storage + i * a.ncols + j;
			
			mpfr_mul(tmp0, *aptr, *yptr, MPFR_RNDD);
			mpfr_add(sum, sum, tmp0, MPFR_RNDD);
		}
		
		if (mpfr_cmp(sum, *bptr) >= 0)
		{
			mpfr_sub(tmp0, *bptr, sum, MPFR_RNDD);
			
			if (mpfr_cmp(tmp0, max) > 0)
			{
				*nrow = i;
				mpfr_set(max, tmp0, MPFR_RNDD);
			}
		}
	}
	
	int rop = mpfr_cmp_ui(max, 0);
	
	mpfr_clear(max);
	mpfr_clear(sum);
	mpfr_clear(tmp0);
	
	return rop;
}

int em_set_r2(mpfr_t* const R2,
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
	
	mpfr_set(*R2, tmp1, MPFR_RNDD);
	
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

int em_set_maxiter(mpz_t* const maxiter,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c,
		mpfr_t R2,
		double eps)
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
	mpfr_mul(tmp0, tmp0, R2, MPFR_RNDD);
	mpfr_log(tmp0, tmp0, MPFR_RNDD);
	
	mpfr_set_d(tmp1, eps, MPFR_RNDD);
	mpfr_abs(tmp1, tmp1, MPFR_RNDD);
	mpfr_log(tmp1, tmp1, MPFR_RNDD);
	
	mpfr_add(tmp0, tmp0, tmp1, MPFR_RNDD);
	
	mpfr_set_ui(tmp1, a.ncols + 1, MPFR_RNDD);
	mpfr_sqr(tmp1, tmp1, MPFR_RNDD);
	mpfr_mul_ui(tmp1, tmp1, 2, MPFR_RNDD);
	mpfr_mul(tmp0, tmp0, tmp1, MPFR_RNDD);
	
	mpfr_get_z(*maxiter, tmp0, MPFR_RNDD);
	
	if (errno == ERANGE)
		return -1;
	
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
		return -1;
	
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
