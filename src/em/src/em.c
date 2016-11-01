#include <math.h>
#include <stdlib.h>
#include <em/em.h>
#include <errno.h>

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

int em_oprimize(mpfr_t* const fx,
		struct mtx* const x,
		struct mtx const a,
		struct mtx const b,
		struct mtx const c)
{
	if (em_check_dimensions(x, a, b, c))
		return -1;
	
	size_t size;
	em_set_size(&size, a, b, c);
	
	mpfr_prec_t const prec = 2 * log2(a.ncols) + 2 * size;
	if (errno == EDOM || errno == ERANGE)
		return -1;
	
	struct mtx y;
	if (mtx_init(&y, a.nrows, 1, prec))
		return -1;
	
	if (mtx_fill(y, 0.f, 0.f))
		return -1;
	
	struct mtx H;
	if (mtx_init(&H, a.ncols, a.ncols, prec))
		return -1;
	if (mtx_fill(H, 0, pow(a.ncols, 2) * pow(2, 2 * size)))
		return -1;
	
	size_t iter;
	
	for (iter = 0; iter < 1; ++iter)
	{
		
	}
	
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
