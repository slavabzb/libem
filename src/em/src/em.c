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
	{
		return -1;
	}
	
	size_t size;
	em_set_size(&size, a, b, c);
	
	size_t iter = 0;
	struct mtx y;
	struct mtx B;
	
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
