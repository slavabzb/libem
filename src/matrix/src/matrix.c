#include <omp.h>
#include <stdlib.h>
#include <matrix/matrix.h>


int mtx_init(struct mtx* const m, size_t rows, size_t columns, mpfr_prec_t prec)
{
	m->nrows = rows;
	m->ncols = columns;
	m->storage = (mpfr_t*) malloc(sizeof(mpfr_t) * m->nrows * m->ncols);
	
	if (NULL == m->storage)
	{
		return -1;
	}
	
	size_t i, j;
	
#pragma omp parallel for private(i, j) schedule(static)
	for(i = 0; i < m->nrows; ++i)
	{
		for(j = 0; j < m->ncols; ++j)
		{
			mpfr_init2(*(m->storage + i * m->ncols + j), prec);
		}
	}

	return 0;
}

int mtx_clear(struct mtx const m)
{
	size_t i, j;
	
#pragma omp parallel for private(i, j) schedule(static)
	for(i = 0; i < m.nrows; ++i)
	{
		for(j = 0; j < m.ncols; ++j)
		{
			mpfr_clear(*(m.storage + i * m.ncols + j));
		}
	}

	free(m.storage);
	
	return 0;
}

int mtx_fprint(FILE* stream, struct mtx const m)
{
	int total = 0;
	int chars = 0;
	
	size_t i, j;
	
	for(i = 0; i < m.nrows; ++i)
	{
		for(j = 0; j < m.ncols; ++j)
		{
			chars = mpfr_fprintf(stream, "%.2Rf ", *(m.storage + i * m.ncols + j));
			
			if (chars < 0)
			{
				return 0;
			}
			
			total += chars;
		}
		
		chars = fprintf(stream, "\n");
		
		if (chars < 0)
		{
			return 0;
		}

		total += chars;
	}
	
	return total;
}

int mtx_mul(struct mtx rop, struct mtx const op1, struct mtx const op2)
{
	if (rop.nrows != op1.nrows || rop.ncols != op2.ncols)
	{
		return -1;
	}
	
	if (op1.ncols != op2.nrows)
	{
		return -1;
	}
	
	mpfr_prec_t const prec = mpfr_get_prec(*rop.storage);
	
	size_t i, j, k;

#pragma omp parallel for shared(rop) private(i, j, k) schedule(static)
	for (i = 0; i < op1.nrows; ++i)
	{
		for (j = 0; j < op2.ncols; ++j)
		{
			mpfr_t* const prop = rop.storage + i * rop.ncols + j;	
			mpfr_set_ui(*prop, 0, MPFR_RNDD);

			for (k = 0; k < op1.ncols; ++k)
			{
				mpfr_t* const pop1 = op1.storage + i * op1.ncols + k;
				mpfr_t* const pop2 = op2.storage + k * op2.ncols + j;

				mpfr_t tmp;
				mpfr_init2(tmp, prec);

				mpfr_mul(tmp, *pop1, *pop2, MPFR_RNDD);
				mpfr_add(*prop, *prop, tmp, MPFR_RNDD);

				mpfr_clear(tmp);
			}
		}
	}

	return 0;
}

int mtx_add(struct mtx rop, struct mtx const op1, struct mtx const op2)
{
	if (op1.nrows != op2.nrows || op1.ncols != op2.ncols)
	{
		return -1;
	}
	
	if (rop.nrows != op1.nrows || rop.ncols != op1.ncols)
	{
		return -1;
	}
	
	size_t i, j;

#pragma omp parallel for shared(rop) private(i, j) schedule(static)
	for (i = 0; i < rop.nrows; ++i)
	{
		for (j = 0; j < rop.ncols; ++j)
		{
			mpfr_add(*(rop.storage + i * rop.ncols + j), *(op1.storage + i * op1.ncols + j), *(op2.storage + i * op2.ncols + j), MPFR_RNDD);
		}
	}

	return 0;
}

int mtx_tr(struct mtx rop, struct mtx const op)
{
	if (rop.nrows != op.ncols || rop.ncols != op.nrows)
	{
		return -1;
	}
	
	size_t i, j;

#pragma omp parallel for shared(rop) private(i, j) schedule(static)	
	for (i = 0; i < rop.nrows; ++i)
	{
		for (j = 0; j < rop.ncols; ++j)
		{
			mpfr_set(*(rop.storage + i * rop.ncols + j), *(op.storage + j * op.ncols + i), MPFR_RNDD);
		}
	}

	return 0;
}
