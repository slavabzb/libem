#include <stdlib.h>
#include <matrix/matrix.h>


int mtx_init(struct matrix* const m, size_t rows, size_t columns, mpfr_prec_t prec)
{
	m->nrows = rows;
	m->ncols = columns;
	m->storage = (mpfr_t*) malloc(sizeof(mpfr_t) * m->nrows * m->ncols);
	
	if (NULL == m->storage)
	{
		return -1;
	}
	
	size_t i, j;
	
	for(i = 0; i < m->nrows; ++i)
	{
		for(j = 0; j < m->ncols; ++j)
		{
			mpfr_init2(*(m->storage + i * m->ncols + j), prec);
		}
	}
	
	return 0;
}

int mtx_clear(struct matrix const* const m)
{
	size_t i, j;
	
	for(i = 0; i < m->nrows; ++i)
	{
		for(j = 0; j < m->ncols; ++j)
		{
			mpfr_clear(*(m->storage + i * m->ncols + j));
		}
	}
	
	free(m->storage);
	
	return 0;
}

int mtx_print(FILE* stream, struct matrix const* const m)
{
	int total = 0;
	int chars = 0;
	
	size_t i, j;
	
	for(i = 0; i < m->nrows; ++i)
	{
		for(j = 0; j < m->ncols; ++j)
		{
			chars = mpfr_fprintf(stream, "%.2Rf ", *(m->storage + i * m->ncols + j));
			
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

int mtx_mul(struct matrix* const rop, struct matrix const* const op1, struct matrix const* const op2)
{
	if (op1->ncols != op2->nrows)
	{
		return -1;
	}
	
	mpfr_prec_t const prec = mpfr_get_prec(*rop->storage);
	
	size_t i, j, k;
	
	for (i = 0; i < op1->nrows; ++i)
	{
		for (j = 0; j < op2->ncols; ++j)
		{
			mpfr_t* const prop = rop->storage + i * rop->ncols + j;	
			mpfr_set_ui(*prop, 0, MPFR_RNDD);
			
			for (k = 0; k < op1->ncols; ++k)
			{
				mpfr_t* const pop1 = op1->storage + i * op1->ncols + k;
				mpfr_t* const pop2 = op2->storage + k * op2->ncols + j;
				
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

int mtx_add(struct matrix* const rop, struct matrix const* const op1, struct matrix const* const op2)
{
	if (op1->nrows != op2->nrows || op1->ncols != op2->ncols)
	{
		return -1;
	}
	
	size_t i, j;
}
