#include <gmp.h>
#include <stdlib.h>
#include <matrix/matrix.h>


int mtx_init(mpq_t** mtx, size_t const nrows, size_t const ncols)
{
	*mtx = (mpq_t*) malloc(sizeof(mpq_t) * nrows * ncols);
	
	if (NULL == *mtx)
	{
		return -1;
	}
	
	int i, j;
	
	for(i = 0; i < nrows; ++i)
	{
		for(j = 0; j < ncols; ++j)
		{
			mpq_init(*(*mtx + i * ncols + j));
		}
	}
	
	return 0;
}

int mtx_clear(mpq_t** mtx, size_t const nrows, size_t const ncols)
{
	int i, j;
	
	for(i = 0; i < nrows; ++i)
	{
		for(j = 0; j < ncols; ++j)
		{
			mpq_clear(*(*mtx + i * ncols + j));
		}
	}
	
	free(*mtx);
	*mtx = NULL;
	
	return 0;
}

int mtx_mul(mpq_t* res, mpq_t* a, mpq_t* b, size_t const arows, size_t const bcols, size_t const cols)
{
	return 0;
}

size_t mtx_print(FILE* stream, mpq_t* mtx, size_t const nrows, size_t const ncols, int base)
{
	int i, j;
	size_t bytes;
	
	for(i = 0; i < nrows; ++i)
	{
		for(j = 0; j < ncols; ++j)
		{
			bytes += mpq_out_str(stream, base, *(mtx + i * ncols + j));
		}
	}
	
	return bytes;
}
