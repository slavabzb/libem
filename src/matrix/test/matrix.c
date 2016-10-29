#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <matrix/matrix.h>

#include <CUnit/Basic.h>

struct matrix ma;
struct matrix mb;
struct matrix mc;

int suite_init_init(void)
{
	ma.storage = NULL;	
	return 0;
}

int suite_init_clean(void)
{
	free(ma.storage);
	return 0;
}

void test_mtx_init(void)
{
	CU_ASSERT(0 == mtx_init(&ma, 2, 3, 10));
	CU_ASSERT(NULL != ma.storage);
}

int suite_clear_init(void)
{
	ma.nrows = 2;
	ma.ncols = 3;
	ma.storage = (mpfr_t*) malloc(sizeof(mpfr_t) * ma.nrows * ma.ncols);	
	return 0;
}

int suite_clear_clean(void)
{
	return 0;
}

void test_mtx_clear(void)
{
	CU_ASSERT(0 == mtx_clear(&ma));
}

int suite_ops_init(void)
{
	mtx_init(&ma, 2, 3, 10);
	mtx_init(&mb, 3, 4, 10);
	mtx_init(&mc, 2, 4, 10);

	size_t i, j;
	
	for (i = 0; i < ma.nrows; ++i)
	{
		for (j = 0; j < ma.ncols; ++j)
		{
			mpfr_set_ui(*(ma.storage + i * ma.ncols + j), i + j, MPFR_RNDD);
		}
	}
	
	for (i = 0; i < mb.nrows; ++i)
	{
		for (j = 0; j < mb.ncols; ++j)
		{
			mpfr_set_ui(*(mb.storage + i * mb.ncols + j), 2*i + j, MPFR_RNDD);
		}
	}
		
	return 0;
}

int suite_ops_clean(void)
{
	return mtx_clear(&ma) || mtx_clear(&mb) || mtx_clear(&mc);
}

void test_mtx_mul(void)
{
	CU_ASSERT(0 != mtx_mul(&mc, &ma, &ma));
	CU_ASSERT(0 == mtx_mul(&mc, &ma, &mb));
	
	size_t i, j, k;
	
	for (i = 0; i < ma.nrows; ++i)
	{
		for (j = 0; j < mb.ncols; ++j)
		{
			mpfr_t* const prop = mc.storage + i * mc.ncols + j;	
			
			mpfr_t tmp1;
			mpfr_init2(tmp1, mpfr_get_prec(*prop));
			mpfr_set_ui(tmp1, 0, MPFR_RNDD);
			
			for (k = 0; k < ma.ncols; ++k)
			{
				mpfr_t* const pop1 = ma.storage + i * ma.ncols + k;
				mpfr_t* const pop2 = mb.storage + k * mb.ncols + j;
				
				mpfr_t tmp2;
				mpfr_init2(tmp2, mpfr_get_prec(*prop));
				
				mpfr_mul(tmp2, *pop1, *pop2, MPFR_RNDD);
				mpfr_add(tmp1, tmp1, tmp2, MPFR_RNDD);
				
				mpfr_clear(tmp2);
			}
	
			CU_ASSERT(0 == mpfr_cmp(*prop, tmp1));
			
			mpfr_clear(tmp1);
		}
	}
}

void test_mtx_add(void)
{
	CU_ASSERT(0 != mtx_add(&mc, &ma, &mb));
}

int main()
{	
	CU_pSuite suite_init = NULL;
	CU_pSuite suite_clear = NULL;
	CU_pSuite suite_ops = NULL;

	if (CUE_SUCCESS != CU_initialize_registry())
	{
		return CU_get_error();
	}

	suite_init = CU_add_suite("init", suite_init_init, suite_init_clean);
	if (NULL == suite_init)
	{
		CU_cleanup_registry();
		return CU_get_error();
	}

	if ((NULL == CU_add_test(suite_init, "mtx_init", test_mtx_init)))
	{
		CU_cleanup_registry();
		return CU_get_error();
	}

	suite_clear = CU_add_suite("clear", suite_clear_init, suite_clear_clean);
	if (NULL == suite_clear)
	{
		CU_cleanup_registry();
		return CU_get_error();
	}
	
	if ((NULL == CU_add_test(suite_clear, "mtx_clear", test_mtx_clear)))
	{
		CU_cleanup_registry();
		return CU_get_error();
	}
	
	suite_ops = CU_add_suite("ops", suite_ops_init, suite_ops_clean);
	if (NULL == suite_ops)
	{
		CU_cleanup_registry();
		return CU_get_error();
	}

	if ((NULL == CU_add_test(suite_ops, "mtx_mul", test_mtx_mul)) ||
			(NULL == CU_add_test(suite_ops, "mtx_add", test_mtx_add)))
	{
		CU_cleanup_registry();
		return CU_get_error();
	}
	
	CU_basic_set_mode(CU_BRM_VERBOSE);
	CU_basic_run_tests();
	CU_cleanup_registry();

	return CU_get_error();
}
