#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <matrix/matrix.h>

#include <CUnit/Basic.h>

struct mtx ma;
struct mtx mb;
struct mtx mc;

struct mtx msquare;

mpfr_t val;
mpfr_t diagval;

double val_d;
double diagval_d;

mpfr_t mulval;

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
	CU_ASSERT(0 == mtx_clear(ma));
}

int suite_ops_init(void)
{
	if (mtx_init(&ma, 2, 3, 10))
	{
		return -1;
	}
	
	if (mtx_init(&mb, 3, 4, 10))
	{
		return -1;
	}
	
	if (mtx_init(&mc, 2, 4, 10))
	{
		return -1;
	}
	
	if (mtx_init(&msquare, 2, 2, 10))
	{
		return -1;
	}
	
	size_t i, j;
	
	for (i = 0; i < ma.nrows; ++i)
	{
		for (j = 0; j < ma.ncols; ++j)
		{
			mpfr_set_ui(*(ma.storage + i * ma.ncols + j), i + j, MPFR_RNDN);
		}
	}
	
	for (i = 0; i < mb.nrows; ++i)
	{
		for (j = 0; j < mb.ncols; ++j)
		{
			mpfr_set_ui(*(mb.storage + i * mb.ncols + j), 2*i + j, MPFR_RNDN);
		}
	}
	
	val_d = 1.f;
	diagval_d = 2.f;
	
	mpfr_init_set_d(val, val_d, MPFR_RNDN);
	mpfr_init_set_d(diagval, diagval_d, MPFR_RNDN);
	
	mpfr_init_set_d(mulval, 3.f, MPFR_RNDN);
	
	return 0;
}

int suite_ops_clean(void)
{
	mpfr_clear(val);
	mpfr_clear(diagval);
	mpfr_clear(mulval);
	
	return mtx_clear(ma) ||
			mtx_clear(mb) ||
			mtx_clear(mc) ||
			mtx_clear(msquare);
}

void test_mtx_fill(void)
{
	CU_ASSERT(0 == mtx_fill(mb, val, diagval));
	
	int i, j;
	
	for (i = 0; i < mb.nrows; ++i)
	{
		for (j = 0; j < mb.ncols; ++j)
		{
			mpfr_t* const ptr = mb.storage + i * mb.ncols + j;
			CU_ASSERT(0 == mpfr_cmp(*ptr, val));
		}
	}
	
	CU_ASSERT(0 == mtx_fill(msquare, val, diagval));
	
	for (i = 0; i < msquare.nrows; ++i)
	{
		for (j = 0; j < msquare.ncols; ++j)
		{
			mpfr_t* const ptr = msquare.storage + i * msquare.ncols + j;
			
			if (i == j)
			{
				CU_ASSERT(0 == mpfr_cmp(*ptr, diagval));
			}
			else
			{
				CU_ASSERT(0 == mpfr_cmp(*ptr, val));
			}
		}
	}
}

void test_mtx_fill_d(void)
{
	CU_ASSERT(0 == mtx_fill_d(mb, val_d, diagval_d));
	
	int i, j;
	
	for (i = 0; i < mb.nrows; ++i)
	{
		for (j = 0; j < mb.ncols; ++j)
		{
			mpfr_t* const ptr = mb.storage + i * mb.ncols + j;
			CU_ASSERT(0 == mpfr_cmp_d(*ptr, val_d));
		}
	}
	
	CU_ASSERT(0 == mtx_fill_d(msquare, val_d, diagval_d));
	
	for (i = 0; i < msquare.nrows; ++i)
	{
		for (j = 0; j < msquare.ncols; ++j)
		{
			mpfr_t* const ptr = msquare.storage + i * msquare.ncols + j;
			
			if (i == j)
			{
				CU_ASSERT(0 == mpfr_cmp_d(*ptr, diagval_d));
			}
			else
			{
				CU_ASSERT(0 == mpfr_cmp_d(*ptr, val_d));
			}
		}
	}
}

void test_mtx_mul(void)
{
	CU_ASSERT(0 != mtx_mul(ma, ma, mb));
	CU_ASSERT(0 != mtx_mul(mc, ma, ma));
	CU_ASSERT(0 == mtx_mul(mc, ma, mb));
	
	size_t i, j, k;
	
	for (i = 0; i < ma.nrows; ++i)
	{
		for (j = 0; j < mb.ncols; ++j)
		{
			mpfr_t* const prop = mc.storage + i * mc.ncols + j;	
			
			mpfr_t tmp1;
			mpfr_init2(tmp1, mpfr_get_prec(*prop));
			mpfr_set_ui(tmp1, 0, MPFR_RNDN);
			
			for (k = 0; k < ma.ncols; ++k)
			{
				mpfr_t* const pop1 = ma.storage + i * ma.ncols + k;
				mpfr_t* const pop2 = mb.storage + k * mb.ncols + j;
				
				mpfr_t tmp2;
				mpfr_init2(tmp2, mpfr_get_prec(*prop));
				
				mpfr_mul(tmp2, *pop1, *pop2, MPFR_RNDN);
				mpfr_add(tmp1, tmp1, tmp2, MPFR_RNDN);
				
				mpfr_clear(tmp2);
			}
	
			CU_ASSERT(0 == mpfr_cmp(*prop, tmp1));
			
			mpfr_clear(tmp1);
		}
	}
}

void test_mtx_mulval(void)
{
	CU_ASSERT(0 != mtx_mulval(mc, mb, mulval));
	CU_ASSERT(0	== mtx_mulval(ma, ma, mulval));
	
	int i, j;
	
	for (i = 0; i < ma.nrows; ++i)
	{
		for (j = 0; j < ma.ncols; ++j)
		{
			mpfr_t tmp;
			mpfr_init2(tmp, 10);
			mpfr_div(tmp, *(ma.storage + i * ma.ncols + j), mulval, MPFR_RNDN);
			CU_ASSERT(0 == mpfr_cmp_ui(tmp, i + j));
			mpfr_clear(tmp);
		}
	}
}

void test_mtx_add(void)
{
	CU_ASSERT(0 != mtx_add(mc, mb, ma));
	CU_ASSERT(0 != mtx_add(mc, mb, mb));
	CU_ASSERT(0 == mtx_add(mb, mb, mb));
	
	size_t i, j;
	
	for (i = 0; i < mb.nrows; ++i)
	{
		for (j = 0; j < mb.ncols; ++j)
		{
			mpfr_t tmp;
			mpfr_init_set_d(tmp, 2 * val_d, MPFR_RNDN);
			CU_ASSERT(0 == mpfr_cmp(*(mb.storage + i * mb.ncols + j), tmp));
			mpfr_clear(tmp);
		}
	}
}

void test_mtx_tr(void)
{
	CU_ASSERT(0 != mtx_tr(mc, ma));
	
	struct mtx tr;
	mtx_init(&tr, ma.ncols, ma.nrows, 10);
	
	CU_ASSERT(0 == mtx_tr(tr, ma));
	
	size_t i, j;
	
	for (i = 0; i < tr.nrows; ++i)
	{
		for (j = 0; j < tr.ncols; ++j)
		{
			CU_ASSERT(0 == mpfr_cmp(*(tr.storage + i * tr.ncols + j), *(ma.storage + j * ma.ncols + i)));
		}
	}
	
	mtx_clear(tr);
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

	if ((NULL == CU_add_test(suite_ops, "mtx_fill", test_mtx_fill)) ||
			(NULL == CU_add_test(suite_ops, "mtx_fill_d", test_mtx_fill_d)) ||
			(NULL == CU_add_test(suite_ops, "mtx_mul", test_mtx_mul)) ||
			(NULL == CU_add_test(suite_ops, "mtx_mulval", test_mtx_mulval)) ||
			(NULL == CU_add_test(suite_ops, "mtx_add", test_mtx_add)) ||
			(NULL == CU_add_test(suite_ops, "mtx_tr", test_mtx_tr)))
	{
		CU_cleanup_registry();
		return CU_get_error();
	}
	
	CU_basic_set_mode(CU_BRM_VERBOSE);
	CU_basic_run_tests();
	CU_cleanup_registry();

	return CU_get_error();
}
