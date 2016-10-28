#include <gmp.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <matrix/matrix.h>

#include <CUnit/Basic.h>

mpq_t* ma = NULL;
mpq_t* mb = NULL;
mpq_t* mc = NULL;

size_t const arows = 2;
size_t const acols = 3;

size_t const brows = 3;
size_t const bcols = 4;

size_t const crows = 2;
size_t const ccols = 4;

int suite_init_init(void)
{
	ma = NULL;
	return 0;
}

int suite_init_clean(void)
{
	free(*ma);
	return 0;
}

void test_mtx_init(void)
{
	CU_ASSERT(0 == mtx_init(&ma, arows, acols));
	CU_ASSERT(NULL != ma);
}

int suite_clear_init(void)
{
	ma = (mpq_t*) malloc(sizeof(mpq_t) * arows * acols);
	mtx_print(stdout, ma, arows, acols, 10);
	return 0;
}

int suite_clear_clean(void)
{
	free(*ma);
	return 0;
}

void test_mtx_clear(void)
{
	CU_ASSERT(0 == mtx_clear(&ma, arows, acols));
	CU_ASSERT(NULL == ma);
}

int suite_ops_init(void)
{
	int ra = mtx_init(&ma, arows, acols);
	int rb = mtx_init(&mb, brows, bcols);
	int rc = mtx_init(&mc, crows, ccols);
	return ra || rb || rc;
}

int suite_ops_clean(void)
{
	int ra = mtx_clear(&ma, arows, acols);
	int rb = mtx_clear(&mb, brows, bcols);
	int rc = mtx_clear(&mc, crows, ccols);
	return ra || rb || rc;
}

void test_mtx_mul(void)
{
	CU_ASSERT(0 == mtx_mul(mc, ma, mb, arows, bcols, acols));
}

void test_mtx_add(void)
{
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
