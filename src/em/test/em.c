#include <matrix/matrix.h>
#include <CUnit/Basic.h>

struct mtx ma;
struct mtx mb;
struct mtx mc;
struct mtx mx;

mpfr_t fx;

mpfr_prec_t const prec = 10;

int suite_init(void)
{
	if (0 != mtx_init(&ma, 2, 3, prec))
	{
		return -1;
	}
	
	if (0 != mtx_init(&mb, 2, 1, prec))
	{
		return -1;
	}
	
	if (0 != mtx_init(&mc, 1, 3, prec))
	{
		return -1;
	}
	
	if (0 != mtx_init(&mx, 2, 1, prec))
	{
		return -1;
	}
	
	mpfr_init2(fx, prec);
	
	return 0;
}

int suite_clean(void)
{
	mpfr_clear(fx);
	
	return mtx_clear(ma)
			|| mtx_clear(mb)
			|| mtx_clear(mc)
			|| mtx_clear(mx);
}

void test_em_optimize(void)
{
	CU_ASSERT(0 == em_oprimize(&fx, &mx, ma, mb, mc));
}

int main()
{
	CU_pSuite suite = NULL;

	if (CUE_SUCCESS != CU_initialize_registry())
	{
		return CU_get_error();
	}

	suite = CU_add_suite("em", suite_init, suite_clean);
	if (NULL == suite)
	{
		CU_cleanup_registry();
		return CU_get_error();
	}

	if ((NULL == CU_add_test(suite, "test_em_optimize", test_em_optimize)))
	{
		CU_cleanup_registry();
		return CU_get_error();
	}

	CU_basic_set_mode(CU_BRM_VERBOSE);
	CU_basic_run_tests();
	CU_cleanup_registry();
	
	return CU_get_error();
}
