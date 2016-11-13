#include <matrix/matrix.h>
#include <CUnit/Basic.h>

struct mtx ma;
struct mtx mb;
struct mtx mc;
struct mtx mx;

mpfr_t fx;
mpz_t niters;

mpfr_prec_t prec = 10;
double const eps = 1e-3;

int suite_init(void)
{
	if (mtx_init(&ma, 3, 2, prec))
		return -1;
	
	mpfr_set_ui(*(ma.storage + 0 * ma.ncols + 0), 1, MPFR_RNDN);
	mpfr_set_ui(*(ma.storage + 0 * ma.ncols + 1), 0, MPFR_RNDN);
	mpfr_set_ui(*(ma.storage + 1 * ma.ncols + 0), 0, MPFR_RNDN);
	mpfr_set_ui(*(ma.storage + 1 * ma.ncols + 1), 2, MPFR_RNDN);
	mpfr_set_ui(*(ma.storage + 2 * ma.ncols + 0), 3, MPFR_RNDN);
	mpfr_set_ui(*(ma.storage + 2 * ma.ncols + 1), 2, MPFR_RNDN);
	
	if (mtx_init(&mb, 3, 1, prec))
		return -1;
	
	mpfr_set_ui(*(mb.storage + 0 * mb.ncols + 0), 180, MPFR_RNDN);
	mpfr_set_ui(*(mb.storage + 1 * mb.ncols + 0), 150, MPFR_RNDN);
	mpfr_set_ui(*(mb.storage + 2 * mb.ncols + 0), 300, MPFR_RNDN);
	
	if (mtx_init(&mc, 1, 2, prec))
		return -1;
	
	mpfr_set_si(*(mc.storage + 0 * mc.ncols + 0), -3, MPFR_RNDN);
	mpfr_set_si(*(mc.storage + 0 * mc.ncols + 1), -5, MPFR_RNDN);
	
	if (mtx_init(&mx, 2, 1, prec))
		return -1;
	
	mpfr_init2(fx, prec);
	mpz_init(niters);
	
	printf("\nma\n");
	mtx_fprint(stdout, ma);
	printf("\nmb\n");
	mtx_fprint(stdout, mb);
	printf("\nmc\n");
	mtx_fprint(stdout, mc);

	printf("\nprec: %ld, eps: %f\n", prec, eps);
	
	return 0;
}

int suite_clean(void)
{
	mpz_clear(niters);
	mpfr_clear(fx);
	
	return mtx_clear(ma)
			|| mtx_clear(mb)
			|| mtx_clear(mc)
			|| mtx_clear(mx);
}

void test_em_optimize(void)
{
	CU_ASSERT(0 == em_optimize(&fx, &mx, &prec, &niters, ma, mb, mc, eps));
	
	mpfr_printf("\nfx: %Rf\n", fx);
	printf("\nmx\n");
	mtx_fprint(stdout, mx);
	gmp_printf("\nprec: %ld\nniters: %Zd\n", prec, niters);
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
