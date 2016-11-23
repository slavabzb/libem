#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <em/em.h>
#include <matrix/matrix.h>

size_t const PREC = 53;

void print_help(char const* const name);

int read_mtx(char const* const file,
		struct mtx ma,
		struct mtx mc,
		struct mtx mb0,
		struct mtx mb1);

int main (int argc, char *argv[]) 
{
	if (argc != 2)
	{
		print_help(argv[0]);
		return -1;
	}
	
	struct mtx ma;
	mtx_init(&ma, 125, 50, PREC);
	
	struct mtx mb0;
	mtx_init(&mb0, ma.nrows, 1, PREC);
	
	struct mtx mb1;
	mtx_init(&mb1, ma.nrows, 1, PREC);
	
	struct mtx mc;
	mtx_init(&mc, 1, ma.ncols, PREC);
	
	struct mtx mx;
	mtx_init(&mx, mc.ncols, mc.nrows, PREC);
	
	if (read_mtx(argv[1], ma, mc, mb0, mb1))
		return -1;
	
	//mtx_fprint(stdout, ma);
	//mtx_fprint(stdout, mc);
	//mtx_fprint(stdout, mb0);
	//mtx_fprint(stdout, mb1);
	
	if (mtx_clear(ma) || mtx_clear(mb1) || mtx_clear(mb0) || mtx_clear(mc) || mtx_clear(mx))
		return -1;
	
	return 0;
}

void print_help(char const* const name)
{
	printf("Usage:\n");
	printf("\t%s FILE\n\n", name);
	printf("FILE - file containing an optimization problem\n");
}

int read_mtx(char const* const file,
		struct mtx ma,
		struct mtx mc,
		struct mtx mb0,
		struct mtx mb1)
{
	FILE* fp = fopen(file, "r");
	
	if (NULL == fp)
	{
		printf("Can't open file '%s'\n", file);
		return -1;
	}
	
	struct mtx mb;
	mtx_init(&mb, 1, 125, PREC);
	
	char buf[LINE_MAX];
	
	fgets(buf, LINE_MAX, fp);
	if (mtx_fscan(fp, mc, "\t"))
		return -1;
	
	fgets(buf, LINE_MAX, fp);
	if (mtx_fscan(fp, mb, "\t"))
		return -1;
	
	if (mtx_tr(mb0, mb))
		return -1;
	
	fgets(buf, LINE_MAX, fp);
	if (mtx_fscan(fp, mb, "\t"))
		return -1;
	
	if (mtx_tr(mb1, mb))
		return -1;
	
	fgets(buf, LINE_MAX, fp);
	if (mtx_fscan(fp, ma, "\t"))
		return -1;
	
	if (fclose(fp))
		return -1;
	
	if (mtx_clear(mb))
		return -1;
	
	return 0;
}
