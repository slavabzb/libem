#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <em/em.h>
#include <matrix/matrix.h>

void print_help(char* name);

int main (int argc, char *argv[]) 
{
	if (argc != 2)
	{
		print_help(argv[0]);
		return -1;
	}
		
	FILE* fp = fopen(argv[1], "r");
	if (NULL == fp)
	{
		printf("Can't open file '%s'\n", argv[1]);
		return -1;
	}
	
	char buf[LINE_MAX];
	
	fgets(buf, LINE_MAX, fp);
	fgets(buf, LINE_MAX, fp);
	
	int const M = atoi(buf);
	
	fgets(buf, LINE_MAX, fp);
	fgets(buf, LINE_MAX, fp);
	fgets(buf, LINE_MAX, fp);
		
	int const N = atoi(buf);
	
	fgets(buf, LINE_MAX, fp);
	fgets(buf, LINE_MAX, fp);
	//fgets(buf, LINE_MAX, fp);
	
	//printf("%s", buf);
	
	struct mtx mc;
	mtx_init(&mc, 1, 50, 10);
	
	if (mtx_fscan(fp, mc, "\t"))
		return -1;
	
	mtx_fprint(stdout, mc);
	
	while (fgets(buf, LINE_MAX, fp))
	{
		printf("%s", buf);
		break;
	}
	
	if (fclose(fp))
		return -1;
	
	if (mtx_clear(mc))
		return -1;
	
	return 0;
}

void print_help(char* name)
{
	printf("Usage:\n");
	printf("\t%s FILE\n\n", name);
	printf("FILE - file containing an optimization problem\n");
}
