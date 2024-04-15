#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
//#include <crtdbg.h>
#include <stdio.h>
#include <string.h>

#include "global.h"


int main(int argc, char *argv[])
{
	double dmin_deg_ratio, dmin_vtx_overlap, dmin_cross_score;
	int nmin_size;

	if(argc!=7)
	{
		printf("Usage:\n");
		printf("\t%s ppi_score_filename min_deg_ratio min_size filter_thres  merger_thres output_file\n", argv[0]);
		return 0;
	}

	dmin_deg_ratio = atof(argv[2]);
	nmin_size = atoi(argv[3]);
	dmin_vtx_overlap = atof(argv[4]);
	dmin_cross_score = atof(argv[5]);

	GenCmplx(argv[1], dmin_deg_ratio, nmin_size, dmin_vtx_overlap, dmin_cross_score, argv[6]);
//
//	_CrtDumpMemoryLeaks();

	return 0;
}


