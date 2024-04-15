#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/timeb.h>

#include "global.h"

//#define EXEC_LOCATION  ""
#define EXEC_LOCATION  ".\\exec\\"
//#define EXEC_LOCATION  "./exec/"

void DelFile(char* szfilename)
{
	char szcmd[500];

	sprintf(szcmd, "del %s\n", szfilename);
	//sprintf(szcmd, "rm %s\n", szfilename);
	printf(szcmd);
	system(szcmd);
}



void StrUpr(char* szname)
{
	int i;

	i = 0;
	while(szname[i]!=0)
	{
		if(szname[i]>='a' && szname[i]<='z')
			szname[i] = szname[i] + 'A'-'a';
		i++;
	}
}


int GetRowNum(char* szfilename)
{
	FILE *fp;
	char ch;
	int num_of_rows;

	num_of_rows = 0;

	fp = fopen(szfilename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szfilename);
		return 0;
	}

	ch = fgetc(fp);
	while(!feof(fp))
	{
		while(!feof(fp) && ch!='\n')
			ch = fgetc(fp);
		num_of_rows++;

		ch = fgetc(fp);
	}
	fclose(fp);

	return num_of_rows;
}

int LoadIds(char* szpid_filename, vector<string> *pvec_ids)
{
	FILE *fp;
	char ch, szpname[100];
	int num_of_ids, nlen;

	pvec_ids->clear();
	num_of_ids = 0;

	fp = fopen(szpid_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szpid_filename);
		return 0;
	}

	ch = fgetc(fp);
	while(!feof(fp))
	{
		nlen = 0;
		while(!feof(fp) && ch!='\n')
		{
			if(ch!='\r')
				szpname[nlen++] = ch;
			ch = fgetc(fp);
		}
		szpname[nlen] = 0;
		
		if(nlen>0)
			StrUpr(szpname);
		pvec_ids->push_back(szpname);
		num_of_ids++;

		ch = fgetc(fp);
	}
	fclose(fp);

	return num_of_ids;
}

int LoadIdMap(char* szpid_filename, map<string, int> *ppid_map)
{
	FILE *fp;
	char ch, szpname[100];
	map<string, int>::iterator  map_it;
	int num_of_ids, nlen;

	ppid_map->clear();
	num_of_ids = 0;

	fp = fopen(szpid_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szpid_filename);
		return 0;
	}

	ch = fgetc(fp);
	while(!feof(fp))
	{
		nlen = 0;
		while(!feof(fp) && ch!='\n')
		{
			if(ch!='\r')
				szpname[nlen++] = ch;
			ch = fgetc(fp);
		}
		szpname[nlen] = 0;
		
		if(nlen>0)
		{
			StrUpr(szpname);
			map_it = ppid_map->find(szpname);
			if(map_it==ppid_map->end())
				(*ppid_map)[szpname] = num_of_ids;
			else
				printf("Error: duplicated id name %s \n", szpname);
		}
		num_of_ids++;

		ch = fgetc(fp);
	}
	fclose(fp);

	printf("#ids: %d\n", num_of_ids);

	return num_of_ids;
}


int GenIdMap(char* szppi_filename, char* szpidmap_filename)
{
	FILE *fp, *fpmap;
	char szpname[100];
	map<string, int> pid_map;
	map<string, int>::iterator map_it;
	int num_of_ids, num_of_ppis;

	fp = fopen(szppi_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szppi_filename);
		return 0;
	}
	fpmap = fopen(szpidmap_filename, "wt");
	if(fpmap==NULL)
	{
		printf("Error: cannot open file %s for write\n", szpidmap_filename);
		return 0;
	}

	num_of_ppis = 0;
	num_of_ids = 0;
	fscanf(fp, "%s", szpname);
	while(!feof(fp))
	{
		StrUpr(szpname);
		map_it = pid_map.find(szpname);
		if(map_it==pid_map.end())
		{
			pid_map[szpname] = num_of_ids;
			num_of_ids++;
			fprintf(fpmap, "%s\n", szpname);
		}

		fscanf(fp, "%s", szpname);
		StrUpr(szpname);
		map_it = pid_map.find(szpname);
		if(map_it==pid_map.end())
		{
			pid_map[szpname] = num_of_ids;
			num_of_ids++;
			fprintf(fpmap, "%s\n", szpname);
		}
		num_of_ppis++;

		fscanf(fp, "%s", szpname);
	}
	fclose(fp);
	fclose(fpmap);

	printf("#PPIs: %d, #distinct proteins: %d\n", num_of_ppis, num_of_ids);

	return num_of_ids;
}

void LoadPPIScoreMatrix(char* szppi_score_filename, char* szprtn_ids_filename, double *pdppiscore_matrix)
{
	FILE *fp;
	char szpname1[100], szpname2[100];
	double dscore, davg_score;
	int npid1, npid2, num_of_ppis;
	map<string, int> pid_map;
	map<string, int>::iterator map_it;

	LoadIdMap(szprtn_ids_filename, &pid_map);

	fp = fopen(szppi_score_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szppi_score_filename);
		return;
	}

	num_of_ppis = 0;
	davg_score = 0;

	fscanf(fp, "%s", szpname1);
	while(!feof(fp))
	{
		fscanf(fp, "%s", szpname2);
		fscanf(fp, "%lf", &dscore);

		StrUpr(szpname1);
		map_it = pid_map.find(szpname1);
		if(map_it!=pid_map.end())
			npid1 = map_it->second;
		else
		{
			npid1 = -1;
			printf("Error: cannot find protein %s in the map\n", szpname1);
		}

		StrUpr(szpname2);
		map_it = pid_map.find(szpname2);
		if(map_it!=pid_map.end())
			npid2 = map_it->second;
		else
		{
			npid2 = -1;
			printf("Error: cannot find protein %s in the map\n", szpname2);
		}

		if(npid1>=0 && npid2>=0)
		{
			if(npid1>npid2)
			{
				if(pdppiscore_matrix[npid1*(npid1+1)/2+npid2]==0)
				{
					pdppiscore_matrix[npid1*(npid1+1)/2+npid2] = dscore;
					num_of_ppis++;
					davg_score += dscore;
				}
				else if(pdppiscore_matrix[npid1*(npid1+1)/2+npid2]!=dscore)
					printf("Error: inconsistent score %f %f\n", pdppiscore_matrix[npid1*(npid1+1)/2+npid2], dscore);
			}
			else if(npid1<npid2)
			{
				if(pdppiscore_matrix[npid2*(npid2+1)/2+npid1]==0)
				{
					pdppiscore_matrix[npid2*(npid2+1)/2+npid1] = dscore;
					num_of_ppis++;
					davg_score += dscore;
				}
				else if(pdppiscore_matrix[npid2*(npid2+1)/2+npid1]!=dscore)
					printf("Error: inconsistent score %f %f\n", pdppiscore_matrix[npid2*(npid2+1)/2+npid1], dscore);
			}
		}
		fscanf(fp, "%s", szpname1);
	}
	fclose(fp);

	printf("distinct #PPIs in file %s: %d\n", szppi_score_filename, num_of_ppis);
	davg_score /= num_of_ppis;
	printf("Average score: %.3f\n", davg_score);
}

void GetPPIs(char* szppi_score_filename, char *szoutput_filename)
{
	FILE *fp, *fpout;
	char szpname1[200], szpname2[200];
	double dscore;
	int num_of_ppis;

	fp = fopen(szppi_score_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szppi_score_filename);
		return;
	}
	fpout = fopen(szoutput_filename, "a+");
	if(fpout==NULL)
	{
		printf("Error: cannot open file %s for write\n", szoutput_filename);
		return;
	}

	num_of_ppis = 0;

	fscanf(fp, "%s", szpname1);
	while(!feof(fp))
	{
		fscanf(fp, "%s", szpname2);
		fscanf(fp, "%lf", &dscore);
		num_of_ppis++;

		fprintf(fpout, "%s\t%s\n", szpname1, szpname2);

		fscanf(fp, "%s", szpname1);
	}
	fclose(fp);
	fclose(fpout);

	printf("Total PPIs: %d\n", num_of_ppis);
}


void ConvertPPI(char* szrawppi_filename, char* szprtn_id_filename, char* szppi_pair_filename, char* szppi_matrix_filename)
{
	char szcmd[500];

	sprintf(szcmd, "%sConvertPPI %s %s %s %s\n", EXEC_LOCATION, szrawppi_filename, szprtn_id_filename, szppi_pair_filename, szppi_matrix_filename);
	printf(szcmd);
	system(szcmd);
}

void quasiCliques(char* szgraph_filename, double dmin_deg_ratio, int nmin_size, int nmax_size, char* szoutput_filename)
{
	char szcmd[500];

	sprintf(szcmd, "%squasiCliques %s %.3f %d %d %s\n", EXEC_LOCATION, szgraph_filename, dmin_deg_ratio, nmin_size, nmax_size, szoutput_filename);
	printf(szcmd);
	system(szcmd);
}

