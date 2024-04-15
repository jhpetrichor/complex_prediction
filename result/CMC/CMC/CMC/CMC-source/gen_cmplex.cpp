#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <cstring>

#include "global.h"

int *gpdiff_set1;
int *gpdiff_set2;

struct QCLQ
{
	int nsize;
	int *pprtn_set;
	double dscore;
	int nmax_overlap_size;
};

int comp_int(const void *e1, const void *e2)
{
	int n1, n2;

	n1 = *(int*)e1;
	n2 = *(int*)e2;

	if(n1<n2)
		return -1;
	else if(n1>n2)
		return 1;
	else 
		return 0;
}

int LoadQClqs(char* szqclq_filename, int num_of_prtns, QCLQ* &pqclqs)
{
	FILE *fp;
	int num_of_qclqs, nclq_no, nsize, i, nprtn_in_qclqs;
	bool *pbprtn_flags;

	num_of_qclqs = GetRowNum(szqclq_filename);
	if(num_of_qclqs==0)
	{
		pqclqs = NULL;
		return 0;
	}



	fp = fopen(szqclq_filename, "rt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for read\n", szqclq_filename);
		return 0;
	}

	pbprtn_flags = new bool[num_of_prtns];
	memset(pbprtn_flags, 0, sizeof(bool)*num_of_prtns);

	pqclqs = new QCLQ[num_of_qclqs];

	nclq_no = 0;
	fscanf(fp, "%d", &nsize);
	while(!feof(fp))
	{
		pqclqs[nclq_no].nsize = nsize;
		pqclqs[nclq_no].pprtn_set = new int[nsize];

		for(i=0;i<nsize;i++)
		{
			fscanf(fp, "%d", &pqclqs[nclq_no].pprtn_set[i]);
			pbprtn_flags[pqclqs[nclq_no].pprtn_set[i]] = true;
		}
		qsort(pqclqs[nclq_no].pprtn_set, pqclqs[nclq_no].nsize, sizeof(int), comp_int);

		nclq_no++;
		fscanf(fp, "%d", &nsize);
	}
	fclose(fp);
	if(nclq_no!=num_of_qclqs)
		printf("Error: inconsistent number of quasi-cliques\n");

	nprtn_in_qclqs = 0;
	for(i=0;i<num_of_prtns;i++)
	{
		if(pbprtn_flags[i])
			nprtn_in_qclqs++;
	}
	delete []pbprtn_flags;

	printf("#quasi-cliques: %d\n", num_of_qclqs);
	printf("#proteins in quasi-cliques: %d\n", nprtn_in_qclqs);

	return num_of_qclqs;
}

double CalcQclqScore(QCLQ *pqclq, double *pdppiscore_matrix)
{
	int i, j, npid1, npid2, nppis;
	double dweight_sum, dscore;

	nppis = 0;
	dweight_sum = 0;
	for(i=0;i<pqclq->nsize;i++)
	{
		npid1 = pqclq->pprtn_set[i];
		for(j=0;j<i;j++)
		{
			npid2 = pqclq->pprtn_set[j];
			if(pdppiscore_matrix[npid1*(npid1+1)/2+npid2]>0)
				nppis++;
			dweight_sum += pdppiscore_matrix[npid1*(npid1+1)/2+npid2];
		}
	}

	//dscore = (double)nppis*2/(pqclq->nsize-1);
	//dscore = (double)nppis*2/(pqclq->nsize*(pqclq->nsize-1));
	dscore = dweight_sum*2/(pqclq->nsize*(pqclq->nsize-1));
	//dscore = dweight_sum*2/(pqclq->nsize-1);
	//dscore = log((double)pqclq->nsize)*dweight_sum*2/(pqclq->nsize*(pqclq->nsize-1));

	return dscore;
}

int comp_qclq_score(const void *e1, const void *e2)
{
	QCLQ *p1, *p2;

	p1 = (QCLQ*)e1;
	p2 = (QCLQ*)e2;

	if(p1->dscore > p2->dscore)
		return -1;
	else if(p1->dscore < p2->dscore)
		return 1;
	else if(p1->nsize > p2->nsize)
		return -1;
	else if(p1->nsize < p2->nsize)
		return 1;
	else
		return 0;
}


int GetDiff(int *pset1, int nlen1, int *pset2, int nlen2, int *pdiff1_2, int* pdiff2_1)
{
	int i, j, k, t;

	i = 0;
	j = 0;
	k = 0;
	t = 0;
	while(i<nlen1 && j<nlen2)
	{
		if(pset1[i]<pset2[j])
			pdiff1_2[k++] = pset1[i++];
		else if(pset1[i]>pset2[j])
			pdiff2_1[t++] = pset2[j++];
		else
		{
			i++;
			j++;
		}
	}
	while(i<nlen1)
		pdiff1_2[k++] = pset1[i++];
	while(j<nlen2)
		pdiff2_1[t++] = pset2[j++];

	if(nlen1-k!=nlen2-t)
		printf("Error: inconsistent number of common items\n");

	return nlen1-k;
}

double CalcOverlapScore(QCLQ *pqclq1, QCLQ* pqclq2, double *pdppiscore_matrix, int &ncmmn_prtns)
{
	int *pdiff_set1, *pdiff_set2, ndiff_size1, ndiff_size2, i, j, num_of_interacts1, num_of_interacts2; 
	double dppi_score, dppi_score_sum1, dppi_score_sum2, doverlap_score;

	pdiff_set1 = gpdiff_set1;
	pdiff_set2 = gpdiff_set2;

	ncmmn_prtns = GetDiff(pqclq1->pprtn_set, pqclq1->nsize, pqclq2->pprtn_set, pqclq2->nsize, pdiff_set1, pdiff_set2);
	ndiff_size1 = pqclq1->nsize-ncmmn_prtns;
	ndiff_size2 = pqclq2->nsize-ncmmn_prtns;
	if(ndiff_size1==0)
	{
		printf("Error: non-maximal clique\n");
		return 0;
	}
	if(ndiff_size2==0)
		return 1;
	if(ncmmn_prtns==0)
		return 0;

/*
	dppi_score_sum1 = 0;
	for(i=0;i<ndiff_size1;i++)
	{
		for(j=0;j<ndiff_size2;j++)
		{
			if(pdiff_set1[i]>pdiff_set2[j] && pdppiscore_matrix[pdiff_set1[i]*(pdiff_set1[i]+1)/2+pdiff_set2[j]])
				dppi_score_sum1 += pdppiscore_matrix[pdiff_set1[i]*(pdiff_set1[i]+1)/2+pdiff_set2[j]];
			else if(pdiff_set1[i]<pdiff_set2[j] && pdppiscore_matrix[pdiff_set2[j]*(pdiff_set2[j]+1)/2+pdiff_set1[i]])
				dppi_score_sum1 += pdppiscore_matrix[pdiff_set2[j]*(pdiff_set2[j]+1)/2+pdiff_set1[i]];
		}
	}
	doverlap_score = dppi_score_sum1/(ndiff_size1*ndiff_size2);
*/


	dppi_score_sum1 = 0;
	num_of_interacts1 = 0;
	for(i=0;i<ndiff_size1;i++)
	{
		for(j=0;j<pqclq2->nsize;j++)
		{
			if(pdiff_set1[i]>pqclq2->pprtn_set[j] && pdppiscore_matrix[pdiff_set1[i]*(pdiff_set1[i]+1)/2+pqclq2->pprtn_set[j]]>0)
				dppi_score = pdppiscore_matrix[pdiff_set1[i]*(pdiff_set1[i]+1)/2+pqclq2->pprtn_set[j]];
			else if(pdiff_set1[i]<pqclq2->pprtn_set[j] && pdppiscore_matrix[pqclq2->pprtn_set[j]*(pqclq2->pprtn_set[j]+1)/2+pdiff_set1[i]]>0)
				dppi_score = pdppiscore_matrix[pqclq2->pprtn_set[j]*(pqclq2->pprtn_set[j]+1)/2+pdiff_set1[i]];
			else
				dppi_score = 0;
			if(dppi_score>0)
				num_of_interacts1++;
			dppi_score_sum1 += dppi_score;
		}
	}
	dppi_score_sum1 = dppi_score_sum1/(ndiff_size1*pqclq2->nsize);

	dppi_score_sum2 = 0;
	num_of_interacts2 = 0;
	for(i=0;i<ndiff_size2;i++)
	{
		for(j=0;j<pqclq1->nsize;j++)
		{
			if(pdiff_set2[i]>pqclq1->pprtn_set[j] && pdppiscore_matrix[pdiff_set2[i]*(pdiff_set2[i]+1)/2+pqclq1->pprtn_set[j]])
				dppi_score = pdppiscore_matrix[pdiff_set2[i]*(pdiff_set2[i]+1)/2+pqclq1->pprtn_set[j]];
			else if(pdiff_set2[i]<pqclq1->pprtn_set[j] && pdppiscore_matrix[pqclq1->pprtn_set[j]*(pqclq1->pprtn_set[j]+1)/2+pdiff_set2[i]])
				dppi_score = pdppiscore_matrix[pqclq1->pprtn_set[j]*(pqclq1->pprtn_set[j]+1)/2+pdiff_set2[i]];
			else
				dppi_score = 0;
			if(dppi_score>0)
				num_of_interacts2++;
			dppi_score_sum2 += dppi_score;
		}
	}
	dppi_score_sum2 = dppi_score_sum2/(ndiff_size2*pqclq1->nsize);

	doverlap_score = sqrt(dppi_score_sum1*dppi_score_sum2);


	return doverlap_score;
}

int Union(int *pset1, int nlen1, int* pset2, int nlen2, int *presult)
{
	int i, j, k;

	i = 0;
	j = 0;
	k = 0;
	while(i<nlen1 && j<nlen2)
	{
		if(pset1[i]<pset2[j])
			presult[k++] = pset1[i++];
		else if(pset1[i]>pset2[j])
			presult[k++] = pset2[j++];
		else
		{
			presult[k++] = pset1[i++];
			j++;
		}
	}
	while(i<nlen1)
		presult[k++] = pset1[i++];
	while(j<nlen2)
		presult[k++] = pset2[j++];

	return k;
}

void GetPrtnCoverStat(int *pprtn_cover_times, int num_of_prtns)
{
	int num_of_prtn_covered, nmax_cover_times, nmin_cover_times, i, *pcovertimes_distr;
	double davg_cover_times;

	num_of_prtn_covered = 0;
	davg_cover_times = 0;
	nmax_cover_times = 0;
	nmin_cover_times = -1;
	for(i=0;i<num_of_prtns;i++)
	{
		if(pprtn_cover_times[i]>0)
		{
			num_of_prtn_covered++;
			davg_cover_times += pprtn_cover_times[i];
			if(nmax_cover_times < pprtn_cover_times[i])
				nmax_cover_times = pprtn_cover_times[i];
			if(nmin_cover_times==-1 || nmin_cover_times > pprtn_cover_times[i])
				nmin_cover_times = pprtn_cover_times[i];
		}
	}
	davg_cover_times /= num_of_prtn_covered;

	pcovertimes_distr = new int[nmax_cover_times+1];
	memset(pcovertimes_distr, 0, sizeof(int)*(nmax_cover_times+1));
	for(i=0;i<num_of_prtns;i++)
		pcovertimes_distr[pprtn_cover_times[i]]++;	

	printf("#number of proteins: %d, #proteins being covered: %d\n", num_of_prtns, num_of_prtn_covered);
	printf("Cover times: avg=%.2f, min=%d, max=%d\n", davg_cover_times, nmin_cover_times, nmax_cover_times);
	for(i=1;i<=nmax_cover_times;i++)
	{
		if(pcovertimes_distr[i]>0)
			printf("%d   %d\n", i, pcovertimes_distr[i]);
	}
	delete []pcovertimes_distr;
}


int MergeCmplx(QCLQ *pqclqs, int num_of_qclqs, double *pdppiscore_matrix, int num_of_prtns, double dfilter_overlap_thres, double dmerge_thres)
{
	int i, j, *punion_set, nunion_len, ncommon_prtns, ncovered_prtns, nmax_cluster_size;
	double dscore;
	int *pprtn_cover_times, ntop_k, num_of_merges, num_of_incmerges, num_of_removals, num_of_remain_clqs, niteration_num;
	bool bmerged;

	printf("\n............ Removing and merging highly overlaped cliques ............\n");
	printf("#total cliques before removing and merging: %d\n", num_of_qclqs);

	ntop_k = 0;
	nmax_cluster_size = 100;

	gpdiff_set1 = new int[num_of_prtns];
	gpdiff_set2 = new int[num_of_prtns];
	pprtn_cover_times = new int[num_of_prtns];

	punion_set = gpdiff_set1;

	num_of_merges = 0;
	num_of_incmerges = 0;
	num_of_removals = 0;
	niteration_num = 1;


	bmerged = true;
	while(bmerged)
	{
		memset(pprtn_cover_times, 0, sizeof(int)*num_of_prtns);
		for(i=0;i<num_of_qclqs;i++)
			pqclqs[i].nmax_overlap_size = 0;

		bmerged = false;
		for(i=0;i<num_of_qclqs;i++)
		{
			if(i>ntop_k && pqclqs[i].nsize>0)
			{
				ncovered_prtns = 0;
				for(j=0;j<pqclqs[i].nsize;j++)
				{
					if(pprtn_cover_times[pqclqs[i].pprtn_set[j]]>=1)
						ncovered_prtns++;
				}
				if((double)pqclqs[i].nmax_overlap_size/pqclqs[i].nsize>=dfilter_overlap_thres && pqclqs[i].nsize-ncovered_prtns<=1)
				{
					delete []pqclqs[i].pprtn_set;
					pqclqs[i].pprtn_set = NULL;
					pqclqs[i].nsize = 0;
					num_of_removals++;
				}
			}
			if(pqclqs[i].nsize>0)
			{
				for(j=i+1;j<num_of_qclqs;j++)
				{
					if(pqclqs[j].nsize>0)
					{
						dscore = CalcOverlapScore(&pqclqs[i], &pqclqs[j], pdppiscore_matrix, ncommon_prtns);
						//if(dmerge_thres<1 && pqclqs[i].nsize<nmax_cluster_size && (double)ncommon_prtns/pqclqs[j].nsize>=dfilter_overlap_thres && dscore>=dmerge_thres)
						if(dmerge_thres<1 && (double)ncommon_prtns/pqclqs[j].nsize>=dfilter_overlap_thres && dscore>=dmerge_thres)
						{
							nunion_len = Union(pqclqs[i].pprtn_set, pqclqs[i].nsize, pqclqs[j].pprtn_set, pqclqs[j].nsize, punion_set);
							num_of_merges++;
							if(nunion_len>pqclqs[i].nsize)
							{
								num_of_incmerges++;
								bmerged = true;
							}
							delete []pqclqs[i].pprtn_set;
							pqclqs[i].nsize = nunion_len;
							pqclqs[i].pprtn_set = new int[nunion_len];
							memcpy(pqclqs[i].pprtn_set, punion_set, sizeof(int)*nunion_len);
							if(j>ntop_k)
							{
								delete []pqclqs[j].pprtn_set;
								pqclqs[j].pprtn_set = NULL;
								pqclqs[j].nsize = 0;
								pqclqs[j].dscore = 0;
								num_of_removals++;
							}
						}
						if(pqclqs[j].nmax_overlap_size<ncommon_prtns)
							pqclqs[j].nmax_overlap_size = ncommon_prtns;
					}
				}
				for(j=0;j<pqclqs[i].nsize;j++)
					pprtn_cover_times[pqclqs[i].pprtn_set[j]]++;
				pqclqs[i].dscore = CalcQclqScore(&pqclqs[i], pdppiscore_matrix);
			}
		}
		num_of_remain_clqs = 0;
		for(i=0;i<num_of_qclqs;i++)
		{
			if(pqclqs[i].nsize>0)
				pqclqs[num_of_remain_clqs++] = pqclqs[i];
		}
		printf("#cliques after %d-th round: %d\n", niteration_num, num_of_remain_clqs);
		num_of_qclqs = num_of_remain_clqs;
		niteration_num++;
	}

	GetPrtnCoverStat(pprtn_cover_times, num_of_prtns);

	delete []gpdiff_set1;
	delete []gpdiff_set2;
	delete []pprtn_cover_times;

	printf("#cliques removed: %d, #merge: %d, #merge resulting in size increase: %d\n", num_of_removals, num_of_merges, num_of_incmerges);


	return num_of_qclqs;
}


void OutputCmplex(QCLQ *pqclqs, int num_of_qclqs, char* szprtn_ids_filename, char* szoutput_filename)
{
	FILE *fp;
	int i, j, ncomplex_no, nmax_cmplx_size;
	vector<string> vec_pids;


	fp = fopen(szoutput_filename, "wt");
	if(fp==NULL)
	{
		printf("Error: cannot open file %s for write\n", szoutput_filename);
		return;
	}

	nmax_cmplx_size = 0;
	ncomplex_no = 0;

	if(num_of_qclqs>0)
	{
		LoadIds(szprtn_ids_filename, &vec_pids);

		for(i=0;i<num_of_qclqs;i++)
		{
			if(pqclqs[i].nsize>0 && pqclqs[i].pprtn_set!=NULL)
			{
				fprintf(fp, "C%d(%d_%.3f): ", ncomplex_no, pqclqs[i].nsize, pqclqs[i].dscore);
				for(j=0;j<pqclqs[i].nsize;j++)
					fprintf(fp, "%s ", vec_pids[pqclqs[i].pprtn_set[j]].c_str());
				fprintf(fp, "\n");
				delete []pqclqs[i].pprtn_set;
				if(nmax_cmplx_size<pqclqs[i].nsize)
					nmax_cmplx_size = pqclqs[i].nsize;
				ncomplex_no++;
			}
		}
		delete []pqclqs;
	}
	fclose(fp);


	printf("\n");
	printf("#clusters: %d\n", ncomplex_no);
	printf("Maximal cluster size: %d\n", nmax_cmplx_size);
}

int GetOverlap(int *pset1, int nlen1, int *pset2, int nlen2)
{
	int i, j, noverlap;

	noverlap = 0;
	i = 0;
	j = 0;
	while(i<nlen1 && j<nlen2)
	{
		if(pset1[i]<pset2[j])
			i++;
		else if(pset1[i]>pset2[j])
			j++;
		else
		{
			i++;
			j++;
			noverlap++;
		}
	}
	return noverlap;
}



struct CLQ_PAIR
{
	int nclq_id1;
	int nclq_id2;
	double doverlap_score;
};

int comp_clq_pair(const void *e1, const void *e2)
{
	CLQ_PAIR *p1, *p2;

	p1 = (CLQ_PAIR*)e1;
	p2 = (CLQ_PAIR*)e2;

	if(p1->doverlap_score > p2->doverlap_score)
		return -1;
	else if(p1->doverlap_score < p2->doverlap_score)
		return 1;
	else if(p1->nclq_id2 > p2->nclq_id2)
		return -1;
	else if(p1->nclq_id2 < p2->nclq_id2)
		return 1;
	else
		return 0;
}


int RmvOverlapCmplx(QCLQ *pqclqs, int num_of_qclqs, double *pdppiscore_matrix, int num_of_prtns, double dfilter_overlap_thres, int nmax_clq_num)
{
	int i, j, noverlap, num_of_clq_pairs, nclq_id, num_of_clq_rmved, ntop_k, *pcovertimes_distr;
	double doverlap_score, davg_cover_times; //dcross_score;
	int *pprtn_cover_times, nmin_prtn_cover_times, num_of_prtn_covered, nmin_cover_times, nmax_cover_times;
	CLQ_PAIR *poverlap_clq_pairs;
	int ncvr_clq_id; //nunion_len, *punion_set;
	bool *pbused_flag;

	if(num_of_qclqs<=nmax_clq_num)
	{
		printf("No cliques are removed\n");
		return num_of_qclqs;
	}

	printf("\n............ Removing highly overlaped cliques ............\n");

	ntop_k = 50;

	nmin_prtn_cover_times = 1;

	//gpdiff_set1 = new int[num_of_prtns];
	//gpdiff_set2 = new int[num_of_prtns];
	//punion_set = gpdiff_set1;

	pprtn_cover_times = new int[num_of_prtns];
	memset(pprtn_cover_times, 0, sizeof(int)*num_of_prtns);

	num_of_clq_pairs = 0;
	for(i=0;i<num_of_qclqs;i++)
	{
		for(j=0;j<pqclqs[i].nsize;j++)
			pprtn_cover_times[pqclqs[i].pprtn_set[j]]++;
		for(j=i+1;j<num_of_qclqs;j++)
		{
			if(1 || pqclqs[j].nsize <= pqclqs[i].nsize)
			{
				noverlap = GetOverlap(pqclqs[i].pprtn_set, pqclqs[i].nsize, pqclqs[j].pprtn_set, pqclqs[j].nsize);
				doverlap_score = (double)noverlap/pqclqs[j].nsize;
				if(doverlap_score>=dfilter_overlap_thres)
					num_of_clq_pairs++;
			}
		}
	}

	poverlap_clq_pairs = new CLQ_PAIR[num_of_clq_pairs];
	num_of_clq_pairs = 0;
	for(i=0;i<num_of_qclqs;i++)
	{
		for(j=i+1;j<num_of_qclqs;j++)
		{
			if(1 || pqclqs[j].nsize <= pqclqs[i].nsize)
			{
				noverlap = GetOverlap(pqclqs[i].pprtn_set, pqclqs[i].nsize, pqclqs[j].pprtn_set, pqclqs[j].nsize);
				doverlap_score = (double)noverlap/pqclqs[j].nsize;
				if(doverlap_score>=dfilter_overlap_thres)
				{
					poverlap_clq_pairs[num_of_clq_pairs].nclq_id1 = i;
					poverlap_clq_pairs[num_of_clq_pairs].nclq_id2 = j;
					poverlap_clq_pairs[num_of_clq_pairs].doverlap_score = doverlap_score;
					num_of_clq_pairs++;
				}
			}
		}
	}
	qsort(poverlap_clq_pairs, num_of_clq_pairs, sizeof(CLQ_PAIR), comp_clq_pair);

	pbused_flag = new bool[num_of_qclqs];
	memset(pbused_flag, 0, sizeof(bool)*num_of_qclqs);

	num_of_clq_rmved = 0;
	for(i=0;i<num_of_clq_pairs;i++)
	{
		nclq_id = poverlap_clq_pairs[i].nclq_id2;
		ncvr_clq_id = poverlap_clq_pairs[i].nclq_id1;
		if(pqclqs[nclq_id].nsize>0 && pqclqs[ncvr_clq_id].nsize>0)// && !pbused_flag[nclq_id])
		{
			/*
			dcross_score = CalcOverlapScore(&pqclqs[ncvr_clq_id], &pqclqs[nclq_id], pdppiscore_matrix, noverlap, dfilter_overlap_thres);
			//merging if possible
			if(dcross_score>=dmerge_thres)
			{
				for(j=0;j<pqclqs[ncvr_clq_id].nsize;j++)
					pprtn_cover_times[pqclqs[ncvr_clq_id].pprtn_set[j]]--;
				nunion_len = Union(pqclqs[ncvr_clq_id].pprtn_set, pqclqs[ncvr_clq_id].nsize, pqclqs[nclq_id].pprtn_set, pqclqs[nclq_id].nsize, punion_set);
				delete []pqclqs[ncvr_clq_id].pprtn_set;
				pqclqs[ncvr_clq_id].nsize = nunion_len;
				pqclqs[ncvr_clq_id].pprtn_set = new int[nunion_len];
				memcpy(pqclqs[ncvr_clq_id].pprtn_set, punion_set, sizeof(int)*nunion_len);
				for(j=0;j<pqclqs[ncvr_clq_id].nsize;j++)
					pprtn_cover_times[pqclqs[ncvr_clq_id].pprtn_set[j]]++;
				bmerged = true;
			}
			*/

			//removing
			if(nclq_id>=ntop_k || poverlap_clq_pairs[i].doverlap_score>=0.9)
			{
				for(j=0;j<pqclqs[nclq_id].nsize;j++)
				{
					if(pprtn_cover_times[pqclqs[nclq_id].pprtn_set[j]]<=nmin_prtn_cover_times)
						break;
				}
				if(j>=pqclqs[nclq_id].nsize)
				{
					pbused_flag[ncvr_clq_id] = true;
					for(j=0;j<pqclqs[nclq_id].nsize;j++)
						pprtn_cover_times[pqclqs[nclq_id].pprtn_set[j]]--;
					delete []pqclqs[nclq_id].pprtn_set;
					pqclqs[nclq_id].nsize = 0;
					num_of_clq_rmved++;
					if(num_of_qclqs-num_of_clq_rmved<=nmax_clq_num)
						break;
				}
			}
		}
	}
	delete []poverlap_clq_pairs;
	//delete []gpdiff_set1;
	//delete []gpdiff_set2;
	delete []pbused_flag;


	num_of_prtn_covered = 0;
	davg_cover_times = 0;
	nmax_cover_times = 0;
	nmin_cover_times = -1;
	for(i=0;i<num_of_prtns;i++)
	{
		if(pprtn_cover_times[i]>0)
		{
			num_of_prtn_covered++;
			davg_cover_times += pprtn_cover_times[i];
			if(nmax_cover_times < pprtn_cover_times[i])
				nmax_cover_times = pprtn_cover_times[i];
			if(nmin_cover_times==-1 || nmin_cover_times > pprtn_cover_times[i])
				nmin_cover_times = pprtn_cover_times[i];
		}
	}
	davg_cover_times /= num_of_prtn_covered;
	pcovertimes_distr = new int[nmax_cover_times+1];
	memset(pcovertimes_distr, 0, sizeof(int)*(nmax_cover_times+1));
	for(i=0;i<num_of_prtns;i++)
		pcovertimes_distr[pprtn_cover_times[i]]++;	

	printf("Total cliques: %d, #cliques removed: %d\n", num_of_qclqs, num_of_clq_rmved);
	printf("#number of proteins: %d, #proteins being covered: %d\n", num_of_prtns, num_of_prtn_covered);
	printf("Cover times: avg=%.2f, min=%d, max=%d\n", davg_cover_times, nmin_cover_times, nmax_cover_times);
	for(i=1;i<=nmax_cover_times;i++)
	{
		if(pcovertimes_distr[i]>0)
			printf("%d   %d\n", i, pcovertimes_distr[i]);
	}
	delete []pprtn_cover_times;
	delete []pcovertimes_distr;


	return num_of_qclqs-num_of_clq_rmved;
}



void GenCmplx(char* szppi_score_filename, double dmin_deg_ratio, int nmin_size, double dfilter_overlap_thres, double dmin_cross_score, char* szoutput_filename)
{
	char szppi_filename[200], szprtn_ids_filename[200];
	QCLQ *pqclqs;
	int num_of_qclqs, num_of_prtns, i;
	vector<string> vecpids;
	double *pdppiscore_matrix, dclq_mine_time, dtotal_run_time;
	struct timeb start, clq_mine_start, end;

	DelFile("ppi.temp.txt");
	DelFile("prtn.temp.ids");
	DelFile("ppi.pairs.txt");
	DelFile("ppi.matrix.txt");
	DelFile("qclq.txt");

	ftime(&start);

	sprintf(szppi_filename, "ppi.temp.txt");
	GetPPIs(szppi_score_filename, szppi_filename);

	sprintf(szprtn_ids_filename, "prtn.temp.ids");
	num_of_prtns = GenIdMap(szppi_filename, szprtn_ids_filename);
	printf("#proteins: %d\n", num_of_prtns);

	ConvertPPI(szppi_filename, szprtn_ids_filename, "ppi.pairs.txt", "ppi.matrix.txt");

	printf("\n--------------------- Mining cliques --------------------- \n");
	ftime(&clq_mine_start);
	quasiCliques("ppi.matrix.txt", dmin_deg_ratio, nmin_size, num_of_prtns, "qclq.txt");
	ftime(&end);
	dclq_mine_time = end.time-clq_mine_start.time+(double)(end.millitm-clq_mine_start.millitm)/1000;
	printf("Maximal clique mining time: %.3f\n", dclq_mine_time);
	printf("--------------------- Mining finished --------------------- \n\n");

	num_of_qclqs = LoadQClqs("qclq.txt", num_of_prtns, pqclqs);

	if(num_of_qclqs>0)
	{
		pdppiscore_matrix = new double[num_of_prtns*(num_of_prtns+1)/2];
		memset(pdppiscore_matrix, 0, sizeof(double)*num_of_prtns*(num_of_prtns+1)/2);
		LoadPPIScoreMatrix(szppi_score_filename, szprtn_ids_filename, pdppiscore_matrix);

		for(i=0;i<num_of_qclqs;i++)
			pqclqs[i].dscore = CalcQclqScore(&pqclqs[i], pdppiscore_matrix);
		qsort(pqclqs, num_of_qclqs, sizeof(QCLQ), comp_qclq_score);

		//if(dmin_cross_score>1)
		//	num_of_qclqs = RmvOverlapCmplx(pqclqs, num_of_qclqs, pdppiscore_matrix, num_of_prtns, dfilter_overlap_thres, (int)dmin_cross_score);
		if(dfilter_overlap_thres<1 || dmin_cross_score<1)
			num_of_qclqs = MergeCmplx(pqclqs, num_of_qclqs, pdppiscore_matrix, num_of_prtns, dfilter_overlap_thres, dmin_cross_score);
			
		delete []pdppiscore_matrix;
	}

	OutputCmplex(pqclqs, num_of_qclqs, szprtn_ids_filename, szoutput_filename);

	ftime(&end);
	dtotal_run_time = end.time-start.time+(double)(end.millitm-start.millitm)/1000;
	printf("Total runing time: %.3f\n", dtotal_run_time);

}



