#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <iterator>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>
#include <ctime>

using namespace std;

const long double Threshold = 1e-6;

struct InteractionRaw {
	string protein1, protein2;
};

struct Interaction {
	int protein_id1, protein_id2;
};

struct PPINetwork {
	int total_protein;
	string *protein_name;

	int **adjacency_matrix;
	vector<int> *neighbor_list;

	int *vertex_degree;
	int **num_common_neighbor;
};

struct Complex {
	vector<int> core;
	vector<int> attachment;
	vector<int> whole;
	long double score;
};

struct ValueInfo {
	long double value;
	int index;
};

long double *log_factor;
/*****************************************************************************************************************************************/
/* ------------------------------------------------------------------------------------------------------------------------------------ */
void PrintVector(const vector<int> &v)
{
	ostream_iterator<int> output(cout, " ");
	copy(v.begin(), v.end(), output);
	cout << endl;
}

int CountCommon(const vector<int> &v1, const vector<int> &v2)
{
	vector<int> v(v1.size() + v2.size());
	merge(v1.begin(), v1.end(), v2.begin(), v2.end(), v.begin());
	sort(v.begin(), v.end());
	int total_common = 0;
	for (int i = 1; i < (int)v.size(); ++i)
		if (v[i] == v[i-1])
			++ total_common;
	return total_common;
} 

long double LogCombination(int n, int r)
{
	if (r+r > n) 
		r = n-r;
	if (n < 0 || r < 0)
		return 0;
	return log_factor[n] - log_factor[n-r] - log_factor[r];
}

void Flagging(bool *flag, const vector<int> &set)
{
	for (int i = 0; i < (int)set.size(); ++i)
		flag[set[i]] = true;
}

/* ------------------------------------------------------------------------------------------------------------------------------------ */ 
void ReadPPI(const char *filename, vector<InteractionRaw> &raw_interactions)
{
	FILE *fpi = fopen(filename, "r");

	char protein_name1[128];
	char protein_name2[128];

	raw_interactions.resize(0);
	while (fscanf(fpi, "%s%s", protein_name1, protein_name2) != EOF)
	{
		InteractionRaw e;
		e.protein1 = protein_name1;
		e.protein2 = protein_name2;
		raw_interactions.push_back(e);
	}

	fclose(fpi);
}

/* ------------------------------------------------------------------------------------------------------------------------------------ */
void UpdateMapping(const string name, int &id, map<string, int> &name_to_id, map<int, string> &id_to_name)
{
	if (name_to_id[name] == 0)
	{
		++ id;
		name_to_id[name] = id;
		id_to_name[id-1] = name;
	}
}

int SetMapping(const vector<InteractionRaw> &raw_interactions, map<string, int> &name_to_id, map<int, string> &id_to_name)
{
	name_to_id.clear();
	id_to_name.clear();
	int current_id = 0;
	for (int i = 0; i < (int)raw_interactions.size(); ++i)
	{
		UpdateMapping(raw_interactions[i].protein1, current_id, name_to_id, id_to_name);
		UpdateMapping(raw_interactions[i].protein2, current_id, name_to_id, id_to_name);
	}
	for (int i = 0; i < current_id; ++i)
		-- name_to_id[id_to_name[i]];

	return current_id;
}

void RecordName(map<int, string> &id_to_name, int max_id, string *name)
{
	for (int i = 0; i < max_id; ++i)
		name[i] = id_to_name[i];
}

void MapInteractions(const vector<InteractionRaw> &raw_interactions, map<string, int> &name_to_id, vector<Interaction> &interactions)
{
	interactions.resize(0);
	for (int i = 0; i < (int)raw_interactions.size(); ++i)
	{
		Interaction ita;
		ita.protein_id1 = name_to_id[raw_interactions[i].protein1];
		ita.protein_id2 = name_to_id[raw_interactions[i].protein2];
		interactions.push_back(ita);
	}
}

void PrepareAdjacencyMatrix(const vector<Interaction> &interactions, int matrix_size, int **matrix)
{
	for (int i = 0; i < matrix_size; ++i)
		fill_n(matrix[i], matrix_size, 0);

	for (int i = 0; i < (int)interactions.size(); ++i)
	{
		matrix[interactions[i].protein_id1][interactions[i].protein_id2] = 1;
		matrix[interactions[i].protein_id2][interactions[i].protein_id1] = 1;
	}
}

void PrepareNeighborhoodList(const vector<Interaction> &interactions, int total_list, vector<int> *neighbor_list)
{
	for (int i = 0; i < total_list; ++i)
		neighbor_list[i].resize(0);

	for (int i = 0; i < (int)interactions.size(); ++i)
	{
		neighbor_list[interactions[i].protein_id1].push_back(interactions[i].protein_id2);
		neighbor_list[interactions[i].protein_id2].push_back(interactions[i].protein_id1);
	}
}

void SetVertexDegree(const PPINetwork &ppi_network, int *vertex_degree)
{
	for (int i = 0; i < ppi_network.total_protein; ++i)
		vertex_degree[i] = (int)ppi_network.neighbor_list[i].size();
}

void CountCommonNeighbors(const PPINetwork &ppi_network, int **num_common_neighbors)
{
	for (int i = 0; i < ppi_network.total_protein; ++i)
		for (int j = i; j < ppi_network.total_protein; ++j)
			num_common_neighbors[i][j] = num_common_neighbors[j][i] = CountCommon(ppi_network.neighbor_list[i], ppi_network.neighbor_list[j]);
}

void CreateNetwork(const vector<InteractionRaw> &raw_interactions, PPINetwork &ppi_network)
{
	map<string, int> name_to_id;
	map<int, string> id_to_name;
	ppi_network.total_protein = SetMapping(raw_interactions, name_to_id, id_to_name);

	ppi_network.protein_name = new string[ppi_network.total_protein];
	RecordName(id_to_name, ppi_network.total_protein, ppi_network.protein_name);

	vector<Interaction> interactions;
	MapInteractions(raw_interactions, name_to_id, interactions);

	ppi_network.adjacency_matrix = new int *[ppi_network.total_protein];
	for (int i = 0; i < ppi_network.total_protein; ++i)
		ppi_network.adjacency_matrix[i] = new int[ppi_network.total_protein];
	PrepareAdjacencyMatrix(interactions, ppi_network.total_protein, ppi_network.adjacency_matrix);

	ppi_network.neighbor_list = new vector<int>[ppi_network.total_protein];
	PrepareNeighborhoodList(interactions, ppi_network.total_protein, ppi_network.neighbor_list);

	ppi_network.vertex_degree = new int[ppi_network.total_protein];
	SetVertexDegree(ppi_network, ppi_network.vertex_degree);

	ppi_network.num_common_neighbor = new int *[ppi_network.total_protein];
	for (int i = 0; i < ppi_network.total_protein; ++i)
		ppi_network.num_common_neighbor[i] = new int[ppi_network.total_protein];
	CountCommonNeighbors(ppi_network, ppi_network.num_common_neighbor);
}

/* ------------------------------------------------------------------------------------------------------------------------------------ */
void InitializeLogFactor(int maximum)
{
	log_factor[0] = 0.0;
	for (int i = 1; i < maximum; ++i)
	       log_factor[i] = log_factor[i-1] + log((long double)i);	
}

/* ------------------------------------------------------------------------------------------------------------------------------------ */
long double CalculatePValue(const PPINetwork &ppi_network, int u, int v)
{
	int d1 = (int)ppi_network.vertex_degree[u];
	int d2 = (int)ppi_network.vertex_degree[v];
	int i = ppi_network.adjacency_matrix[u][v];
	int m = ppi_network.num_common_neighbor[u][v];
	int N = ppi_network.total_protein;

	long double value = 0.0;
	for (int j = i; j <= 1; ++j)
		for (int k = m; k <= min(d1, d2)-j; ++k)
			value += exp(LogCombination(N-2, k) + LogCombination(N-2-k, d1-k-j) + LogCombination(N-2+j-d1, d2-j-k) -
			       LogCombination(N-2, d1-1) - LogCombination(N-2, d2-1) - log(1 + (N-1-d1) * (N-1-d2) / (d1*d2)));

	return log(value);
}

long double EvaluatePPIPValue(const PPINetwork &ppi_network, int i, int j)
{
	return (i != j ? CalculatePValue(ppi_network, i, j) : 0.0);
}

bool PValueCompare(const ValueInfo &x, const ValueInfo &y)
{
	if (x.value == y.value)
		return x.index < y.index;
	return x.value < y.value;
}

void InitSingleCore(int n, vector< vector<int> > &cores, int *core_id)
{
	cores.resize(0);
	for (int i = 0; i < n; ++i)
	{
		vector<int> single(1);
		single[0] = i;
		cores.push_back(single);
		core_id[i] = i;
	}
}

void ExpandCores(const PPINetwork &ppi_network, int core_size, vector< vector<int> > &cores, int *core_id, ValueInfo **pvalue_info)
{
	bool *flag = new bool[ppi_network.total_protein];
	for (int i = 0; i < (int)cores.size(); ++i)
	{
		if (cores[i].size() > 0)
		{
			fill_n(flag, ppi_network.total_protein, false);
			flag[cores[i][0]] = true;
			for (int j = 0; j < core_size-1; ++j)
				flag[pvalue_info[cores[i][0]][j].index] = true;

			bool status = true;
			for (int j = 0; j < core_size-1 && status; ++j)
			{
				int pid = pvalue_info[cores[i][0]][j].index;
				for (int k = 0; k < core_size-1 && status; ++k)
					if (!flag[pvalue_info[pid][k].index])
					       status = false;
			}
			if (status)
			{
				for (int j = (int)cores[i].size()-1; j < core_size-1; ++j)
				{
					int pid = pvalue_info[cores[i][0]][j].index;
					int cid = core_id[pid];
					if (cores[cid].size() > 0)
						cores[cid].resize(0);
					cores[i].push_back(pid);
					core_id[pid] = core_id[i];
				}
			}
		}
	}
	delete[] flag;
}

void GenerateCores(const PPINetwork &ppi_network, vector< vector<int> > &cores)
{
	const int MaxCoreSize = min(10, ppi_network.total_protein);
	int *core_id = new int[ppi_network.total_protein];

	ValueInfo **pvalue_info = new ValueInfo *[ppi_network.total_protein];
	for (int i = 0; i < ppi_network.total_protein; ++i)
		pvalue_info[i] = new ValueInfo[ppi_network.total_protein];

	for (int i = 0; i < ppi_network.total_protein; ++i)
	{
		for (int j = 0; j < ppi_network.total_protein; ++j)
		{
			pvalue_info[i][j].index = j;
			pvalue_info[i][j].value = EvaluatePPIPValue(ppi_network, i, j);
		}
		sort(pvalue_info[i], pvalue_info[i]+ppi_network.total_protein, PValueCompare);
	}

	InitSingleCore(ppi_network.total_protein, cores, core_id);
	for (int i = 2; i <= MaxCoreSize; ++i)
	{
		printf("\t\tExpand Core (Size = %d)\n", i);
		ExpandCores(ppi_network, i, cores, core_id, pvalue_info);
	}

	for (int i = 0; i < ppi_network.total_protein; ++i)
		delete[] pvalue_info[i];
	delete[] pvalue_info;

	delete[] core_id;
}

/* ------------------------------------------------------------------------------------------------------------------------------------ */
vector<int> GetAllNeighbors(const vector<int> &vertex, vector<int> *neighbor_list, int max_id)
{
	vector<int> v;

	bool *flag = new bool[max_id];
	fill_n(flag, max_id, false);
	Flagging(flag, vertex);
	v.resize(0);
	for (int i = 0; i < (int)vertex.size(); ++i)
	{
		for (int j = 0; j < (int)neighbor_list[vertex[i]].size(); ++j)
		{
			if (!flag[neighbor_list[vertex[i]][j]])	
			{
				v.push_back(neighbor_list[vertex[i]][j]);
				flag[neighbor_list[vertex[i]][j]] = true;	
			}
		}
	}	

	delete[] flag;

	return v;
}

void AddAttachment(const PPINetwork &ppi_network, const vector<int> &core, vector<int> &attachment)
{
	if (core.size() <= 0)
	{
		return;
	}

	if (core.size() == 1)
	{
		attachment = ppi_network.neighbor_list[core[0]];
		return;
	}

	attachment.resize(0);
	vector<int> neighbor_union = GetAllNeighbors(core, ppi_network.neighbor_list, ppi_network.total_protein);
	for (int i = 0; i < (int)neighbor_union.size(); ++i)
	{
		int total = 0;
		for (int j = 0; j < (int)core.size(); ++j)
			total += ppi_network.adjacency_matrix[neighbor_union[i]][core[j]];
		if (total + total > (int)core.size())
			attachment.push_back(neighbor_union[i]);
	}
}

void GenerateAttachments(const PPINetwork &ppi_network, vector< vector<int> > &cores, vector< vector<int> > &attachments)
{ 
	attachments.resize(ppi_network.total_protein);
	for (int i = 0; i < (int)ppi_network.total_protein; ++i)
		attachments[i].resize(0);

	bool *flag = new bool[ppi_network.total_protein];
	fill_n(flag, ppi_network.total_protein, false);
	for (int i = 0; i < (int)ppi_network.total_protein; ++i)
	{
		if (cores[i].size() >= 2)
		{
			AddAttachment(ppi_network, cores[i], attachments[i]);
			vector<int> neighbor_union = GetAllNeighbors(cores[i], ppi_network.neighbor_list, ppi_network.total_protein);
			Flagging(flag, neighbor_union);
		}
	}

	for (int i = 0; i < (int)ppi_network.total_protein; ++i)
	{
		if (cores[i].size() == 1)
		{
			if (flag[cores[i][0]])
				cores[i].resize(0);
			else
				AddAttachment(ppi_network, cores[i], attachments[i]);
		}
	}

	delete[] flag;
}

/* ------------------------------------------------------------------------------------------------------------------------------------ */
int CountDegree(const vector<int> &vertex_set, vector<int> *neighbor_list)
{
	int total_degree = 0;

	for (int i = 0; i < (int)vertex_set.size(); ++i)
		total_degree += (int)neighbor_list[vertex_set[i]].size();
	
	return total_degree;
}

int CountEdge(const vector<int> &vertex_set, int **adjacency_matrix)
{
	int total_edge = 0;

	for (int i = 0; i < (int)vertex_set.size(); ++i)
		for (int j = i+1; j < (int)vertex_set.size(); ++j)
			total_edge += adjacency_matrix[vertex_set[i]][vertex_set[j]];

	return total_edge;
}

int CountEdge(const vector<int> &vertex_set1, const vector<int> &vertex_set2, int **adjacency_matrix)
{
	int total_edge = 0;

	for (int i = 0; i < (int)vertex_set1.size(); ++i)
		for (int j = 0; j < (int)vertex_set2.size(); ++j)
			total_edge += adjacency_matrix[vertex_set1[i]][vertex_set2[j]];

	return total_edge;
}

long double LogCoreProb(int ic, int m, int ia, int c, int dc, int N, int dmin)
{
	int cc = c*(c-1)/2;
	int maxsize_num, maxsize_den, index;

	maxsize_num = 0;
	for (int _ic = ic; _ic <= cc; ++_ic)
		for (int _m = m; _m <= (dc-2*_ic)/dmin; ++_m)
			for (int _ia = max(dmin*_m, ia); _ia <= min(dc - 2*_ic, _m*c); ++_ia)
				++ maxsize_num;

	long double *numerator = new long double[maxsize_num];
	index = 0;
	for (int _ic = ic; _ic <= cc; ++_ic)
		for (int _m = m; _m <= (dc-2*_ic)/dmin; ++_m)
			for (int _ia = max(dmin*_m, ia); _ia <= min(dc - 2*_ic, _m*c); ++_ia)
				numerator[index++] = LogCombination(cc, _ic) + LogCombination(N-c, _m) +
				       	LogCombination(N-c-_m, dc-2*_ic-_ia) + LogCombination(_ia-dmin*_m+_m-1, _m-1);

	maxsize_den = 0;
	for (int _ic = 0; _ic <= cc; ++_ic)
		for (int _m = 0; _m <= (dc-2*_ic)/dmin; ++_m)
			for (int _ia = dmin * _m; _ia <= min(dc - 2*_ic, _m*c); ++_ia)
				++ maxsize_den;

	long double *denominator= new long double[maxsize_den];
	index = 0;
	for (int _ic = 0; _ic <= cc; ++_ic)
		for (int _m = 0; _m <= (dc-2*_ic)/dmin; ++_m)
			for (int _ia = dmin * _m; _ia <= min(dc - 2*_ic, _m*c); ++_ia)
				denominator[index++] = LogCombination(cc, _ic) + LogCombination(N-c, _m) + 
					LogCombination(N-c-_m, dc-2*_ic-_ia) + LogCombination(_ia-dmin*_m +_m-1, _m-1);

	long double mean_num = 0.0;
	for (int i = 0; i < maxsize_num; ++i)
		mean_num += numerator[i];
	mean_num /= maxsize_num;

	long double mean_den = 0.0;
	for (int i = 0; i < maxsize_den; ++i)
		mean_den += denominator[i];
	mean_den /= maxsize_den;

	long double num = 0.0;
	for (int i = 0; i < maxsize_num; ++i)
		num += exp(numerator[i] - mean_num);
	num = log(num) + mean_num;

	long double den = 0.0;
	for (int i = 0; i < maxsize_den; ++i)
		den += exp(denominator[i] - mean_den);
	den = log(den) + mean_den;

	delete[] numerator;
	delete[] denominator;
	return num - den;
}

bool CoreScoreCompare(const ValueInfo &x, const ValueInfo &y)
{
	return x.value < y.value;
}

bool CheckComplex(int u, vector< vector<int> > &cores, vector< vector<int> > &attachments, bool *core_flag, bool *attachment_flag)
{
	if (cores[u].size() == 0)
		return false;

	for (int i = 0; i < (int)cores[u].size(); ++i)
		if (!attachment_flag[cores[u][i]])
			return true;
	return false;
}

void RemoveNoisyCores(const PPINetwork &ppi_network, vector< vector<int> > &cores, vector< vector<int> > &attachments)
{
	ValueInfo *core_score_info = new ValueInfo[ppi_network.total_protein];

	for (int i = 0; i < ppi_network.total_protein; ++i)
	{
		if (cores[i].size() >= 2)
		{
			int c = (int)cores[i].size();
			int dc = CountDegree(cores[i], ppi_network.neighbor_list);
			int ic = CountEdge(cores[i], ppi_network.adjacency_matrix);
			int ia = CountEdge(cores[i], attachments[i], ppi_network.adjacency_matrix);
			int m = (int)attachments[i].size();
			int dmin = (c+1) - (int)(c+1)/2;
			int N = (int)ppi_network.total_protein;
		
			core_score_info[i].value = LogCoreProb(ic, m, ia, c, dc, N, dmin);
		}
		else if (cores[i].size() >= 1)
		{
			core_score_info[i].value = log((long double)ppi_network.vertex_degree[cores[i][0]]);
		}

		core_score_info[i].index = i;
	}

	sort(core_score_info, core_score_info+ppi_network.total_protein, CoreScoreCompare);

	bool *core_flag = new bool[ppi_network.total_protein];
	bool *attachment_flag = new bool[ppi_network.total_protein];

	fill_n(core_flag, ppi_network.total_protein, false);
	fill_n(attachment_flag, ppi_network.total_protein, false);
	for (int i = 0; i < ppi_network.total_protein; ++i)
	{
		int u = core_score_info[i].index;
		if (CheckComplex(u, cores, attachments, core_flag, attachment_flag))
		{
			Flagging(core_flag, cores[u]);
			Flagging(attachment_flag, attachments[u]);
		}
		else
		{
			cores[u].resize(0);
			attachments[u].resize(0);
		}
	}

	delete[] core_flag;
	delete[] attachment_flag;

	delete[] core_score_info;
}

/* ------------------------------------------------------------------------------------------------------------------------------------ */
void ConstructComplex(const vector< vector<int> > &cores, const vector< vector<int> > &attachments, vector<Complex> &complexes)
{
	complexes.resize(0);
	for (int i = 0; i < (int)cores.size(); ++i)
	{
		if (cores[i].size() > 0)
		{
			Complex tmp;

			tmp.core = cores[i];
			tmp.attachment = attachments[i];
			tmp.whole.resize(cores[i].size() + attachments[i].size());
			merge(cores[i].begin(), cores[i].end(), attachments[i].begin(), attachments[i].end(), tmp.whole.begin());

			complexes.push_back(tmp);
		}
	}
}

/* ------------------------------------------------------------------------------------------------------------------------------------ */
long double RankScore(int i, int q, int N, int deg0)
{
	int qq = q * (q - 1) / 2;
	int I = min((2*i+deg0)/2, qq);
	long double *x = new long double[I+1];
	long double *y = new long double[I+1];

	long double max_x = 0.0;
	for (int j = i; j <= I; ++j)
	{
		x[j] = LogCombination(qq, j) + LogCombination(q*(N-q), 2*i+deg0-2*j);
		max_x = max(max_x, x[j]);
	}

	long double max_y = 0.0;
	for (int k = 0; k <= I; ++k)
	{
		y[k] = LogCombination(qq, k) + LogCombination(q*(N-q), 2*i+deg0-2*k);
		max_y = max(max_y, y[k]);
	}

	long double delta_x = max_x - 500;
	long double delta_y = max_y - 500;

	long double num = 0.0;
	for (int j = i; j <= I; ++j)
		num += exp(x[j] - delta_x);
	num = log(num) + delta_x;

	long double den = 0.0;
	for (int k = 0; k <= I; ++k)
		den += exp(y[k] - delta_y);
	den = log(den) + delta_y;

	delete[] x;
	delete[] y;

	return num - den;
}

bool ComplexCompare(const Complex &c1, const Complex &c2)
{
	return c1.score < c2.score;
}

void RankComplex(const PPINetwork &ppi_network, vector<Complex> &complexes)
{
	for (int u = 0; u < (int)complexes.size(); ++u)
	{
		int i = CountEdge(complexes[u].whole, ppi_network.adjacency_matrix);
		int q = (int)complexes[u].whole.size();
		int N = (int)ppi_network.total_protein;
		int deg0 = CountDegree(complexes[u].whole, ppi_network.neighbor_list) - 2*i;

		complexes[u].score = RankScore(i, q, N, deg0);
	}

	sort(complexes.begin(), complexes.end(), ComplexCompare);
}

/* ------------------------------------------------------------------------------------------------------------------------------------ */
void PredictComplex(const PPINetwork &ppi_network, vector<Complex> &complexes)
{
	const int MaxFactor = 10000000;
	log_factor = new long double[MaxFactor];	
	InitializeLogFactor(MaxFactor);

	printf("\tGenerating Cores ...\n");
	vector< vector<int> > raw_cores;
	GenerateCores(ppi_network, raw_cores);
	printf("\tCores Generated\n");

	vector< vector<int> > raw_attachments;
	GenerateAttachments(ppi_network, raw_cores, raw_attachments);
	printf("\tAttachments Added\n");

	RemoveNoisyCores(ppi_network, raw_cores, raw_attachments);
	printf("\tNoisy Cores Filtered\n");

	ConstructComplex(raw_cores, raw_attachments, complexes);
	printf("\tComplexes Constructed\n");

	RankComplex(ppi_network, complexes);
	printf("\tComplexes Ranked\n");

	delete[] log_factor;
}

/* ------------------------------------------------------------------------------------------------------------------------------------ */
void PrintComponent(FILE *fpo, const PPINetwork &ppi_network, const vector<int> &comp, const char *comp_title)
{
	vector<string> tmp;
	tmp.resize(0);
	for (int i = 0; i < (int)comp.size(); ++i)
		tmp.push_back(ppi_network.protein_name[comp[i]]);
	sort(tmp.begin(), tmp.end());

	fprintf(fpo, "%s:", comp_title);
	for (int i = 0; i < (int)comp.size(); ++i)
		fprintf(fpo, " %s", tmp[i].c_str());
	fprintf(fpo, "\n");
}

void PrintResult(const char *filename, const PPINetwork &ppi_network, const vector<Complex> &complexes)
{
	FILE *fpo = fopen(filename, "wr");

	for (int i = 0; i < (int)complexes.size(); ++i)
	{
		fprintf(fpo, "Complex %d\n", i+1);
		PrintComponent(fpo, ppi_network, complexes[i].core, "Core");
		PrintComponent(fpo, ppi_network, complexes[i].attachment, "Attachment");
		fprintf(fpo, "\n");
	}
	fclose(fpo);
}

/* ------------------------------------------------------------------------------------------------------------------------------------ */
void DestroyNetwork(PPINetwork &ppi_network)
{
	delete[] ppi_network.protein_name;
	for (int i = 0; i < ppi_network.total_protein; ++i)
		delete[] ppi_network.adjacency_matrix[i];
	delete[] ppi_network.adjacency_matrix;
	delete[] ppi_network.neighbor_list;
	delete[] ppi_network.vertex_degree;
	for (int i = 0; i < ppi_network.total_protein; ++i)
		delete[] ppi_network.num_common_neighbor[i];
	delete[] ppi_network.num_common_neighbor;
}

void ComplexPrediction(const char *infile, const char *outfile)
{
	printf("Loading PPI ...\n");
	vector<InteractionRaw> raw_interactions;
	ReadPPI(infile, raw_interactions);
	printf("PPI Loaded\n");

	printf("Constructing PPI Network ...\n");
	PPINetwork ppi_network;
	CreateNetwork(raw_interactions, ppi_network);
	printf("PPI Network Constructed\n"); 

	printf("Predicting Complexes ...\n");
	vector<Complex> complexes;
	PredictComplex(ppi_network, complexes);
	printf("Complexes Predicted\n");

	printf("Printing Result ...\n");
	PrintResult(outfile, ppi_network, complexes);
	printf("Done.\n");

	DestroyNetwork(ppi_network);
}

/* ------------------------------------------------------------------------------------------------------------------------------------ */
int main(int argc, char *argv[])
{
	if (argc == 3)
	{
		clock_t begin, end;

		begin = clock();
		ComplexPrediction(argv[1], argv[2]);
		end = clock();

		printf("Total Time: %.2lf(s)\n", (double)(end - begin)/CLOCKS_PER_SEC);
	}
	else 
	{
		printf("USAGE: Core PPI-file output-file\n");
	}

	return 0;
}
