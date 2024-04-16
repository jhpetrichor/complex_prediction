#include <map>
#include <string>
#include <vector>

using namespace std;

#pragma warning (disable:4996)

void DelFile(char* szfilename);
int GetRowNum(char* szfilename);

int GenIdMap(char* szppi_filename, char* szpidmap_filename);
int LoadIds(char* szpid_filename, vector<string> *pvec_ids);
int LoadIdMap(char* szpid_filename, map<string, int> *ppid_map);

void LoadPPIScoreMatrix(char* szppi_score_filename, char* szprtn_ids_filename, double *pdppiscore_matrix);


void GetPPIs(char* szppi_score_filename, char *szoutput_filename);
void ConvertPPI(char* szrawppi_filename, char* szprtn_id_filename, char* szppi_pair_filename, char* szppi_matrix_filename);

void quasiCliques(char* szgraph_filename, double dmin_deg_ratio, int nmin_size, int nmax_size, char* szoutput_filename);

void GenCmplx(char* szppi_score_filename, double dmin_deg_ratio, int nmin_size, double dmin_vtx_overlap, double dmin_cross_score, char* szoutput_filename);

