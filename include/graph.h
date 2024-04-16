#ifndef GRAPH_H
#define GRAPH_H

#include "bio_information.h"
#include <cmath>
#include <iostream>
#include <queue>
#include <fstream>
using namespace std;
// const int MAXN = 20010,MAXP = 20;

#define MAXN 20010
#define MAXP 20

struct PPI;

extern map<int, string> Id2protein;
extern map<string, int> Protein2id;
extern int Protein_count;

struct Interaction {
    int proteina, proteinb;
    double interaction;
    double balanced_interaction;

    Interaction(int _proteina = 0, int _proteinb = 0, double _interaction = 0, double _balanced_interaction = 0);
    bool operator<(const Interaction&) const;
};

struct PPI {
    vector<int> protein;
    vector<Interaction> interaction;
    bool mp[MAXP][MAXP]; 
};

struct Result {
    double cohesion;
    vector<int> protein;

    bool operator<(const Result& b) const;
};

// add new
void delete_interaction_by_pearson(PPI&, BioInformation&);
void jiaquan2(PPI&, BioInformation&);
void weight_interaction(PPI&, BioInformation&);
void calculate_weight(PPI&);
void get_adjacency_list(const PPI& Current_ppi, map<int, set<int>>& adjacency); 
void get_hocn(PPI&, const map<int, set<int>>&, vector<vector<double>>&);

void get_common_number(PPI&, const map<int, set<int>>&, vector<vector<int>>&);
void get_quadrangle_weight(PPI&, const map<int, set<int>>&, vector<vector<double>>&);
void max_sms(PPI&, const map<int, set<int>>&, vector<vector<int>>&);

// add new

// Reading weighted-PPI from Ppidata_file
void read_proteins_with_weight(PPI&, string);
// Reading PPI from Ppidata_file
void read_proteins(PPI&, string);
void calculate_weight(PPI&);
// This function calculates the balanced weight of interaction
void get_balanced_interaction(PPI&, double);
// This is a data structure, called Disjoint Set Union to check if two proteins are in the same set
int getfa(int, int fa[]);
// This function is to split huge ppi into small ppi
// 分解较大的 PPI 为小型 PPI
// 这里没有递归实现
void split_ppi(queue<PPI>&, vector<PPI>&);
// caculate the similarity of two complexs
double calculate_similarity(vector<int>&, vector<int>&);
// This function caculate the cohesion of one complex
double calculate_complex_cohesion(PPI&, vector<int>&);
// This function determines whether the connected subset is a protein complex and removes protein complexes that are too similar.
void update_result(Result&, vector<Result>&, double);
// This function gets the corresponding serial number of the protein in the current PPI
int get_Current_ppi_protein(PPI&, int);
// enumerates the connected subset of each small PPIN
void get_complexs(PPI&, vector<Result>&, double);
//This function divides the PPIN into smaller PPIN and calls the function to solve each smaller PPIN
vector <Result> get_result(PPI&, double);
// This function writes predicted protein complexes in result
void write_proteins(vector<Result>, string);
void print_information(string, string, double);

#endif
