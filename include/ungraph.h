#ifndef COMPLEX_PREDICT_UNGRAPH_H
#define COMPLEX_PREDICT_UNGRAPH_H

#include "config.h"
#include "bio_information.h"
#include "dag.h"

#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>

class Protein{
public:
    int id;
    std::string protein_name;
    std::set<Protein*> neighbor;   // 指向邻居蛋白质
public:
    Protein(int, std::string);
    void add_neighbor(Protein* protein);
    void remove_neighbor(Protein* _protein);
    int degree();
};

class Edge {
public:
    Edge(Edge *pEdge);

    Protein* node_a;
    Protein* node_b;
    double weight;
    double balanced_weight;
    int visited_count;   // 访问次数

    Edge(Protein* node_a, Protein* node_b, double _weight = 0.0, double _balanced_weight = 0.0, int _visited_count = 0);
    bool operator<(const Edge&)const;
};

struct SubPPI;

class UnGraph {
public:
    std::vector<Protein*> ID2Protein;
    std::map<Protein*, int> Protein2ID;
    std::map<std::string, int> protein_name_id;
    std::set<Protein*> proteins;
    std::vector<Edge*> edges;
    std::vector<std::vector<bool>> connected;    // 存储个节点是否直接相连

    std::map<std::set<std::string>, int> Edge2ID;


    UnGraph();
    ~UnGraph();
    Edge* getEdge(const Protein* protein1, const Protein* protein2);
    void weight_by_go_term(BioInformation& bio, DAG& dag);
    void calculate_balanced_weight();

    // ewca 相关算法
    void calculate_structure_similarty(vector<vector<double>>&);
    void get_common_neighbor_size(vector<vector<int>>&);
    void get_JCS(vector<vector<double>>&, const vector<vector<int>>&);
    void get_CNS(vector<vector<double>>&, const vector<vector<double>>&);

    // cns相关算法
    void calculate_attraction(vector<double>& attractions);
    void calculate_average_attraction(vector<double>& av_attractions);

    // 计算转移概率
    void calculate_walk_probability(vector<vector<double>>& probability);


    // BOPS相关算法
    static int get_fa(int fa[], int x);
    void split_graph(const UnGraph& origin_ppi, queue<SubPPI>& ppi_queue, vector<SubPPI>& splited_ppi);
//    void
private:
    void read_edge_list(std::string, std::set<std::string>&, std::vector<std::string>&);
    void add_edge(Protein*, Protein*);

    static bool compare_pairs(const pair<Edge*, int>& pair1, const pair<Edge*, int>& pair);
};

struct SubPPI {
    vector<Protein*> proteins;
    vector<Edge*>     edges;    // 边对应的访问次数

    static bool CompareByVisitedCount(const Edge* edge1, const Edge edge2) {
        return edge1->visited_count < edge2.visited_count;
    }
};

//
//struct CompareEdgeByVisitedCount{
//    bool operator()(const Edge*& e1, const Edge*& e2) const {
//        return e1 < e2;
//    }
//};


//typedef std::set<Protein*> Complex;

#endif //COMPLEX_PREDICT_UDGRAPH_H
