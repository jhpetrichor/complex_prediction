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
    double weight;
public:
    Protein(int, std::string, double _weight = 0.0);
    void add_neighbor(Protein* protein);
    void remove_neighbor(Protein* _protein);
    int degree() const;
    static bool ProteinCompareByWeight(const Protein* p1, const Protein* p2);
};

class Edge {
public:
    explicit Edge(Edge *pEdge);

    Protein* node_a;
    Protein* node_b;
    double weight;
    double balanced_weight;
    double visited_count;   // 访问次数

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
    /**
     * @description: 从文件中读取文件
     * @param {string} ppi_file PPI连边文件
     */    
    explicit UnGraph(string ppi_file);
    /**
     * @description: 从节点和连边列表创建图
     * @param {&&} set_proteins: 蛋白质列表
     * @param {&&} list_edges: PPI连边列表 每相邻两个蛋白质为一组边 [0, 1] [2, 3]
     */    
    UnGraph(std::set<std::string>&& set_proteins, std::vector<std::string>&& list_edges);
    UnGraph() = default;
    ~UnGraph();
    void display() const;
    Edge* getEdge(const Protein* protein1, const Protein* protein2);

    double agglomeration_coefficient(const vector<Protein*>& nodes);
    void weight_by_go_term(BioInformation& bio, DAG& dag);
    void calculate_balanced_weight();

    // ewca 相关算法
    __attribute__((unused)) void calculate_structure_similarty(vector<vector<double>>&);
    void get_common_neighbor_size(vector<vector<int>>&);
    void get_JCS(vector<vector<double>>&, const vector<vector<int>>&);
    void get_CNS(vector<vector<double>>&, const vector<vector<double>>&);

    // cns相关算法
    void calculate_attraction(vector<double>& attractions);
    void calculate_average_attraction(vector<double>& av_attractions);

    // 计算转移概率
    void calculate_walk_probability(vector<vector<double>>& probability);


    // BOPS相关算法
//    static int get_fa(int fa[], int x);
    void split_graph(queue<SubPPI>& ppi_queue, vector<SubPPI>& splited_ppi);
    int find_parent(int protein, map<int, int>& parent);

    // 计算节点权重
    vector<double> calculate_protein_weight();

//    void
private:
    void read_edge_list(std::string, std::set<std::string>&, std::vector<std::string>&);
    void add_edge(Protein*, Protein*);

    static bool compare_pairs(const pair<Edge*, int>& pair1, const pair<Edge*, int>& pair);
};

struct SubPPI {
    vector<Protein*> proteins;
    vector<Edge*>     edges;    // 边对应的访问次数
    set<Edge*>        set_edges;

    static bool CompareByVisitedCount(const Edge* edge1, const Edge* edge2) {
        return edge1->visited_count < edge2->visited_count;
    }

//    void remove_edge(Protein* a, Protein* b);
};
#endif //COMPLEX_PREDICT_UDGRAPH_H
