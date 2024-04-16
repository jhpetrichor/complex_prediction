#ifndef __DAG_H__
#define __DAG_H__

#include "config.h"


#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <queue>
#include <vector>
#include <stack>
#include <utility>
#include <algorithm>


enum Relation {
    IS_A,
    PART_OF,
    NONE,
};

// 每一个几点指向自己的祖先
class Node {
public:
    int id;
    std::string go;
    std::set<Node*> ancestor;

    explicit Node(int id) : id(id) {}
    Node(int id, std::string go): id(id), go(std::move(go)) {}
    void add(Node*);
    void print()const;
};

class DAG {
public:
    std::map<std::pair<std::string, std::string>, double> Similarity;
    std::vector<Node*> nodes;
    std::vector<std::vector<Relation>> relation;
    std::map<std::string, int> GO2ID;
public:
    // 需要传入两个参数，第一个是用来初始化所有的GOterm
    // 并且固定所有的GOTerm对应的参数 
    DAG();
    ~DAG();
    Node* addNode(int id, std::string&);
    void addEdge(int, int, Relation);
    void add_edges(std::vector<std::vector<std::string>>&, Relation);
    void print() const;
    double similarity(const std::string&, const std::string&);
//    std::set<Node*> get_all_ancestor(Node*);
    static void get_all_ancestors(Node*, std::set<Node*>&);
    static void read_ancestor_child(const std::string&, std::vector<std::vector<std::string>>&);
    static void read_all_go_terms(const std::string&, std::set<std::string>&);
    void calculate_SValue(Node*, std::map<Node*, double>&);
    // 查找a到b的全部路径
    void find_all_paths(Node* a, Node* b, std::vector<Node*>& path, std::vector<std::vector<Node*>>& all_path);
    static void print_path(const std::vector<std::vector<Node*>>&);
    double get_similarity_go_gos_by_max(const std::string&, const std::set<std::string>&);
    double get_similarity_protein(const std::set<std::string>& gos1, const std::set<std::string>& gos2);
private:
    double similarity(int, int);
    // 定义go --> gos的相似性
};

#endif // __DAG_H__
