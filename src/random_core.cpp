#include "../include/ungraph.h"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <ctime>
#include <set>
#include <algorithm>
#include "../include/tools.h"
#include <cstdlib>

class RandomWalk {
public:
    UnGraph* graph;
    std::unordered_map<int, std::vector<int>> adjList;
//    std::mt19937 rng;   // 随机数生成
    std::map<pair<int, int>, double> walk_probability;
public:
    RandomWalk(UnGraph* g){
        graph = g;
//        rng = std::mt19937(std::time(0));
        srand(static_cast<unsigned int>(time(NULL)));
        for(auto& e: graph->edges) {
            adjList[e->node_a->id].push_back(e->node_b->id);
            adjList[e->node_b->id].push_back(e->node_a->id);
        }


    }

    // 选择下一个节点
    Protein* choose_next_node(Protein* current_node, const vector<vector<double>>& probability) {
        // 获取当前节点的邻居
        vector<Protein*> neighbor(current_node->neighbor.begin(), current_node->neighbor.end());
        // 记录概率
        vector<double> p(neighbor.size(), 0.0);
        p[0] = probability[current_node->id][neighbor[0]->id];
        for(int i = 1; i < neighbor.size(); ++i) {
            p[i] = probability[current_node->id][neighbor[i]->id] + p[i-1];
        }

        // 生成随机数
        double r = static_cast<double>(std::rand()) / RAND_MAX;
        // 选择节点
        for(int i = 0; i < p.size(); ++i) {
            if(r < p[i]) {
                return neighbor[i];
            }
        }
        return neighbor[neighbor.size() - 1];
    }

    // 在邻域为2的范围内随机游走
    void walk_N2(int startNode, const vector<vector<double>>& probability) {
        std::vector<int> walk;   // 记录路径
        walk.push_back(startNode);

        // 获取二阶邻域内的蛋白质
        set<Protein*> neighbors_n2;
        for(auto& n1: graph->ID2Protein[startNode]->neighbor) {
            neighbors_n2.insert(n1);
            for(auto& n2: n1->neighbor) {
                neighbors_n2.insert(n2);
            }
        }
        // 选择下一个节点
        Protein* current_node = graph->ID2Protein[startNode];
        if(current_node == nullptr) {
            cout << "nullptr!" << endl;
            return;
        } else {
            current_node->weight += 1.0;
        }
        for(int step = 0; step < neighbors_n2.size() * 5; ++step) {
            Protein* nextNode = choose_next_node(current_node, probability);
            if(nextNode == nullptr) {
                cout << "next_node: nullptr" << endl;
                break;
            }
//            Edge* e = graph->getEdge(current_node, nextNode);
//            if(e != nullptr) {
//                e->visited_count += (double)1 / neighbors_n2.size();
//            }

//            std::cout << "Step " << step + 1 << ": Node " << current_node->protein_name << std::endl;
            if(!neighbors_n2.count(current_node)) {  // 游走到了二阶邻域之外
                current_node = graph->ID2Protein[startNode];
                continue;
            } else {
                nextNode->weight += 1;
            }

            current_node = nextNode;
        }
    }
};

bool is_safe(const Protein* protein, const UnGraph& graph, const set<Protein*>& clique) {
    if(clique.empty()) return true;
    for(auto& p: clique) {
        if(!graph.connected[p->id][protein->id]) {
            return false;
        }
    }
    return true;
}

void find_cliques_n1_helper(UnGraph& graph, std::set<Protein*>& candidate_clique, std::set<Protein*>& candidates,
                            std::set<std::set<Protein*>>& result) {
    if (candidates.empty()) {
        // 找到一个最大团
        result.insert(candidate_clique);
        return;
    }

    for (const auto& candidate : candidates) {
        // 检查 candidate 是否与 candidate_clique 中的节点都相连
        bool is_clique = true;
        for (const auto& node : candidate_clique) {
            if (graph.connected[candidate->id][node->id]) {
                is_clique = false;
                break;
            }
        }

        if (is_clique) {
            // 将 candidate 加入候选团
            candidate_clique.insert(candidate);

            // 获取新的候选节点
            std::set<Protein*> new_candidates;
            for (const auto& neighbor : candidates) {
                if (graph.getEdge(candidate, neighbor) != nullptr) {
                    new_candidates.insert(neighbor);
                }
            }

            // 递归寻找最大团
            find_cliques_n1_helper(graph, candidate_clique, new_candidates, result);

            // 回溯，将 candidate 移出候选团
            candidate_clique.erase(candidate);
        }
    }
}



std::set<std::set<Protein*>> find_cliques_n1(UnGraph& graph, Protein* protein) {
    std::set<std::set<Protein*>> result;

    // 获取节点的一阶邻居
    std::set<Protein*> first_neighbors = protein->neighbor;

    // 构建初始候选团
    std::set<Protein*> candidate_clique;
    candidate_clique.insert(protein);

    // 递归寻找最大团
    find_cliques_n1_helper(graph, candidate_clique, first_neighbors, result);

    return std::move(result);
}




//set<Protein*> find_all_neighor(set<Protein*>& clique) {
//    set<Protein*> all_neighbor;
//    for(auto& member: clique) {
//        for(auto& neighbor: member->neighbor) {
//            if(clique.count(neighbor) == 0) {
//                all_neighbor.insert(neighbor);
//            }
//        }
//    }
//    return std::move(all_neighbor);
//}
//
//void findMaxClique(Protein* protein, const UnGraph& graph, set<Protein*>& clique, set<Protein*>&  max_clique, set<set<Protein*>> cliques) {
//    if(clique.size() >= 3) {
//        cliques.insert(clique);
//    }
//    if(clique.size() > max_clique.size()) {
//        max_clique = clique;
//    }
//
//    if (is_safe(protein, graph, clique)) {
//        // 将节点v添加到当前团中
//        clique.insert(protein);
//    }
//    auto all_neighbor = find_all_neighor(clique);
//    for(auto& neighbor: all_neighbor) {
//        findMaxClique(neighbor, graph, clique, max_clique, cliques);
//        clique.erase(neighbor);
//        findMaxClique(, graph, clique);
//    }
//
//    // 不选择节点v，继续搜索下一个节点
//    findMaxClique(v + 1, graph, clique);
//}
//
//
//vector<Protein*> getMaxClique(const UnGraph& graph) {
//    vector<Protein*> clique;  // 当前团
//    vector<Protein*>  max_clique;  // 当前最大的团
//    set<vector<Protein*>> cliques;   // 过程中产生的团   // 只保留大于 3
//
//    findMaxClique(0, graph, clique, max_clique, cliques);
////    return maxClique;
//}

void write(UnGraph& graph, vector<Protein*>& ) {
    set<set<Protein*>> complexes;

}

int main() {
    UnGraph g(COLLINS_PPI);
    RandomWalk rw(&g);
    BioInformation bio;
    DAG dag;
//    for(auto& p: g.proteins) {
//        std::cout << p->protein_name << ": " << p->weight << endl;
//    }

//    set<set<Protein*>> all_cliques;
//    int j = 0;
//    for(auto& protein: g.proteins) {
//        std::cout << ++j << protein->protein_name << endl;
//        auto cliques = find_cliques_n1(g, protein);
//        for(auto& clique: cliques) {
//            all_cliques.insert(clique);
//        }
//    }
//    int i = 0;
//    for(auto& clique: all_cliques) {
//        std::cout << ++i << "\t";
//        for(auto& p: clique) {
//            std::cout << p->protein_name << "\t";
//        }
//        std::cout << endl;
//    }
    g.weight_by_go_term(bio, dag);
    vector<vector<double>> probability;
    // 计算平衡系数之后再游走
//    g.calculate_balanced_weight();
    // 游走
    g.calculate_walk_probability(probability);
    for(auto& node: g.proteins) {
        rw.walk_N2(node->id, probability);
    }
    // 分解
//    for(auto& p: g.proteins) {
//        std::cout << p->protein_name << ": " << p->weight << endl;
//    }
    auto protein_weight = g.calculate_protein_weight();
    for(auto& protein: g.proteins) {
        protein->weight = protein->weight / protein_weight[protein->id];
    }
//    for(auto& p: g.proteins) {
//        std::cout << p->protein_name << ": " << p->weight << endl;
//    }
    vector<Protein*> proteins(g.proteins.begin(), g.proteins.end());
    std::sort(proteins.begin(), proteins.end(), Protein::ProteinCompareByWeight);
    set<set<Protein*>> all_cliques;
    int j = 0;
    for(auto& protein: proteins) {
        if(j == 100) break;
        std::cout << ++j << protein->protein_name << endl;
        auto cliques = find_cliques_n1(g, protein);
        for(auto& clique: cliques) {
            all_cliques.insert(clique);
        }
    }
    int i = 0;
    for(auto& clique: all_cliques) {
        std::cout << ++i << "\t";
        for(auto& p: clique) {
            std::cout << p->protein_name << "\t";
        }
        std::cout << endl;
    }
    return 0;
}


