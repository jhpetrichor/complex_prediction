#include "../include/ungraph.h"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <ctime>
#include <set>
#include <algorithm>
#include "../include/tools.h"

class RandomWalk {
public:
    UnGraph* graph;
    std::unordered_map<int, std::vector<int>> adjList;
    std::mt19937 rng;   // 随机数生成
    std::map<pair<int, int>, double> walk_probability;
public:
    RandomWalk(UnGraph* g){
        graph = g;
        rng = std::mt19937(std::time(0));
        for(auto& e: graph->edges) {
            adjList[e->node_a->id].push_back(e->node_b->id);
            adjList[e->node_b->id].push_back(e->node_a->id);
        }


    }

    std::vector<int> walk(int startNode, int walkLength) {
        std::vector<int> walk;
        walk.push_back(startNode);
        int currentNode = startNode;
        for (int i = 1; i < walkLength; ++i) {
            std::vector<int> neighbors = adjList[currentNode];
            if (!neighbors.empty()) {
                std::uniform_int_distribution<int> dist(0, neighbors.size() - 1);
                int nextNode = neighbors[dist(rng)];
                walk.push_back(nextNode);
                currentNode = nextNode;
            } else {
                break; // 如果当前节点没有邻居，则停止游走
            }
        }
        return walk;
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
    std::vector<int> walk_N2(int startNode, const vector<vector<double>>& probability) {
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
            return walk;
        }
        for (int step = 0; step < neighbors_n2.size(); ++step) {
            // 选择下一个节点
            Protein* nextNode = choose_next_node(current_node, probability);
            if(nextNode == nullptr) {
                cout << "next_node: nullptr" << endl;
                break;
            }
            Edge* e = graph->getEdge(current_node, nextNode);
            if(e != nullptr) {
                e->visited_count += 1;
            }
            // 更新当前节点
            current_node = nextNode;
            // 输出结果

            std::cout << "Step " << step + 1 << ": Node " << current_node->protein_name << std::endl;
            if(!neighbors_n2.count(current_node)) break;
        }
    }
};


class CommunityDetection {
private:
    UnGraph& graph;
    RandomWalk& rw;

public:
    CommunityDetection(UnGraph& g, RandomWalk& walk) : graph(g), rw(walk) {}

    // 通过随机游走进行社区检测
    std::unordered_map<int, int> detectCommunities(int numWalks, int walkLength) {
        std::unordered_map<int, int> communities;
        std::unordered_map<int, std::vector<int>> nodeWalks;

        // 执行多次随机游走
        for (int i = 0; i < numWalks; ++i) {
            int startNode = rand() % rw.adjList.size();
            std::vector<int> walk = rw.walk(startNode, walkLength);
            for (int node : walk) {
                nodeWalks[node].push_back(i);
            }
        }

        // 根据游走路径推断社区结构
        for (const auto& pair : nodeWalks) {
            int mostCommonWalk = 0;
            int maxCount = 0;
            for (int walk : pair.second) {
                int count = std::count(nodeWalks[pair.first].begin(), nodeWalks[pair.first].end(), walk);
                if (count > maxCount) {
                    maxCount = count;
                    mostCommonWalk = walk;
                }
            }
            communities[pair.first] = mostCommonWalk;
        }

        return communities;
    }
};


int main() {
    UnGraph g;
    RandomWalk rw(&g);
    CommunityDetection cd(g, rw);
    BioInformation bio;
    DAG dag;

    g.weight_by_go_term(bio, dag);
    vector<vector<double>> probability;
    g.calculate_walk_probability(probability);
    for(auto& node: g.proteins) {
        rw.walk_N2(node->id, probability);
    }

    SubPPI subppi;
    subppi.edges = g.edges;
    vector<Protein*> proteins(g.proteins.begin(), g.proteins.end());
    subppi.proteins = proteins;
//    subppi.edges[3]->visited_count = 10;
//    subppi.edges[45]->visited_count = 100;
//    subppi.edges[43]->visited_count = 90;

    queue<SubPPI> q;
    q.push(subppi);
    vector<SubPPI> s;
    if(!q.empty()) {
        g.split_graph(g, q, s);
    }
    int i = 0;
    std::cout << "split: "<< q.size() << endl;
//    while(!q.empty()) {
//        SubPPI p = q.front();
//        q.pop();
//        std::cout << ++i << " proteins: " << p.proteins.size() << "\t edges: " << p.edges.size() << endl;
//    }
//    vector<vector<double>> probability;
//    g.calculate_walk_probability(probability);
//    std::cout << "YOL146W:" << endl;
//    Protein* node = g.ID2Protein[g.protein_name_id["YOL146W"]];
//    double sum = 0.0;
//    for(auto& neighbor: node->neighbor) {
//        std::cout << "YOL146W --> " << neighbor->protein_name << " : " << probability[node->id][neighbor->id] << "\t" << probability[neighbor->id][node->id] << endl;
//        sum += probability[node->id][neighbor->id];
//    }
//
//    rw.walk_N2(node->id, probability);

//    set<Protein*> nei_n2;
//    for(auto& n: node->neighbor) {
//        nei_n2.insert(n);
//        for(auto& n2: n->neighbor) {
//            nei_n2.insert(n2);
//        }
//    }
//    for(auto& n: nei_n2) {
//        std::cout << n->protein_name << "\t";
//    }
//    std::cout << std::endl;
//
//    vector<double> av_attractions;
//    g.calculate_average_attraction(av_attractions);
//    std::cout << node->protein_name << "\t" << av_attractions[node->id] << endl;
//    for(auto& neighbor: node->neighbor) {
//        std::cout << neighbor->protein_name << "\t" << av_attractions[neighbor->id] << endl;
//        sum += probability[node->id][neighbor->id];
//    }

    return 0;
}