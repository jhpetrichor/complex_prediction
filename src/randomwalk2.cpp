#include "../include/ungraph.h"
#include "../include/tools.h"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <ctime>
#include <set>
#include <algorithm>
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
        for (int step = 0; step < neighbors_n2.size() * 5; ++step) {
            // 选择下一个节点
            Protein* nextNode = choose_next_node(current_node, probability);
            if(nextNode == nullptr) {
                cout << "next_node: nullptr" << endl;
                break;
            }
            Edge* e = graph->getEdge(current_node, nextNode);
            if(e != nullptr) {
                e->visited_count += (double)1 / neighbors_n2.size();
            }

            std::cout << "Step " << step + 1 << ": Node " << current_node->protein_name << std::endl;
            if(!neighbors_n2.count(current_node)) {
//                current_node = graph->ID2Protein[startNode];
                continue;
            }

            current_node = nextNode;
        }
    }
};

int find_parent(int x, map<int, int>& parent){
    int a = x;
    while (x != parent[x])
    {
        x = parent[x];
    }
    while (a != parent[a])
    {
        int z = a;
        a = parent[a];
        parent[z] = x;
    }
    return x;
}

void split_graph(UnGraph& graph, queue<SubPPI>& ppi_queue, vector<SubPPI>& splited_ppi) {
    SubPPI current_ppi = ppi_queue.front();
    ppi_queue.pop();
    if (current_ppi.proteins.size() <= 20) {
        splited_ppi.emplace_back(current_ppi);
        return;
    } else if(current_ppi.proteins.size() > 50  && current_ppi.proteins.size() <= 80){  // 优化了对中大型蛋白质复合物的识别
        if(graph.agglomeration_coefficient(current_ppi.proteins) >= 0.3) {
            splited_ppi.emplace_back(current_ppi);
        }
    }

    map<int, int> parent;
    for (int i = 0; i < current_ppi.proteins.size(); i++) {
        parent[current_ppi.proteins[i]->id] = current_ppi.proteins[i]->id;
    }

    sort(current_ppi.edges.begin(), current_ppi.edges.end(), SubPPI::CompareByVisitedCount);
//    for(auto& e: current_ppi.edges) {
//        std::cout << e->node_a->protein_name << "\t" << e->node_b->protein_name << "\t" << e->visited_count << endl;
//    }

    int count = 0;
    for (int i = current_ppi.edges.size() - 1; i >= 0; --i) {
        Edge* edge = current_ppi.edges[i];
        int proteina = find_parent(edge->node_a->id, parent);
        int proteinb = find_parent(edge->node_b->id, parent);

        if (proteina == proteinb)
            continue;

        parent[proteina] = proteinb;
        current_ppi.set_edges.insert(edge);
        count += 1;
        if(count == current_ppi.edges.size() - 2) {
            break;
        }
    }
    int location = -1;
    while (1) {
        bool success = false;
        SubPPI new_ppi;
        set<Protein*> protein_set;

        for(int i = location + 1; i < current_ppi.proteins.size(); ++i) {
            int protein_a = current_ppi.proteins[i]->id;
            if(find_parent(protein_a, parent) == protein_a) {
                location = i;
                success = true;
                break;
            }
        }

        if (!success)
            break;

        for (int i = 0;i < current_ppi.proteins.size();i++)
        {
            if (find_parent(current_ppi.proteins[i]->id, parent) == current_ppi.proteins[location]->id)
            {
                new_ppi.proteins.push_back(current_ppi.proteins[i]);
                protein_set.insert(current_ppi.proteins[i]);
            }
        }

        vector<vector<bool>> visited(5000, vector<bool>(5000, false));

        for (Edge* edge : current_ppi.set_edges) {
            if (protein_set.count(edge->node_a) && protein_set.count(edge->node_b)) {
                if (!visited[edge->node_a->id][edge->node_b->id]) {
                    new_ppi.edges.push_back(edge);
                    visited[edge->node_a->id][edge->node_b->id] = true;
                    visited[edge->node_b->id][edge->node_a->id] = true;
                }
            }
        }

        ppi_queue.push(new_ppi);
    }
}

void write(vector<SubPPI>& splitted_ppi) {
    vector<set<Protein*>> complexes;
    for(auto& s:splitted_ppi) {
        if (s.proteins.size() >= 3) {
            set<Protein *> complex(s.proteins.begin(), s.proteins.end());
            complexes.emplace_back(complex);
        } else if (s.proteins.size() == 2) {
            set<Protein *> complex(s.proteins.begin(), s.proteins.end());
            for (auto &p: s.proteins) {
                for (auto &nei: p->neighbor) {
                    complex.insert(nei);
                }
            }
            complexes.emplace_back(complex);
        } else if (s.proteins.size() == 1) {
            set<Protein*> complex;
            complex.insert(s.proteins[0]);
            for(auto& nei: s.proteins[0]->neighbor) {
                complex.insert(nei);
            }
            complexes.emplace_back(complex);
        }
    }
    ofstream file("/home/jh/code/complex_predict/evaluation/resdd.txt");
    for(auto& complex: complexes) {
        for(auto& p: complex) {
            file << p->protein_name << "\t";
        }
        file << endl;
    }
    file.close();
}

int main() {
    UnGraph g(COLLINS_PPI);
    RandomWalk rw(&g);
    BioInformation bio;
    DAG dag;
    g.weight_by_go_term(bio, dag);
    vector<vector<double>> probability;
    // 计算平衡系数之后再游走
    g.calculate_balanced_weight();
    // 游走
    g.calculate_walk_probability(probability);
    for(auto& node: g.proteins) {
        rw.walk_N2(node->id, probability);
    }
    // 分解
    queue<SubPPI> queue_ppi;
    vector<Edge*> edges(g.edges.begin(), g.edges.end());
    vector<Protein*> proteins(g.proteins.begin(), g.proteins.end());
    SubPPI sub_ppi;
    sub_ppi.proteins = std::move(proteins);
    sub_ppi.edges = std::move(edges);
    queue_ppi.push(sub_ppi);
    vector<SubPPI> splitted_ppi;
    while(!queue_ppi.empty()) {
        split_graph(g, queue_ppi, splitted_ppi);
    }
    std::cout << splitted_ppi.size() << endl;
    for(auto& s: splitted_ppi) {
        std::cout << s.proteins.size() << "\t" << s.edges.size() << endl;
    }
    write(splitted_ppi);
    return 0;
}


