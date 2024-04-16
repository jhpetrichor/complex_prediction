#include "../include/ungraph.h"
#include "../include/bops.h"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <ctime>
#include <set>
#include <algorithm>
#include "../include/tools.h"
#include <cstdlib>
#include <thread>
#include <mutex>

int parent[10000];

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

//    std::vector<int> walk(int startNode, int walkLength) {
//        std::vector<int> walk;
//        walk.push_back(startNode);
//        int currentNode = startNode;
//        for (int i = 1; i < walkLength; ++i) {
//            std::vector<int> neighbors = adjList[currentNode];
//            if (!neighbors.empty()) {
//                std::uniform_int_distribution<int> dist(0, neighbors.size() - 1);
//                int nextNode = neighbors[dist(rng)];
//                walk.push_back(nextNode);
//                currentNode = nextNode;
//            } else {
//                break; // 如果当前节点没有邻居，则停止游走
//            }
//        }
//        return walk;
//    }

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

//            std::cout << "Step " << step + 1 << ": Node " << current_node->protein_name << std::endl;
            if(!neighbors_n2.count(current_node)) {
                current_node = graph->ID2Protein[startNode];
                continue;
            }

            current_node = nextNode;
        }
    }
};

int find(int x)
{//This is a data structure, called Disjoint Set Union to check if two proteins are in the same set
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

void split_graph(queue<SubPPI>& ppi_queue, vector<SubPPI>& splited_ppi) {
    SubPPI current_ppi = ppi_queue.front();
    ppi_queue.pop();
    if (current_ppi.proteins.size() <= 20) {
        splited_ppi.emplace_back(current_ppi);
        return;
    }
    // Initializes disjoint set union with protein ids
    for (int i = 0; i < current_ppi.proteins.size(); i++) {
        parent[current_ppi.proteins[i]->id] = current_ppi.proteins[i]->id;
    }

    sort(current_ppi.edges.begin(), current_ppi.edges.end(), SubPPI::CompareByVisitedCount);

    int count = 0;
    for (Edge* edge : current_ppi.edges) {
        int proteina = find(edge->node_a->id);
        int proteinb = find(edge->node_b->id);

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
            if(find(protein_a) == protein_a) {
                location = i;
                success = true;
                break;
            }
        }

        if (!success)
            break;

        for (int i = 0;i < current_ppi.proteins.size();i++)
        {
            if (find(current_ppi.proteins[i]->id) == current_ppi.proteins[location]->id)
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
//        if(new_ppi.proteins.size() >= 3){
//            ppi_queue.push(new_ppi);
//            std::cout << new_ppi.proteins.size() << "\t" << new_ppi.edges.size() << endl;
//        }
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
    ofstream file("~/code/complex_predict/evaluation/resdd.txt");
    for(auto& complex: complexes) {
        for(auto& p: complex) {
            file << p->protein_name << "\t";
        }
        file << endl;
    }
    file.close();
}

void process_dataset(string ppi_file, string result_file, string temp, int i) {
    UnGraph g(ppi_file);
    RandomWalk rw(&g);
    BioInformation bio;
    DAG dag;
    g.weight_by_go_term(bio, dag);
    vector<vector<double>> probability;
    g.calculate_walk_probability(probability);
    for(auto& node: g.proteins) {
        rw.walk_N2(node->id, probability);
    }
    string temp_ppi = "/home/jh/code/complex_predict/" + temp + to_string(i);

    ofstream file(temp_ppi);
    if(!file.is_open()) {
        cout << "Failed to open file!" << endl;
    }
    for(auto& e: g.edges) {
        file << e->node_a->protein_name << "\t" << e->node_b->protein_name << "\t" << e->visited_count << endl;
    }
    file.close();
    file.flush();
    vector<Result> complexes = bops(temp_ppi);
    write_proteins(complexes, result_file + to_string(i));
}

int main() {
//    vector<string> dataset{COLLINS_PPI, GAVIN_PPI, KROGAN_CORE_PPI, KROGAN_EXTENDED_PPI, BIOGRID_PPI};
//    vector<string> result {COLLINS_RESULT, GAVIN_RESULT, KROGAN_CORE_RESULT, KROGAN_EXTENDED_RESULT, BIOGRID_RESULT};
//    vector<string> temp{"collins", "gavin", "krogan_core", "krogan_extended", "biogrid"};
//    vector<string> dataset{GAVIN_PPI};
//    vector<string> result {GAVIN_RESULT};
//    vector<string> temp{"gavin"};

//    std::vector<std::thread> threads;
//    pthread_attr_t attr;
//    pthread_attr_init(&attr);
//    size_t stackSize = 500 * 1024 * 1024;  // 设置线程堆栈大小为500MB
//    pthread_attr_setstacksize(&attr, stackSize);
//    for(int i = 0; i < dataset.size(); ++i) {
//        threads.emplace_back(process_dataset, dataset[i], result[i], temp[i]);
//    }
//    for(auto& thread: threads){
//        thread.join();
//    }
    for(int i = 0; i < 20; ++i) {
//        process_dataset(COLLINS_PPI, COLLINS_RESULT, "collins", i);
        process_dataset(BIOGRID_PPI, BIOGRID_RESULT, "biogrid", i);
    }
    return 0;
}

