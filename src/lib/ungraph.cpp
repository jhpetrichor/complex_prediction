#include "../../include/ungraph.h"

Protein::Protein(int idx, std::string _protein_name) {
    id = idx;
    protein_name = std::move(_protein_name);
}

void Protein::add_neighbor(Protein* _protein) {
    if(!_protein) {
        std::cerr << "Protein not founded!" << std::endl;
        return;
    }

    if(_protein == this) return;
    neighbor.insert(_protein);
}

void Protein::remove_neighbor(Protein* _protein) {
    if(!_protein) {
        std::cerr << "Protein not founded!" << std::endl;
        return;
    }
    neighbor.erase(_protein);
}

int Protein::degree() {
    return neighbor.size();
}

Edge::Edge(Protein* _node_a, Protein* _node_b, double _weight, double _balanced_weight, int _visited_count) {
    if(!_node_a || !_node_b) {
        std::cerr << "Node not found!" << std::endl;
        return;
    }
    node_a = _node_a;
    node_b = _node_b;
    weight = _weight;
    balanced_weight = _balanced_weight;
    visited_count = _visited_count;
}

bool Edge::operator<(const Edge& other) const {
    return balanced_weight < other.balanced_weight;
}

Edge::Edge(Edge *pEdge) {
    node_a = pEdge->node_a;
    node_b = pEdge->node_b;
    weight = pEdge->weight;
    balanced_weight = pEdge->balanced_weight;
    visited_count = pEdge->visited_count;
}

UnGraph::UnGraph() {
    std::set<std::string> protein_list;
    std::vector<std::string> edge_list;    // 记录边
    read_edge_list(PPI_FILE, protein_list, edge_list);

    ID2Protein.resize(protein_list.size());
    connected.resize(protein_list.size(), std::vector<bool>(protein_list.size(), false));
    
    // 构建节点
    int idx = 0;
    for(auto& p: protein_list) {
        protein_name_id[p] = idx;
        Protein* new_protein = new Protein(idx, p);
        proteins.insert(new_protein);
        Protein2ID[new_protein] = idx;
        ID2Protein[idx] = new_protein;
        idx +=1;
    }

    // 添加边
    for(int i= 0; i < edge_list.size(); i+=2) {
        Protein* protein_a = ID2Protein[protein_name_id[edge_list[i]]];
        Protein* protein_b = ID2Protein[protein_name_id[edge_list[i+1]]];
        connected[protein_a->id][protein_b->id] = true;
        protein_a->add_neighbor(protein_b);
        protein_b->add_neighbor(protein_a);
        Edge *e = new Edge(protein_a, protein_b);
        std::set<std::string> e_set = {protein_a->protein_name, protein_b->protein_name};
        Edge2ID[e_set] = edges.size();
        edges.emplace_back(e);
    }
    std::cout << "Node: " << proteins.size() << std::endl;
    std::cout << "edge: " << edges.size() << std::endl;
}

UnGraph::~UnGraph() {
    for(auto& _protein: proteins) {
        if (_protein) {
            delete _protein;
        }
    }
    for(auto& _edge: edges) {
        if(_edge){
            delete _edge;
        }
    }
}

void UnGraph::read_edge_list(std::string file_path, std::set<std::string>& proeins_list, std::vector<std::string>& edge_list) {
    std::fstream file(file_path);
    if(!file.is_open()) {
        std::cout << "Failed to open file! " << file_path;
        exit(1);
    }
    std::string line;
    while(std::getline(file, line)) {
        std::istringstream iss(line);
        std::string protein_name;
        while(iss >> protein_name) {
            proeins_list.insert(protein_name);
            edge_list.emplace_back(protein_name);
        }
    }
    file.close();
}

void UnGraph::add_edge(Protein* protein_a, Protein* protein_b) {
    if(!protein_a || !protein_b) {
        std::cerr << "Protein not founded!" << std::endl;
        return;
    }
    if(protein_a == protein_b) {
        return;
    }
    protein_a->add_neighbor(protein_b);
     protein_b->add_neighbor(protein_a);
}

Edge* UnGraph::getEdge(const Protein* protein_a, const Protein* protein_b) {
    std::set<std::string> e = {protein_a->protein_name, protein_b->protein_name};
    auto edge_it = Edge2ID.find(e);
    if(edge_it == Edge2ID.end()) {
        return nullptr;
    } else {
        return edges[edge_it->second];
    }
}

// 使用go term爲連邊添加權重
void UnGraph::weight_by_go_term(BioInformation& bio, DAG& dag) {
    std::cout << "weight_ppi_by_go_terms" << std::endl;
    for(int i = 0; i < edges.size(); ++i) {
        std::string protein_1 = edges[i]->node_a->protein_name;
        std::string protein_2 = edges[i]->node_b->protein_name;
        std::set<std::string> go_term_1 = bio.go_slim[protein_1];
        std::set<std::string> go_term_2 = bio.go_slim[protein_2];
        double weight = dag.get_similarity_protein(go_term_1, go_term_2);
        edges[i]->weight = weight;
        edges[i]->balanced_weight = weight;
    }
}

// 計算平衡後的加權係數
void UnGraph::calculate_balanced_weight() {
    map <int,double> Sum;
    for(int i = 0;i < edges.size();i++) {
        Sum[edges[i]->node_a->id] += edges[i]->weight;
        Sum[edges[i]->node_b->id] += edges[i]->weight;
    }
    for(int i = 0;i < edges.size();i++) {
        if(edges[i]->weight - 0 <= 0.0000001) {
            edges[i]->balanced_weight = 0.0;
            continue;
        }
        edges[i]->balanced_weight = pow(edges[i]->weight, BALANCED_INDEX) \
            / pow(Sum[edges[i]->node_a->id], BALANCED_INDEX - 1 )
                                    + pow(edges[i]->weight, BALANCED_INDEX) \
            / pow(Sum[edges[i]->node_b->id],BALANCED_INDEX - 1);
    }
}

void UnGraph::calculate_structure_similarty(vector<vector<double>>& ss_weight) {
    ss_weight.resize(ID2Protein.size(), vector<double>(ID2Protein.size(), 0.0));
    vector<vector<int>> common_neighbor_size;
    get_common_neighbor_size(common_neighbor_size);
    // 计算每条边的jcs
    vector<vector<double>> JCS;
    get_JCS(JCS, common_neighbor_size);
    // 计算cns
    vector<vector<double>> CNS;
    get_CNS(CNS, JCS);
    for(auto it1 = proteins.begin(); it1 != proteins.end(); ++it1) {
        for(auto it2 = std::next(it1); it2 != proteins.end(); ++it2) {
            double ss_w = (JCS[(*it1)->id][(*it2)->id] + CNS[(*it1)->id][(*it2)->id]) / (common_neighbor_size[(*it1)->id][(*it2)->id] + 1);
            ss_weight[(*it1)->id][(*it2)->id] = ss_w;
            ss_weight[(*it2)->id][(*it1)->id] = ss_w;
        }
    }
}

void UnGraph::get_common_neighbor_size(vector<vector<int>>& common_size) {
    common_size.resize(ID2Protein.size(), vector<int>(ID2Protein.size(), 0));
    for(auto it1 = proteins.begin(); it1 != proteins.end(); ++it1) {
        for(auto it2 = std::next(it1); it2 != proteins.end(); ++it2) {
            vector<Protein*> common;
            set_intersection((*it1)->neighbor.begin(), (*it1)->neighbor.end(),
                             (*it2)->neighbor.begin(), (*it2)->neighbor.end(),
                             inserter(common, common.begin()));
            common_size[(*it1)->id][(*it2)->id] = common.size();
            common_size[(*it2)->id][(*it1)->id] = common.size();
        }
    }
}

void UnGraph::get_JCS(vector<vector<double>>& JCS, const vector<vector<int>>& common_size) {
    JCS.resize(ID2Protein.size(), vector<double>(ID2Protein.size(), 0));
    for(auto it1 = proteins.begin(); it1 != proteins.end(); ++it1) {
        for(auto it2 = std::next(it1); it2 != proteins.end(); ++it2) {
            double jcs = (double)common_size[(*it1)->id][(*it2)->id] / ((*it1)->neighbor.size() + (*it2)->neighbor.size() - common_size[(*it1)->id][(*it2)->id]);
            JCS[(*it1)->id][(*it2)->id] = jcs;
            JCS[(*it2)->id][(*it1)->id] = jcs;
        }
    }
}

void UnGraph::get_CNS(vector<vector<double>>& CNS, const vector<vector<double>>& JCS) {
    CNS.resize(ID2Protein.size(), vector<double>(ID2Protein.size(), 0));
    for(auto it1 = proteins.begin(); it1 != proteins.end(); ++it1) {
        for(auto it2 = std::next(it1); it2 != proteins.end(); ++it2) {
            vector<Protein*> common;
            set_intersection((*it1)->neighbor.begin(), (*it1)->neighbor.end(),
                             (*it2)->neighbor.begin(), (*it2)->neighbor.end(),
                             inserter(common, common.begin()));
            double cns = 0.0;
            for(auto* neighbor: common) {
                cns += JCS[(*it1)->id][neighbor->id] + JCS[(*it2)->id][neighbor->id];
            }
            CNS[(*it1)->id][(*it2)->id] = cns;
            CNS[(*it2)->id][(*it1)->id] = cns;
        }
    }
}

// 高度节点必然大
void UnGraph::calculate_attraction(vector<double>& attractions) {
    attractions.resize( ID2Protein.size(), 0.0);

    vector<double> weight_node(ID2Protein.size(), 0.0);   // 每个节点对应的质量
    for(auto& e: edges) {
        weight_node[e->node_b->id] += e->balanced_weight;
        weight_node[e->node_a->id] += e->balanced_weight;
    }
    // 计算每一条边上面的吸引力GMm / (1 - sim)
    for(auto& e: edges) {
        double d = 2.0 - e->balanced_weight;  // 相异度   避免出现除以零的情况
        double attraction = weight_node[e->node_a->id] * weight_node[e->node_b->id] / pow(d, 2);
        attractions[e->node_a->id] += attraction;
        attractions[e->node_b->id] += attraction;
    }
}

void UnGraph::calculate_average_attraction(vector<double>& attractions) {
    attractions.resize( ID2Protein.size(), 0.0);

    vector<double> weight_node(ID2Protein.size(), 0.0);   // 每个节点对应的质量
    for(auto& e: edges) {
        weight_node[e->node_b->id] += e->balanced_weight;
        weight_node[e->node_a->id] += e->balanced_weight;
    }
    // 计算每一条边上面的吸引力GMm / (1 - sim)
    for(auto& e: edges) {
        double d = 1.1 - e->balanced_weight;  // 相异度   避免出现除以零的情况
        double attraction = weight_node[e->node_a->id] * weight_node[e->node_b->id] / pow(d, 2);
        attractions[e->node_a->id] += attraction;
        attractions[e->node_b->id] += attraction;
    }

    for(int i = 0; i < attractions.size(); ++i) {
        attractions[i] /= ID2Protein[i]->neighbor.size();
//        std::cout << ID2Protein[i]->protein_name << "\t" << attractions[i] / ID2Protein[i]->neighbor.size() << "\t" << attractions[i] << endl;
    }

//    for(int i = 0; i <attractions.size(); ++i) {
//        std::cout << ID2Protein[i]->protein_name << "\t" << attractions[i] << endl;
//    }
}

void UnGraph::calculate_walk_probability(vector<vector<double>>& probability) {
    probability.resize(ID2Protein.size(), vector<double>(ID2Protein.size(), 0.0));

    // 计算每个节点的权重
    vector<double> node_weight(ID2Protein.size(), 0.0);
    for(int i = 0; i < edges.size(); ++i) {
        node_weight[edges[i]->node_a->id] += edges[i]->balanced_weight;
        node_weight[edges[i]->node_b->id] += edges[i]->balanced_weight;
    }
    for(auto& e: edges) {
        // a ---> b 的概率
        probability[e->node_a->id][e->node_b->id] = e->balanced_weight / node_weight[e->node_a->id];
        probability[e->node_b->id][e->node_a->id] = e->balanced_weight / node_weight[e->node_b->id];
    }

}

// 其实可以用map<int, int>代替
// 将与X联通的节点直接与根节点相连
int UnGraph::get_fa(int fa[], int x) {
    int temp = x;
    while(x != fa[x]) {
        x = fa[x];
    }

    while( temp != fa[temp]) {
        int z = temp;
        temp = fa[temp];
        fa[z] = x;
    }

    return x;
}

void UnGraph::split_graph(const UnGraph& origin_ppi, queue<SubPPI>& ppi_queue, vector<SubPPI>& splited_ppi) {
    SubPPI current_ppi = ppi_queue.front();
    ppi_queue.pop();

    if(current_ppi.proteins.size() <= COMPLEX_MAX_SIZE) {
        splited_ppi.emplace_back(current_ppi);
        return;
    }

    int fa[origin_ppi.proteins.size()];
    for(int i = 0; i <= origin_ppi.proteins.size(); ++i) {
        fa[i] = i;
    }

    int add_edge_count = 0;  // 记录添加的边的数量
    std::sort(current_ppi.edges.begin(), current_ppi.edges.end(), SubPPI::CompareByVisitedCount);
//    for(auto& e: current_ppi.edges) {
//        cout << e->node_a->protein_name << "\t" << e->node_b->protein_name << "\t" << e->visited_count << endl;
//    }
    for(int i = current_ppi.edges.size() - 1; i >= 0; --i) {
        if(get_fa(fa, current_ppi.edges[i]->node_a->id) == get_fa(fa, current_ppi.edges[i]->node_b->id)) {
            continue;
        }
        fa[get_fa(fa, current_ppi.edges[i]->node_a->id)] = fa[get_fa(fa, current_ppi.edges[i]->node_b->id)];
        add_edge_count += 1;

        if(add_edge_count == current_ppi.edges.size() - 2) {
            break;
        }
    }

    int current_location = 0;
    while(true) {
        bool success = false;
        SubPPI new_ppi;
        set<Protein*> proteins;

        for (int i = current_location + 1; i <current_ppi.proteins.size(); i++) {
            // 找到一个根节点为自身的蛋白质，并记录该节点的id为location
            if (get_fa(fa, current_ppi.proteins[i]->id) == current_ppi.proteins[i]->id) {
                current_location = i;
                success = true;
                break;
            }
        }

        if (success == false)
            break;

        for (int i = 0; i < current_ppi.proteins.size(); ++i) {
            if (get_fa(fa, current_ppi.proteins[i]->id) == current_ppi.proteins[current_location]->id) {
                new_ppi.proteins.push_back(current_ppi.proteins[i]);
                proteins.insert(current_ppi.proteins[i]);
            }
        }

        for (int i = 0;i < current_ppi.edges.size();i++) {
            if (proteins.count(current_ppi.edges[i]->node_a) && proteins.count(current_ppi.edges[i]->node_b)) {
                new_ppi.edges.push_back(current_ppi.edges[i]);
            }
        }
        ppi_queue.push(new_ppi);
    }
}

bool UnGraph::compare_pairs(const pair<Edge*, int>& pair1, const pair<Edge*, int>& pair2) {
    return pair1.second < pair2.second;
}


//int main() {
//    UnGraph ug;
//    for(auto& node: ug.proteins) {
//        std::cout << node->protein_name << node->id << std::endl;
//    }
//    std::cout << ug.proteins.size() << std::endl;
//
//
//    return 0;
//}