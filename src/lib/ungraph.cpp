#include "../../include/ungraph.h"

Protein::Protein(int idx, std::string _protein_name, double _weight) {
    id = idx;
    protein_name = std::move(_protein_name);
    weight = _weight;
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

int Protein::degree() const {
    return (int)neighbor.size();
}

bool Protein::ProteinCompareByWeight(const Protein* p1, const Protein* p2){
    return p1->weight < p2->weight;
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

UnGraph::UnGraph(string ppi_file) {
    std::set<std::string> protein_list;
    std::vector<std::string> edge_list;    // 记录边
    read_edge_list(ppi_file, protein_list, edge_list);

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

// 直接算两个或者更多的节点
double UnGraph::agglomeration_coefficient(const vector<Protein*>& nodes){
    int edge_count = 0;
    for(auto node1 =  nodes.begin(); node1 != nodes.end(); node1++) {
        for(auto node2 = std::next(node1); node2 != nodes.end(); node2++) {
            if(getEdge(*node1, *node2)) {
                edge_count +=1 ;
            }
        }
    }
    return (double)edge_count * 2 / (nodes.size() * (nodes.size() - 1));
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

__attribute__((unused)) void UnGraph::calculate_structure_similarty(vector<vector<double>>& ss_weight) {
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
    }
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

int UnGraph::find_parent(int protein, map<int, int>& parent) {
    int root = protein;
    // find root
    while (root != parent[root]) {
        root = parent[root];
    }
    // path compression
    while (protein != root) {
        int next = parent[protein];
        parent[protein] = root;
        protein = next;
    }
    return root;
}

vector<double> UnGraph::calculate_protein_weight() {
    vector<double> protein_weight(ID2Protein.size(), 0.0);

    for(auto& e: edges) {
        protein_weight[e->node_a->id] += e->balanced_weight;
        protein_weight[e->node_b->id] += e->balanced_weight;
    }

    return std::move(protein_weight);
}

void UnGraph::split_graph(queue<SubPPI>& ppi_queue, vector<SubPPI>& splited_ppi) {
    SubPPI current_ppi = ppi_queue.front();
    ppi_queue.pop();
    if (current_ppi.proteins.size() <= 20) {
        splited_ppi.emplace_back(current_ppi);
        return;
    }

    map<int, int> parent;
    for(auto & protein : current_ppi.proteins) {
        parent[protein->id] = protein->id;
    }

    sort(current_ppi.edges.begin(), current_ppi.edges.end(), SubPPI::CompareByVisitedCount);

    int count = 0;
    for (Edge* edge : current_ppi.edges) {
        int proteina = UnGraph::find_parent(edge->node_a->id, parent);
        int proteinb = UnGraph::find_parent(edge->node_b->id, parent);

        if (proteina == proteinb)
            continue;

        parent[proteina] = proteinb;
        count += 1;
        if(count == current_ppi.proteins.size() - 2) {
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
            if(UnGraph::find_parent(protein_a, parent) == protein_a) {
                location = i;
                success = true;
                break;
            }
        }

        if (!success)
            break;

        for (int i = 0;i < current_ppi.proteins.size();i++)
        {
            if (UnGraph::find_parent(current_ppi.proteins[i]->id, parent) == current_ppi.proteins[location]->id)
            {
                new_ppi.proteins.push_back(current_ppi.proteins[i]);
                protein_set.insert(current_ppi.proteins[i]);
            }
        }

        for (Edge* edge : current_ppi.edges) {
            if (protein_set.count(edge->node_a) && protein_set.count(edge->node_b)) {
                new_ppi.edges.push_back(edge);
            }
        }
//        for(int i = 0; i < current_ppi.edges.size() - 1; ++i) {
//            if (protein_set.count(current_ppi.edges[i]->node_a) && protein_set.count(current_ppi.edges[i]->node_b)) {
//                new_ppi.edges.push_back(current_ppi.edges[i]);
//            }
//        }
        if(new_ppi.proteins.size() >= 3){
            ppi_queue.push(new_ppi);
            std::cout << new_ppi.proteins.size() << "\t" << new_ppi.edges.size() << endl;
        }
    }
}

bool UnGraph::compare_pairs(const pair<Edge*, int>& pair1, const pair<Edge*, int>& pair2) {
    return pair1.second > pair2.second;
}