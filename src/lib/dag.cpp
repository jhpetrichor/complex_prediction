#include "../../include/dag.h"

void Node::add(Node* _ancestor) {
    ancestor.insert(_ancestor);
}


void Node::print() const {
    std::cout << go << ": " << ancestor.size() << " ";
    for(auto _ancestor: ancestor) {
        std::cout << _ancestor-> go << " ";
    }
    std::cout << std::endl;
}

DAG::DAG() {
    std::set<std::string> go_terms;
    read_all_go_terms(PATH_GO_TERMS, go_terms);
    std::cout << "go_terms.size " << go_terms.size() << std::endl;
    nodes.resize(go_terms.size() + 1);
    // 初始化节点
    int idx = 1;
    for(auto& g: go_terms) {
        Node* newNode = new Node(idx, g);
        GO2ID.insert(std::make_pair(g, idx));
        nodes[idx] = newNode;
        idx += 1;
    }
    relation.resize(go_terms.size() + 1, std::vector<Relation>(go_terms.size() + 1, Relation::NONE));

    // add edges
    std::vector<std::vector<std::string>> is_a;
    std::vector<std::vector<std::string>> part_of;
    read_ancestor_child(PATH_IS_A, is_a);
    read_ancestor_child(PATH_PART_OF, part_of);
    add_edges(is_a, Relation::IS_A);
    add_edges(part_of, Relation::PART_OF);
}

DAG::~DAG() {
    for(auto& node: nodes) {
        delete node;
    }
}

Node *DAG::addNode(int id, std::string& go) {
    Node* newNode = new Node(id, go);
    nodes.emplace_back(newNode);
    return newNode;
}

// 添加一条边
void DAG::addEdge(int from, int to, Relation _relation)  {
    if(from == 0 || to == 0) return;
    Node* fromNode = nodes[from];
    Node* toNode = nodes[to];

    if(fromNode && toNode) {
        fromNode->add(toNode);
        relation[from][to] = _relation;
        relation[to][from] = _relation;
    } else {
        std::cout << from << " " << to << std::endl;
        std::cerr << "Node not found!" << std::endl;
    }
}

// 批量添加
void DAG::add_edges(std::vector<std::vector<std::string>>& ancestor_child, Relation _relation) {
    for(auto& item: ancestor_child) {
        if(item.size() < 2) continue;
        for(int i = 1; i < item.size(); ++i) {
            if(GO2ID[item[i]] == 0){
                continue;
            }
            addEdge(GO2ID[item[i]], GO2ID[item[0]], _relation);
        }
    }
}

void DAG::print() const  {
    for(auto& node: nodes) {
        node->print();
    }
}

// 计算两个节点的相似度
//double DAG::similarity(int a, int b) {
//    if (a == b) return 1;
//
//    std::set<Node*> aChildren = nodes[a]->ancestor;
//    std::set<Node*> bChildren = nodes[b]->ancestor;
//
//    int commonChildren = 0;
//    for (Node* child : aChildren) {
//        if (bChildren.find(child) != bChildren.end()) {
//            commonChildren++;
//        }
//    }
//
//    int totalChildren = aChildren.size() + bChildren.size();
//    if (totalChildren == 0) return 0.0; // 避免除以零
//
//    return static_cast<double>(2 * commonChildren) / totalChildren;
//}

// 根据相互关系添加权重，更新直接祖先的s_value
void DAG::calculate_SValue(Node* node, std::map<Node*, double>& s_value){
    if(!node) {
        return;
    }
    for(auto& _ancestor: node->ancestor) {
        if(!s_value.count(_ancestor)) {
            continue;
        }
        double temp_s_value = 0.0;
        switch(relation[node->id][_ancestor->id]) {
            case Relation::IS_A:
                temp_s_value = s_value[node] * 0.8;
                break;
            case Relation::PART_OF:
                temp_s_value = s_value[node] * 0.6;
                break;
            case Relation::NONE:
                temp_s_value = 0.0;
        }
        s_value[_ancestor] = std::max(temp_s_value, s_value[_ancestor]);
    }
}

void DAG::find_all_paths(Node* a, Node* b, std::vector<Node*>& path, std::vector<std::vector<Node*>>& all_path) {
    path.emplace_back(a);
    if(a == b) {
        all_path.emplace_back(path);
    } else {
        for(Node* _ancestor: a->ancestor) {
            // 祖先节点不在当前路径中
            if(!std::count(path.begin(), path.end(), _ancestor)) {
                find_all_paths(_ancestor, b, path, all_path);
            }
        }
    }

    // 回溯时移除当前节点
    path.pop_back();
}

void DAG::print_path(const std::vector<std::vector<Node*>>& all_path) {
    for(auto& p: all_path) {
        std::cout << p[0]->go;
        for(int i = 0; i < p.size() - 1; ++i) {
            std::cout << " --> " << p[i]->go;
        }
        std::cout << " --> "<< p[p.size()-1]->go << std::endl;
    }
}

double DAG::similarity(int a, int b) {
    // 如果是同一GO term相似性为 1
    if (a == b) return 1;

    // 分别获取两个GO term的共同祖先
    Node *node_a = nodes[a];
    Node *node_b = nodes[b];
    std::set<Node *> ancestor_a;
    std::set<Node *> ancestor_b;
    get_all_ancestors(node_a, ancestor_a);
    get_all_ancestors(node_b, ancestor_b);
    std::set<Node *> common_ancestor;
    std::set_intersection(ancestor_a.begin(), ancestor_a.end(),
                          ancestor_b.begin(), ancestor_b.end(),
                          std::inserter(common_ancestor, common_ancestor.begin()));
    // 没有共同祖先
    // 即 a 不是 b 的祖先， b也不是a的祖先，同时, a, b没有共同祖先
    if (common_ancestor.empty() && !ancestor_a.count(node_b) && !ancestor_b.count(node_a)) {
        return 0.0;
    }

    std::set<Node*> all_ancestor;
    std::set_union(ancestor_a.begin(), ancestor_a.end(),
                          ancestor_b.begin(), ancestor_b.end(),
                          std::inserter(all_ancestor, all_ancestor.begin()));

    // 初始化 all_S-value
    std::map<Node*, double> all_S_value;
    for(auto& item: all_ancestor) {
        all_S_value[item] = 0.0;
    }

    // 需要迭代更新, 需分情况讨论
    if (ancestor_b.count(node_a)) { // a 是 b的祖先
        common_ancestor.insert(node_a);
        all_S_value[node_b] = 1.0;
        calculate_SValue(node_b, all_S_value);
        std::stack<Node*> stack;
        stack.push(node_b);
        if(!stack.empty()) {
            Node* current = stack.top();
            stack.pop();
            calculate_SValue(current, all_S_value);
            for(auto& current_ancestor: current->ancestor) {
                stack.push(current_ancestor);
            }
        }
        // 需要查找a--->b的所有路线
        std::vector<Node*> path;
        std::vector<std::vector<Node*>> all_path;
        find_all_paths(node_b, node_a, path, all_path);
        //计算公共
        double sum_1 = 0.0; // 记录分子
        for(auto& common_go: common_ancestor) {
            sum_1 += all_S_value[common_go];
        }
        double sum_2 = sum_1; // 分母
        for(auto& _path: all_path) {
            for(int i = 0; i < _path.size() - 1; ++i) {
                sum_2 += all_S_value[_path[i]];
            }
        }
        return sum_1 / sum_2;
    } else if (ancestor_a.count(node_b)) { // b 是 a的祖先
        common_ancestor.insert(node_b);
        all_S_value[node_a] = 1.0;
        calculate_SValue(node_a, all_S_value);
        std::stack<Node*> stack;
        stack.push(node_a);
        if(!stack.empty()) {
            Node* current = stack.top();
            stack.pop();
            calculate_SValue(current, all_S_value);
            for(auto& current_ancestor: current->ancestor) {
                stack.push(current_ancestor);
            }
        }

        // 需要查找a--->b的所有路线
        std::vector<Node*> path;
        std::vector<std::vector<Node*>> all_path;
        find_all_paths(node_a, node_b, path, all_path);
        //计算公共
        double sum_1 = 0.0; // 记录分子
        for(auto& common_go: common_ancestor) {
            sum_1 += all_S_value[common_go];
        }
        double sum_2 = sum_1; // 分母
        for(auto& _path: all_path) {
            for(int i = 0; i < _path.size() - 1; ++i) {
                sum_2 += all_S_value[_path[i]];
            }
        }
        return sum_1 / sum_2;
    } else { // a 和 b不清楚什么关系（或者说DAG中，ab没有联通）
        all_S_value[node_b] = 1.0;
        all_S_value[node_a] = 1.0;
        calculate_SValue(node_a, all_S_value);
        calculate_SValue(node_b, all_S_value);
        std::stack<Node*> stack;
        stack.push(node_a);
        stack.push(node_b);
        if(!stack.empty()) {
            Node* current = stack.top();
            stack.pop();
            calculate_SValue(current, all_S_value);
            for(auto& current_ancestor: current->ancestor) {
                stack.push(current_ancestor);
            }
        }
        double sum_1 = 0.0;
        for(auto& common_go: common_ancestor) {
            sum_1 += all_S_value[common_go];
        }
        double sum_2 = sum_1;
        std::vector<Node*> path;
        std::vector<std::vector<Node*>> all_path;
        find_all_paths(node_a, node_b, path, all_path);

        std::set<Node*> s1;   // a的祖先但不是b的祖先
        std::set<Node*> s2;  // b的祖先但不是a的祖先
        std::set_difference(ancestor_a.begin(), ancestor_a.end(),
                            ancestor_b.begin(), ancestor_b.end(),
                            std::inserter(s1, s1.begin()));
        std::set_difference(ancestor_b.begin(), ancestor_b.end(),
                            ancestor_a.begin(), ancestor_a.end(),
                            std::inserter(s2, s2.begin()));
        s1.insert(node_a);
        s2.insert(node_b);
        std::set<Node*> ss;   // 表示通往共同祖先所经历的祖先
        for(auto& ii: s1) {
            for(auto& jj: common_ancestor) {
                std::vector<Node*> path1;
                std::vector<std::vector<Node*>> all_paths1;
                find_all_paths(ii, jj, path1, all_paths1);
                for(auto& p: all_paths1) {
                    ss.insert(p.begin(), p.end() - 1);
                }
            }
        }

        for(auto& ii: s2) {
            for(auto& jj: common_ancestor) {
                std::vector<Node*> path1;
                std::vector<std::vector<Node*>> all_paths1;
                find_all_paths(ii, jj, path1, all_paths1);
                for(auto& p: all_paths1) {
                    ss.insert(p.begin(), p.end() - 1);
                }
            }
        }
        for(auto& p: ss) {
            sum_2 += all_S_value[p];
        }

        return sum_1 / sum_2;
    }
}

double DAG::similarity(const std::string& go1, const std::string& go2) {
    auto it = Similarity.find(std::pair<std::string, std::string>{go1, go2});
    if(it != Similarity.end()) {
        return it->second;
    } else {
        double sim = similarity(GO2ID[go1], GO2ID[go2]);
        Similarity.insert(std::make_pair(std::pair<std::string, std::string>{go1, go2}, sim));
        Similarity.insert(std::make_pair(std::pair<std::string, std::string>{go2, go1}, sim));
        return sim;
    }
}

void DAG::read_ancestor_child(const std::string& file_path, std::vector<std::vector<std::string>>& ancestor_child) {
    std::fstream file(file_path);
    if(!file.is_open()) {
        std::cerr << "Failed to open file! " << file_path << std::endl;
    }
    std::string line;
    while(getline(file, line)) {
        std::vector<std::string> _ancestor_child;
        std::istringstream iss(line);
        std::string go_term;
        while(iss >> go_term) {
            _ancestor_child.emplace_back(go_term);
        }
        ancestor_child.emplace_back(_ancestor_child);
    }
    file.close();
}

void DAG::read_all_go_terms(const std::string& file_path, std::set<std::string>& go_terms) {
    std::fstream file(file_path);
    if(!file.is_open()) {
        std::cerr << "Node not found!" << std::endl;
    }
    std::string line;
    while(getline(file, line)) {
        std::istringstream iss(line);
        std::string go_term;
        while(iss >> go_term) {
            go_terms.insert(go_term);
        }
    }
    file.close();
}

void DAG::get_all_ancestors(Node* node, std::set<Node*>& all_ancestor) {
    std::stack<Node*> stack;

    if (node != nullptr) {
        stack.push(node);
        while (!stack.empty()) {
            Node* current = stack.top();
            stack.pop();
            for (Node* ancestor : current->ancestor) {
                if (all_ancestor.insert(ancestor).second) { // Insert only if it's a new ancestor
                    stack.push(ancestor); // Push the ancestor onto the stack to explore its ancestors
                }
            }
        }
    }
}

void write_go_term(const DAG& dag) {
    std::ofstream file("./go_term_go.txt");
    for(auto& it: dag.GO2ID){
        file << it.first << "\t" << it.second << std::endl;
    }
}


void write_nodes(const DAG& dag) {
    std::ofstream file("./go_term_nodes.txt");
    for(auto& it: dag.nodes){
        file << it->go << "\t" << it->id << std::endl;
    }
}

// 计算两个蛋白质的相互作用, 需要传入两个蛋白质的GO term
double DAG::get_similarity_protein(const std::set<std::string>& gos1, const std::set<std::string>& gos2) {
    double sum_sim = 0.0;
    for(auto& g1:gos1) {
        sum_sim += get_similarity_go_gos_by_max(g1, gos2);
    }
    for(auto& g2: gos2) {
        sum_sim += get_similarity_go_gos_by_max(g2, gos1);
    }
    return sum_sim / (gos1.size() + gos2.size());
}

double DAG::get_similarity_go_gos_by_max(const std::string& go, const std::set<std::string>& gos) {
    double sim = 0.0;
    for(auto& g: gos) {
        double temp_sim = similarity(go, g);
        if(temp_sim > sim){
            sim = temp_sim;
        }
    }
    return sim;
}
// 到此为止， 有DAG基本完成
// 需要构造相似性计算的方法
//int main() {
//    DAG dag;
////    std::set<Node*> all_ancestor;
////    DAG::get_all_ancestors(dag.nodes[3], all_ancestor);
////    std::cout << "node3:";
////    dag.nodes[3]->print();
////    std::cout << "ancestors: " << std::endl;
////    for(auto& ances: all_ancestor) {
////        ances->print();
////    }
////    std::cout << "id: " <<dag.nodes[dag.GO2ID["GO:005085"]]->id << std::endl;
////    std::vector<Node*> path;
////    std::vector<std::vector<Node*>> allPaths;
////    dag.find_all_paths(dag.nodes[dag.GO2ID["GO:0000041"]], dag.nodes[dag.GO2ID["GO:0030001"]], path, allPaths);
////    std::cout << allPaths.size() << std::endl;
////    DAG::print_path(allPaths);
//    std::cout << "Node size: " << dag.nodes.size() << std::endl;
//    std::cout << "sim: " << dag.similarity("GO:0008556", "GO:0005267") << std::endl;
//    std::set<std::string> v1 = {"GO:1990904", "GO:0006396", "GO:0003735", "GO:0005739", "GO:0032543", "GO:0004525", "GO:0005840", "GO:0005762"};
//    std::set<std::string> v2 = {"GO:1990904", "GO:0003735", "GO:0019843", "GO:0005739", "GO:0032543", "GO:0006412", "GO:0005840", "GO:0005762"};
//    std::cout << dag.get_smilarity_protein(v1, v2);
//    return 0;
//}
