#include "../include/ungraph.h"

void read_complex(vector<set<string>>& complexes, set<string>& proteins) {
    fstream file("/home/jh/code/complex_predict/dataset/Yeast/complex/yeast_complex.txt");
    if(!file.is_open()) {
        std::cerr << "Failed to open file! " << endl;
        exit(1);
    }

    string line;
    while(std::getline(file, line)) {
        istringstream iss(line);
        set<string> complex;
        string protein;
        while(iss >> protein) {
            if(proteins.count(protein)) {
                complex.insert(protein);
            }
        }
        complexes.emplace_back(complex);
    }
}

int test(UnGraph& g, set<string>& complex) {
    int edge_count = 0;
    for(auto p1 = complex.begin(); p1 != complex.end(); p1++) {
        Protein* pp1 = g.ID2Protein[g.protein_name_id[*p1]];
        if(!pp1) continue;
        for(auto p2 = std::next(p1); p2 != complex.end(); p2++) {
            Protein* pp2 = g.ID2Protein[g.protein_name_id[*p2]];
            if(!pp2) continue;
            Edge* e = g.getEdge(pp1, pp2);
            if(e) {
                edge_count += 1;
            }
        }
    }
    double score = (double)edge_count * 2/ (complex.size() * (complex.size() - 1));
    std::cout << "edge:" << edge_count << "\t" << complex.size() << "\t" << score << endl;
    return score;
}

map<int, double> calculate_coff(UnGraph& g, vector<vector<set<string>>>& size_complex) {
    map<int, double> coff;

    for(auto& count: size_complex) {
        if(count.size() ==0 ) continue;

        double sum = 0;
        for(auto& complex: count) {
            sum += test(g, complex);
        }
        coff[count[0].size()] = (double)sum / count.size();
    }

    return std::move(coff);
}
//

int main() {
    UnGraph graph;
    set<string> proteins;
    for(auto& protein: graph.proteins) {
        proteins.insert(protein->protein_name);
    }
    vector<set<string>> complexes;
    read_complex(complexes, proteins);
    std::cout << "proteins: " << proteins.size() << endl;
    std::cout << "Complexes: " << complexes.size() << endl;

    vector<vector<set<string>>> size_complex(100);

    for(auto& complex: complexes) {
        if(complex.size() > 100) continue;
        size_complex[complex.size()].emplace_back(complex);
    }

    auto coff = calculate_coff(graph, size_complex);
    for(auto& it: coff) {
        cout << it.first << "\t" << it.second << endl;
    }


    return 0;
}