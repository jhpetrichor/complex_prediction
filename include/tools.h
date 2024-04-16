#ifndef COMPLEX_PREDICT_TOOLS_H
#define COMPLEX_PREDICT_TOOLS_H

#include <iostream>
#include <set>

struct CompareSetBeySize {
    bool operator()(const std::set<std::string>& s1, const std::set<std::string>& s2) const {
        return s1.size() < s2.size();
    }
};

void write_protein_complexes_to_file(string file_path, vector<set<string>>& complexes) {
    ofstream file(file_path);
    if(!file.is_open()) {
        std::cerr << "Failed to open file! " << file_path << endl;
        exit(1);
    }

    for(auto& complex: complexes) {
        for(auto& protein: complex){
            file << protein << "\t";
        }
        file << std::endl;
    }
    file.close();
}

void write_protein_complexes_to_file(string file_path, vector<set<Protein*>>& complexes) {
    ofstream file(file_path);
    if(!file.is_open()) {
        std::cerr << "Failed to open file! " << file_path << endl;
        exit(1);
    }

    for(auto& complex: complexes) {
        for(auto& protein: complex){
            file << protein->protein_name << "\t";
        }
        file << std::endl;
    }
    file.close();
}

void read_complex(set<set<string>>& complexes) {
    fstream file(COMPLEX_FILE);
    if(!file.is_open()) {
        std::cerr << "Failed to open file! " << COMPLEX_FILE << std::endl;
        exit(1);
    }

    string line;
    while(getline(file, line)) {
        istringstream iss(line);
        string protein;
        set<string> complex;
        while(iss >> protein) {
            complex.insert(protein);
        }
        complexes.insert(complex);
    }
    file.close();
}

bool complex_mathched(const set<Protein*>& a, const set<Protein*>& b) {
    set<Protein*> common;
    set_intersection(a.begin(), a.end(), b.begin(), b.end(), inserter(common, common.begin()));
    return (double)common.size() / std::max(a.size(), b.size()) > 0.65;
}

// 大于阈值的则不应该保留
void update_complexes(vector<set<Protein*>>& complexes, set<Protein*>& complex) {
    for(auto& c: complexes) {
        if(complex_mathched(c, complex)) {
            return;
        }
    }
    complexes.emplace_back(std::move(complex));
}


#endif //COMPLEX_PREDICT_TOOLS_H
