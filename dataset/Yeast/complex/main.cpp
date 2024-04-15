// 整理protein complex 数据集
// 将CYC2008 和 sgd数据集合二为一
// 删除大小为1的蛋白质复合物
#include <iostream>
#include <set>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

struct CompareSetSize {
    bool operator()(set<basic_string<char>> set1, set<basic_string<char>> set2) const {
        return set1.size() > set2.size();
    }
};


void read_complex(string file_path, set<set<string>>& complexes) {
    fstream file(file_path);
    string line;
    while(getline(file, line)) {
        istringstream iss(line);
        set<string> complex;
        string protein;
        while(iss >> protein) {
            complex.insert(protein);
        }
        if(complex.size() >= 2) complexes.insert(complex);
    }
    file.close();
}

void write_complex(string file_path, vector<set<string>>& complexes) {
    ofstream  file(file_path);
    for(auto& complex: complexes) {
        for(auto& protein: complex) {
            file << protein << "\t";
        }
        file << endl;
    }
}

int main(){
    set<set<string>> complexes;
    read_complex("./gavin.txt", complexes);
    read_complex("./collins.txt", complexes);
    read_complex("./krogan_core.txt", complexes);
    read_complex("./krogan_extented.txt", complexes);
    std::vector<std::set<string>> sortedSets(complexes.begin(), complexes.end());
    std::sort(sortedSets.begin(), sortedSets.end(), CompareSetSize());
    write_complex("./my_complex.txt", sortedSets);

    return 0;
}