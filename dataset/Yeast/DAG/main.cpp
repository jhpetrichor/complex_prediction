#include <iostream>
#include <set>
#include <fstream>
#include <ostream>
#include <sstream>
#include <vector>

using namespace std;

void read_GO(string path, vector<vector<string>>& gos) {
    fstream file(path);
    string line;
    while(getline(file, line)) {
        vector<string> go;
        istringstream iss(line);
        string item;
        while(iss >> item) {
            go.emplace_back(item);
        }
        if(go.size() > 1) {
            gos.emplace_back(go);
        }
    }
    file.close();
}

void write(string path, vector<vector<string>>& gos) {
    ofstream file(path);
    for(auto& go: gos) {
        for(auto& item: go) {
            file << item << "\t";
        }
        file << "\n";
    }
    file.close();
}

int main(int argc, char** argv) {
    string in_path = argv[1];
    string out_put = argv[2];
    vector<vector<string>> gos;
    read_GO(in_path, gos);
    write(out_put, gos);
    return 0;
}
