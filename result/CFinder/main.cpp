#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

int main(int argc, char** argv) {
    vector<vector<string>> complexes;
    fstream in_file(argv[1]);
    if(!in_file.is_open()) {
        cout << "error opening file! " << endl;
        exit(1);
    }
    string line;
    while(getline(in_file, line)) {
        if(line[0] == '#') continue;
        istringstream  iss(line);
        string protein;
        vector<string> complex;
        while(iss >> protein) {
            complex.emplace_back(protein);
        }
        complexes.emplace_back(complex);
    }
    // write
    ofstream out_file(argv[2]);
    if(!out_file.is_open()) {
        std::cout << "error opening file!" << std::endl;
        exit(1);
    }
    for(auto& complex: complexes) {
        for(int i = 1; i < complex.size(); ++i) {
            out_file << complex[i] << "\t";
        }
        out_file << endl;
    }

    return 0;
}