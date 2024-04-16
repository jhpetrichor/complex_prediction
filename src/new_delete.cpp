#include "../include/ungraph.h"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <ctime>
#include <set>
#include <algorithm>
#include "../include/tools.h"
#include <cstdlib>


int main() {
    UnGraph g(GAVIN_PPI);
    BioInformation bio;
    DAG dag;
    // 使用go-term加权
    g.weight_by_go_term(bio, dag);

    vector<set<Protein*>> temp_result;
    for(auto& protein: g.proteins) {
        set<Protein*> complex(protein->neighbor.begin(), protein->neighbor.end());
        complex.insert(protein);
        update_complexes(temp_result, complex);
    }
    std::cout << "temp_result.size: " << temp_result.size() << std::endl;
    write_protein_complexes_to_file("./res.txt", temp_result);


    return 0;
}


