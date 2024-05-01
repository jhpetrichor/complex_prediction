/*
 * @Description: 
 * @Author: jh
 * @Date: 2024-04-30 17:11:48
 * @LastEditors: Please set LastEditors
 * @LastEditTime: 2024-05-01 21:38:38
 */

#include "gene_express.h"
#include "ungraph.h"

int main() {
    UnGraph g(COLLINS_PPI);
    std::cout << "ddd" << std::endl;
    GeneExpress gene_express(GENE_EXPRESSION);
    std::cout << "ddd" << std::endl;
    auto dpins = gene_express.build_dynamic_PPI(&g, DPIN_MEHTOD::THREE_SIGMA);
    std::cout << "ddd" << std::endl;
    std::cout << "dpins.size(): " << dpins.size() << std::endl;
    for(int i = 0; i < dpins.size(); ++i) {
        std::cout << "proteins: " << dpins[i].ID2Protein.size() << "\tedges: " << dpins[i].edges.size();
    }

    // std::set<string> nodes {"111", "222", "333", "444", "555"};
    // std::vector<string> edges{"111", "222", "333", "444", "555", "4444"};
    // UnGraph g(std::move(nodes), std::move(edges));
    // g.display();


    return 0;
}