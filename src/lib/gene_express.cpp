/*
 * @Description: 
 * @Author: jh
 * @Date: 2024-04-30 17:11:48
 * @LastEditTime: 2024-05-01 23:00:41
 */
#include "gene_express.h"
#include "ungraph.h"

#include <cmath>

#include <fstream>
#include <istream>
#include <sstream>
#include <algorithm>

GeneExpress::GeneExpress(std::string file_path) {
    read_gene_express(file_path);
}

GeneExpress::~GeneExpress() {
    gene_express.clear();
}

void GeneExpress::read_gene_express(std::string& file_path) {
    std::fstream file(file_path);
    if(!file.is_open()) {
        std::cerr << "<ERROR> Failed to open file! " << file_path << std::endl;
        exit(1); 
    }
    std::string line;
    while(std::getline(file, line)) {
        std::istringstream iss(line);
        std::string protein_name;
        iss >> protein_name;
        std::vector<double> temp_express(36, 0.0);
        for(int i = 0; i < 36; ++i) {
            iss >> temp_express[i];
        }
        gene_express.insert(std::move(std::make_pair(std::move(protein_name), std::move(temp_express))));
    }

    file.close();
}

std::map<std::string, double> GeneExpress::activate_by_three_sigma() {
    std::map<std::string, double> result;
    // 计算mean and 。。。
    for(auto& it: gene_express) {
        // update average gene expression
        double temp_sum = 0.0;
        for(auto& express: it.second) {
            temp_sum += express;
        }
        double mean = temp_sum / it.second.size();

        // update gene expression variance
        double temp_variance = 0.0;
        for(auto& express: it.second) {
            temp_variance += std::pow(express - mean, 2);
        }
        double variance = temp_variance / (it.second.size() - 1);
        // std::cout << "variance: " << variance << std::endl;
        // update activate threshold
        double active_threshold = active_three_sigma(mean, variance);
        result.insert(std::make_pair(it.first, active_threshold));
    }
    return std::move(result);
}

double GeneExpress::active_three_sigma(double mean, double varience) {
    const double f = 1.0 / (1.0 + varience);
    const double s_varience = std::pow(varience, 0.5);
    return mean + 3  * s_varience * (1.0 - f);
}


// need todo()!
std::map<std::string, double> GeneExpress::active_by_top() {
    std::map<std::string, double> result;
    for(auto& it: gene_express) {
        std::vector<double> temp_expression(it.second.begin(), it.second.end());
        std::sort(temp_expression.begin(), temp_expression.end());
        // return temp_expression[10];
        result.insert(std::make_pair(it.first, temp_expression[12]));
    }

    // for(auto& it: gene_express) {
    //     std::cout << it.first << std::endl;
    //     for(auto value: it.second) {
    //         std::cout << value << "\t";
    //     }
    //     std::cout << std::endl;
    // }
    return std::move(result);
}

std::vector<UnGraph> GeneExpress::build_dynamic_PPI(const UnGraph* g, DPIN_MEHTOD method) {
    switch(method) {
        case DPIN_MEHTOD::THREE_SIGMA: {
            auto active = activate_by_three_sigma();
            return std::move(build_dynamic_PPI_by_active(g, active, 36));
            break;
        }
        case DPIN_MEHTOD::TOP: {
            auto active = active_by_top();
            return std::move(build_dynamic_PPI_by_active(g, active, 36));
            break;
        }
        case DPIN_MEHTOD::TIME:
            
            break;
    }
}

std::vector<UnGraph> GeneExpress::build_dynamic_PPI_by_active(const UnGraph* g, std::map<std::string, double>& active, int count) {
    std:vector<UnGraph> dpins(count);
    // update proteins
    std::vector<std::set<std::string>> set_proteins(count);
    std::vector<std::vector<std::string>> list_edges(count);
    for(auto& protein: g->ID2Protein) {
        auto it = gene_express.find(protein->protein_name);
        // 没有这个基因的表达信息，则在所有自网络中保留
        if(it == gene_express.end()) {
            for(int i = 0; i < count; ++i) {
                set_proteins[i].insert(protein->protein_name);
            }
            continue;
        }
        // 存在基因的表达信息，则需要筛选，满足阈值则保留
        for(int i = 0; i < count; ++i) {
            if(it->second[i] > active[it->first]) {
                set_proteins[i].insert(it->first);
            }
         }
    }

    // update edges 
    for(auto& e: g->edges) {
        // if(set_proteins.count(e->node_a->protein_name) && set_proteins.count(e->node_b->protein_name)) {

        // }
        for(int i = 0; i < count; ++i) {
            if(set_proteins[i].count(e->node_a->protein_name) && set_proteins[i].count(e->node_b->protein_name)) {
                list_edges[i].emplace_back(e->node_a->protein_name);
                list_edges[i].emplace_back(e->node_b->protein_name);
            }
        }
    }
    std::cout << "aaa" << endl;
    for(int i = 0; i < count; ++i) {
        UnGraph g(std::move(set_proteins[i]), std::move(list_edges[i]));
        std::cout << g.proteins.size() << "\t" << g.edges.size() << std::endl;
        dpins[i] = std::move(g);
    }

    return std::move(dpins);
}

std::vector<UnGraph> GeneExpress::build_dynamic_PPI_by_time(const UnGraph* g, std::map<std::string, double>& active, int count) {
    
}
