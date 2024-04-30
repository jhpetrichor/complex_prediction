#include "gene_express.h"

#include <cmath>

#include <fstream>
#include <istream>
#include <sstream>
#include <algorithm>

GeneExpress::GeneExpress(std::string file_path) {
    read_gene_express(file_path);
    calculate_message();
    display();
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

void GeneExpress::calculate_message() {
    // 计算mean and 。。。
    for(auto& it: gene_express) {
        std::vector<double> temp_message;
        // update average gene expression
        double temp_sum = 0.0;
        for(auto& express: it.second) {
            temp_sum += express;
        }
        double mean = temp_sum / it.second.size();
        temp_message.emplace_back(mean);

        // update gene expression variance
        double temp_variance = 0.0;
        for(auto& express: it.second) {
            temp_variance += std::pow(express - mean, 2);
        }
        double variance = temp_variance / (it.second.size() - 1);
        std::cout << "variance: " << variance << std::endl;
        temp_message.emplace_back(variance);
        // update activate threshold
        message.insert(std::move(std::make_pair(it.first, std::move(temp_message))));
    }
}

void GeneExpress::display() const {
    for(auto& it: message) {
        printf("%s, mean: %f, variance: %f\n", it.first.c_str(), it.second[0], it.second[1]);
    }
}
