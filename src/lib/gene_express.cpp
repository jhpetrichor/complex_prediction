#include "gene_express.h"


#include <fstream>
#include <istream>
#include <algorithm>

GeneExpress::GeneExpress(std::string file_path) {
    read_gene_express(file_path);
    calculate_message();
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

    }

    file.close();
}

void GeneExpress::calculate_message() {

}