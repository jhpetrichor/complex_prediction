#ifndef BIO_INFORMATION_H
#define BIO_INFORMATION_H

#include "config.h"
#include <map>
#include <string>
#include <vector>
#include <set>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <algorithm>

using namespace std;

class BioInformation {
public:
    std::map<std::string, std::vector<double>> go_expression;
    std::map<std::string, std::set<std::string>> go_slim;
    std::map<std::string, std::map<std::string, double>> subcelluar;
public:
    BioInformation();
    void read_gene_expression();
    void read_subcellular();
    void read_go_slim();
    double go_expression_pearson(const std::string&, const std::string&);
    double go_expression_spearman(const std::string&, const std::string&);
    double go_slim_slim(const std::string&, const std::string&);
    double go_slim_slim_jaccard(const std::string&, const std::string&);
    void print_go_slim();
};

#endif 
