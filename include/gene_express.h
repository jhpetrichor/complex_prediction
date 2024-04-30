#ifndef __GENE_EXPRESS_H__
#define __GENE_EXPRESS_H__

#include "config.h"

#include <iostream>
#include <map>
#include <vector>

class GeneExpress {
public:
    std::map<std::string, std::vector<double>> gene_express;
    /**
     * message: message[0]: average of expression
     *          message[1]: varience of expression
     *          message[2]: activate threshold of gene
     */
    std::map<std::string, std::vector<double>> message;

    GeneExpress(std::string file_path = GENE_EXPRESSION);
    void read_gene_express(std::string& file_path);
    void calculate_message();
    void display() const;
};


#endif