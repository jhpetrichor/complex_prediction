#include "../../include/graph.h"

#include <algorithm>
#include <cstdio>

using namespace std;

map<int, string> Id2protein;
map<string, int> Protein2id;
int Protein_count = 0;

Interaction::Interaction(int _proteina, int _proteinb, double _interaction, double _balanced_interaction) {
    proteina = _proteina;
    proteinb = _proteinb;
    interaction = _interaction;
    balanced_interaction = _balanced_interaction;
}

bool Interaction::operator<(const Interaction& b) const {
	return balanced_interaction < b.balanced_interaction;
}

bool Result::operator<(const Result& b) const {
	return cohesion > b.cohesion;
}

// new add
void delete_interaction_by_pearson(PPI& current_ppi, BioInformation& bio) {
    std::vector<int> delete_arr;
    for(int i = 0; i < current_ppi.interaction.size(); ++i) {
        std::string protein1 = Id2protein[current_ppi.interaction[i].proteina];
        std::string protein2 = Id2protein[current_ppi.interaction[i].proteinb];
        double gene_expression_slim = bio.go_expression_spearman(protein1, protein2);
        if(gene_expression_slim < -0.25) {
            // std::cout << "[Delete Interaction]: " << protein1 << " " << protein2 << " "<< gene_expression_slim << endl;
            delete_arr.emplace_back(i);
        }
    }

    // delete
    for(int i = delete_arr.size() - 1; i >= 0; --i) {
        current_ppi.interaction.erase(current_ppi.interaction.begin() + delete_arr[i]);
    }

    std::cout << "[delete number]: " << delete_arr.size() << std::endl;
}

void jiaquan2(PPI& current_ppi, BioInformation& bio) {
    for(int i = 0; i < current_ppi.interaction.size(); ++i) {
        std::string protein1 = Id2protein[current_ppi.interaction[i].proteina];
        std::string protein2 = Id2protein[current_ppi.interaction[i].proteinb];
        double gene_expression_slim = bio.go_expression_spearman(protein1, protein2);
        current_ppi.interaction[i].interaction += gene_expression_slim;
    }
}

void weight_interaction(PPI& current_ppi, BioInformation& bio) {
    std::cout << "weight_interaction" << endl;
    for(int i = 0; i < current_ppi.interaction.size(); ++i) {
        std::cout << i << endl;
        std::string protein1 = Id2protein[current_ppi.interaction[i].proteina];
        std::string protein2 = Id2protein[current_ppi.interaction[i].proteinb];
        std::cout << protein1 << "\t" << protein2 << std::endl;
        double weight = bio.go_slim_slim_jaccard(protein1, protein1);
        // std::cout << "weight: " << weight << endl;
        current_ppi.interaction[i].interaction = bio.go_slim_slim_jaccard(protein1, protein2);
    }
}

void get_adjacency_list(const PPI& Current_ppi, map<int, set<int>>& adjacency) {
	for(const auto interaction: Current_ppi.interaction) {
		adjacency[interaction.proteina].insert(interaction.proteinb);
		adjacency[interaction.proteinb].insert(interaction.proteina);
	}
}

void get_hocn(PPI& Current_ppi, const map<int, set<int>>& adjacency, vector<vector<double>>& hocn) {
    // 计算protein的结构相似性
    vector<vector<double>> jcs(Id2protein.size()+1, vector<double>(Id2protein.size()+1, 0.0));
    vector<vector<int>> common_protein(Id2protein.size()+1, vector<int>(Id2protein.size()+1, 0)); // 两个蛋白质的共同邻居的大小
    for(int i = 0; i < Current_ppi.interaction.size(); ++i) {
        int protein1 = Current_ppi.interaction[i].proteina;
        int protein2 = Current_ppi.interaction[i].proteinb;

        // 不可能没有
        auto protein1_nei = adjacency.find(protein1)->second;
        auto protein2_nei = adjacency.find(protein2)->second;
        set<int> common_nei;
        std::set_intersection(protein1_nei.begin(), protein1_nei.end(), 
                              protein2_nei.begin(), protein2_nei.end(),
                            std::inserter(common_nei, common_nei.begin()));
        if(common_nei.empty()) {
            std::cout << "common_nei is empty!" << std::endl;
            common_protein[protein1][protein2] = 0;
            common_protein[protein2][protein1] = 0;
            jcs[protein1][protein2] = 0.0;
            jcs[protein2][protein1] = 0.0;
            continue;
        }
        common_protein[protein1][protein2] = common_nei.size();
        common_protein[protein2][protein1] = common_nei.size();

        int union_nei_size = protein1_nei.size() + protein2_nei.size() - common_nei.size();
        double j = (double)common_nei.size() / (double)union_nei_size;
        jcs[protein1][protein2] = j;
        jcs[protein2][protein1] = j;
    }

    for(int i = 0; i < Current_ppi.interaction.size(); ++i){
        double cns_score = 0.0;
        
        int protein1 = Current_ppi.interaction[i].proteina;
        int protein2 = Current_ppi.interaction[i].proteinb;
        auto protein1_nei = adjacency.find(protein1)->second;
        auto protein2_nei = adjacency.find(protein2)->second;

        set<int> common_nei;
        std::set_intersection(protein1_nei.begin(), protein1_nei.end(), 
                              protein2_nei.begin(), protein2_nei.end(),
                            std::inserter(common_nei, common_nei.begin()));

        for(auto item: common_nei) {
            cns_score +=  jcs[protein1][item] * jcs[item][protein2];
        }

        double hocn_score = (jcs[protein1][protein2] + cns_score) / (double)(common_protein[protein1][protein2] + 1);
        hocn[protein1][protein2] = hocn_score;
        hocn[protein2][protein1] = hocn_score;       
    }
}

void get_common_number(PPI& current_ppi, const map<int, set<int>>& adjacency, vector<vector<int>>& common_number) {
    for(int i = 0; i < current_ppi.interaction.size(); ++i){
        int protein1 = current_ppi.interaction[i].proteina;
        int protein2 = current_ppi.interaction[i].proteinb;

        const set<int> nei_1 = adjacency.find(protein1)->second;
        const set<int> nei_2 = adjacency.find(protein2)->second;

        set<int> common_protein;
        set_intersection(nei_1.begin(), nei_1.end(),
                         nei_2.begin(), nei_2.end(),
                         inserter(common_protein, common_protein.begin()));

        common_number[protein1][protein2] = common_protein.size();
        common_number[protein2][protein1] = common_protein.size();
    }
}

void get_quadrangle_weight(PPI& current_ppi, const map<int, set<int>>& adjacency, vector<vector<double>>& weight) {
    std::cout << "get_common_number" << std::endl;
    vector<vector<int>> common_number(Id2protein.size() + 1, vector<int>(Id2protein.size() + 1, 0));
    get_common_number(current_ppi, adjacency, common_number);
    
    // 对每一条变加权
    std::cout << "weightting ..." << std::endl;
    for(int i = 0; i < current_ppi.interaction.size(); ++i){
        set<int> v1;
        set<int> v2; 
        
        int protein1 = current_ppi.interaction[i].proteina;
        int protein2 = current_ppi.interaction[i].proteinb;
        const set<int> nei_1 = adjacency.find(protein1)->second;
        const set<int> nei_2 = adjacency.find(protein2)->second;

        for(const auto xx: nei_1){
            for(const auto yy: nei_2) {
                if(xx == protein2 || yy == protein1) continue;
                if(xx != yy &&adjacency.find(xx)->second.count(yy)){
                    v1.insert(xx);
                    v2.insert(yy);
                }
            }
        }

        if(v1.empty() || v2.empty()) {
            weight[protein1][protein2] = 0.0;
            weight[protein2][protein1] = 0.0;
            continue;
        }
        int common_nei_number = 0;
        for(auto node: v1) {
            common_nei_number += common_number[node][protein2];
        }
        for(auto node: v2){
            common_nei_number += common_number[node][protein1];
        }

        double weight_score = (double)common_nei_number / (double)(v1.size() +v2.size());
        weight[protein1][protein2] = weight_score;
        weight[protein2][protein1] = weight_score;
    }
}

// 归一化加权
void max_sms(PPI& current_ppi, const map<int, set<int>>& adjacency, vector<vector<double>>& weight) {
    // 计算cn
    vector<vector<double>> cn(Id2protein.size() + 1, vector<double>(Id2protein.size() + 1, 0.0));
    for(int i = 1; i <= Id2protein.size(); ++i) {
        for(int j = i+1; j <= Id2protein.size(); ++j) {
            const auto nei_1 = adjacency.find(i)->second;
            const auto nei_2 = adjacency.find(j)->second;

            set<int> common_nei;
            set<int> union_nei;
            set_intersection(nei_1.begin(), nei_1.end(),
                             nei_2.begin(), nei_2.end(),
                             inserter(common_nei, common_nei.begin()));
            set_union(nei_1.begin(), nei_1.end(),
                             nei_2.begin(), nei_2.end(),
                             inserter(union_nei, union_nei.begin()));
            
            cn[i][j] = (double)common_nei.size() / (double)union_nei.size();
            cn[j][i] = (double)common_nei.size() / (double)union_nei.size();
        }
    }

    
    for(int i = 0; i < current_ppi.interaction.size(); ++i){
        double score = 0.0;
        int protein1 = current_ppi.interaction[i].proteina;
        int protein2 = current_ppi.interaction[i].proteinb;

        const set<int> nei_1 = adjacency.find(protein1)->second;
        const auto nei_2 = adjacency.find(protein2)->second;

        for(auto xx: nei_1) {
            for(auto yy: nei_2) {
                if(xx != yy && adjacency.find(xx)->second.count(yy)){
                    score = max(score, cn[xx][protein1] * cn[yy][protein2]);
                } 
            }
        }
        weight[protein1][protein2] = score;
        weight[protein2][protein1] = score;
    }
}

// Reading PPI from Ppidata_file
void read_proteins_with_weight(PPI &Current_ppi, string Ppidata_file) {
    ifstream fin(Ppidata_file);
    if(!fin.is_open()) {
        cout << "--- [ERROR]: Open PPI_file failed! ---" << Ppidata_file << endl;
        exit(1);
    }
    
    string Protein_a,Protein_b;
    double Protein_interaction;
    set <int> Protein_set;

    while(fin >> Protein_a >> Protein_b >>  Protein_interaction) {
	    if(Protein2id.count(Protein_a) == 0) {
            Protein2id[Protein_a] = ++Protein_count;
            Id2protein[Protein_count] = Protein_a;
        }

        if(Protein2id.count(Protein_b) == 0) {
            Protein2id[Protein_b] = ++Protein_count;
            Id2protein[Protein_count] = Protein_b;
        }

        if(Protein_set.find(Protein2id[Protein_a]) == Protein_set.end()) {
            Protein_set.insert(Protein2id[Protein_a]);
            Current_ppi.protein.push_back(Protein2id[Protein_a]);
        }
        if(Protein_set.find(Protein2id[Protein_b]) == Protein_set.end()) {
            Protein_set.insert(Protein2id[Protein_b]);
            Current_ppi.protein.push_back(Protein2id[Protein_b]);
        }
        Current_ppi.interaction.push_back(Interaction(Protein2id[Protein_a],Protein2id[Protein_b],Protein_interaction));
    }

    fin.close();
}

void read_proteins(PPI& Current_ppi, string Ppidata_file) {
    ifstream fin(Ppidata_file);
    if(!fin.is_open()) {
        cout << "--- [ERROR]: Open PPI_file failed! --- " << Ppidata_file << endl;
        exit(1);
    }
    
    string Protein_a,Protein_b;
    set<int> Protein_set;

    while(fin >> Protein_a >> Protein_b) {
	    if(Protein2id.count(Protein_a) == 0) {
            Protein2id[Protein_a] = ++Protein_count;
            Id2protein[Protein_count] = Protein_a;
        }

        if(Protein2id.count(Protein_b) == 0) {
            Protein2id[Protein_b] =   ++Protein_count;
            Id2protein[Protein_count] = Protein_b;
        }

        if(Protein_set.find(Protein2id[Protein_a]) == Protein_set.end()) {
            Protein_set.insert(Protein2id[Protein_a]);
            Current_ppi.protein.push_back(Protein2id[Protein_a]);
        }
        if(Protein_set.find(Protein2id[Protein_b]) == Protein_set.end()) {
            Protein_set.insert(Protein2id[Protein_b]);
            Current_ppi.protein.push_back(Protein2id[Protein_b]);
        }
        Current_ppi.interaction.push_back(Interaction(Protein2id[Protein_a],Protein2id[Protein_b]));
    }

    fin.close();
}

// This function calculates the balanced weight of interaction
void get_balanced_interaction(PPI& Current_ppi, double Balanced_index) {
    map <int,double> Sum;
	for (int i = 0;i < Current_ppi.interaction.size();i++) {    
	    Sum[Current_ppi.interaction[i].proteina] += Current_ppi.interaction[i].interaction;
	    Sum[Current_ppi.interaction[i].proteinb] += Current_ppi.interaction[i].interaction;
	}
	for (int i = 0;i < Current_ppi.interaction.size();i++) {
	    Current_ppi.interaction[i].balanced_interaction = pow(Current_ppi.interaction[i].interaction,Balanced_index) \
            / pow(Sum[Current_ppi.interaction[i].proteina],Balanced_index - 1) 
            + pow(Current_ppi.interaction[i].interaction,Balanced_index) \
            / pow(Sum[Current_ppi.interaction[i].proteinb],Balanced_index - 1);
	}
    return;
}

// This is a data structure, called Disjoint Set Union to check if two proteins are in the same set
// 会更改fa[]的值
int getfa(int x, int fa[]) {
    int a = x;
    while (x != fa[x]) {
        x = fa[x];
    }

    while (a != fa[a]) {
        int z = a;
        a = fa[a];
        fa[z] = x;
    }
    return x;
}

// This function is to split huge ppi into small ppi
// 分解较大的 PPI 为小型 PPI
// 这里没有递归实现
void split_ppi(queue <PPI>& Ppi_queue, vector<PPI>& Splitted_ppi) {
    // 这里的 Current_ppi 是待分解的 PPI, 不是原始的PPI网络
    PPI Current_ppi = Ppi_queue.front();
	Ppi_queue.pop();
	
    // Count 记录添加边的数量
    int Count = 0,location = -1;
    // 满足大小要求， 则将其暂时视为蛋白质复合物
    if (Current_ppi.protein.size() <= MAXP) {
        Splitted_ppi.push_back(Current_ppi);
        return;
    }
    
    int fa[MAXN + 1];
    // 每个节点的根节点为自身
    for (int i = 1; i <= MAXN; i++) {
        fa[i] = i;
    }

    // if Current_ppi.protein.size() > MAXP .... 大小不满去要求, 需要进一步分解
    // 从权重最大的边开始
    sort(Current_ppi.interaction.begin(), Current_ppi.interaction.end());
    for (int i = Current_ppi.interaction.size() - 1; i >= 0; i--) {
        // 属于同一派系， 继续遍历
        if (getfa(Current_ppi.interaction[i].proteina, fa) == getfa(Current_ppi.interaction[i].proteinb, fa)) {
            continue;
        }

        // 不属于同意派系则添加边
        fa[getfa(Current_ppi.interaction[i].proteina, fa)] = getfa(Current_ppi.interaction[i].proteinb, fa);
        Count++;
        
        if (Count == Current_ppi.protein.size() - 2) {
            break;
        }
    }

    // do while !success    
    while(true) {
        bool success = false;
        PPI New_ppi;
        set<int> Protein_set;
        
        for (int i = location + 1; i < Current_ppi.protein.size(); i++) {
            // 找到一个根节点为自身的蛋白质，并记录该节点的id为location
            if (getfa(Current_ppi.protein[i],fa) == Current_ppi.protein[i]) {
                location = i;
                success = true;
                break;
            }
        }

        // 所有蛋白质均被分配
        if (success == false)
            break;

        // insert protein to Protein_set
        // 讲location的派系添加到当前protein_set
        for (int i = 0; i < Current_ppi.protein.size(); i++) {
            if (getfa(Current_ppi.protein[i],fa) == Current_ppi.protein[location]) { 
                New_ppi.protein.push_back(Current_ppi.protein[i]);
                Protein_set.insert(Current_ppi.protein[i]);
            }
        }

        // add Interaction for new_PPI
        for (int i = 0;i < Current_ppi.interaction.size();i++) {
            if (Protein_set.find(Current_ppi.interaction[i].proteina) != Protein_set.end()
                 && Protein_set.find(Current_ppi.interaction[i].proteinb) != Protein_set.end()) {
                New_ppi.interaction.push_back(Current_ppi.interaction[i]);
            }
        }
        Ppi_queue.push(New_ppi);
    }
    return;
}

// caculate the similarity of two complexs
double calculate_similarity(vector<int>& Complexa, vector<int>& Complexb) {
    set<int> Complexb_set;
    int count = 0;
    
    for (int i = 0; i < Complexb.size(); i++) {
        Complexb_set.insert(Complexb[i]);
    }

    for (int i = 0; i < Complexa.size(); i++) {
        if (Complexb_set.find(Complexa[i]) != Complexb_set.end()) {
		    count += 1;
		}
    }

    return (double)count / max(Complexa.size(),Complexb.size());
}

// This function caculate the cohesion of one complex
double calculate_complex_cohesion(PPI& Current_ppi, vector<int>& Complex) {
	set<int> Complex_set;
    map<int,double> Sum;
    map<int,int> Count;
    double Cohesion = 0.0;
    
    for (int i = 0;i < Complex.size();i++) {
        Complex_set.insert(Complex[i]);
    }

    for (int i = 0;i < Current_ppi.interaction.size();i++) {
        if (Complex_set.find(Current_ppi.interaction[i].proteina) == Complex_set.end() || Complex_set.find(Current_ppi.interaction[i].proteinb) == Complex_set.end()) {
            continue;
    	}
        
        Count[Current_ppi.interaction[i].proteina] += 1;
        Count[Current_ppi.interaction[i].proteinb] += 1;
        Sum[Current_ppi.interaction[i].proteina] += Current_ppi.interaction[i].balanced_interaction;
        Sum[Current_ppi.interaction[i].proteinb] += Current_ppi.interaction[i].balanced_interaction;
    }

	for (int i = 0;i < Complex.size();i++) {
	    Cohesion += Sum[Complex[i]] * (Count[Complex[i]] + 1) / Complex.size();
	}

	Cohesion /= Complex.size();
    return Cohesion;
}


// This function determines whether the connected subset is a protein complex and removes protein complexes that are too similar.
void update_result(Result &Complex,vector<Result> &result,double Similarity_threshold) {
    // 如果两个复合物相似度太高，则不保留
    for(int i = 0;i < result.size();i++) {
        if(calculate_similarity(result[i].protein,Complex.protein) >= Similarity_threshold) {
            return;
        }
    }
    result.push_back(Complex);
    
    return;
}

// This function gets the corresponding serial number of the protein in the current PPI
int get_Current_ppi_protein(PPI &Current_ppi,int x) {
    // lower_bound 二分查找
	return lower_bound(Current_ppi.protein.begin(),Current_ppi.protein.end(),x) - Current_ppi.protein.begin();
}


// enumerates the connected subset of each small PPIN
void get_complexs(PPI& Current_ppi, vector<Result>& result, double Similarity_threshold) {
    bool Connectivity[MAXP][MAXP];    
	map<int,bool> Complex_record;    // 记录访问过的蛋白质
	queue<int> Complex_queue;        // 待处理的蛋白质组合
	memset(Connectivity, false, sizeof(Connectivity));

    // 获取邻接表
	sort(Current_ppi.protein.begin(),Current_ppi.protein.end());
	for (int i = 0;i < Current_ppi.interaction.size();i++) {
		int proteina_id = get_Current_ppi_protein(Current_ppi,Current_ppi.interaction[i].proteina);
		int proteinb_id = get_Current_ppi_protein(Current_ppi,Current_ppi.interaction[i].proteinb);
		Connectivity[proteina_id][proteinb_id] = Connectivity[proteinb_id][proteina_id] = true;
	}

	for (int i = 0;i < Current_ppi.protein.size();i++) {
		Complex_queue.push(1 << i);
		Complex_record[1 << i] = true;
	} 

	while (Complex_queue.size() != 0) {
		int Complex_state = Complex_queue.front();
		Complex_queue.pop();
		for (int i = 0;i < Current_ppi.protein.size();i++) {
			if ((Complex_state & (1 << i)) != 0) {
				for (int j = 0;j < Current_ppi.protein.size();j++) {
					if (Connectivity[i][j] == false) {
						continue;
					}

					if ((Complex_state & (1 << j)) != 0) {
						continue;
					}

					if (Complex_record.count(Complex_state | (1 << j)) != 0) {
						continue;
					}

					Complex_record[Complex_state | (1 << j)] = true;
					Complex_queue.push(Complex_state | (1 << j));
				}
			}
		}
	}
	
	vector <int> Complex;
	map <int,bool>::iterator it;
	vector <Result> Complex_result;
	it = Complex_record.begin();
	while (it != Complex_record.end()) {
		Complex.clear();
		for (int i = 0;i < Current_ppi.protein.size();i++) {
			if ((it->first & (1 << i)) != 0) { 
				Complex.push_back(Current_ppi.protein[i]);
			} 
		}
        if (Complex.size() >= 3) {
			Result tmp;
			tmp.protein = Complex;
			tmp.cohesion = calculate_complex_cohesion(Current_ppi,Complex);
			Complex_result.push_back(tmp);
		}

        it++;
	}
	
	if (Complex_result.size() == 0)
		return;
	
    sort(Complex_result.begin(),Complex_result.end());
	for (int i = 0;i < Complex_result.size();i++) {
		update_result(Complex_result[i],result,Similarity_threshold);
	}

	return;
}

//This function divides the PPIN into smaller PPIN and calls the function to solve each smaller PPIN
vector<Result> get_result(PPI &Current_ppi,double Similarity_threshold) {
    queue<PPI> Ppi_queue;
    vector<PPI> Splitting_ppi;
    
    Ppi_queue.push(Current_ppi);
    while (!Ppi_queue.empty()) {
		split_ppi(Ppi_queue,Splitting_ppi);     
    }

    vector<Result> result;
    for (int i = 0;i < Splitting_ppi.size();i++) {
        vector <Result> Temporary_res;
        get_complexs(Splitting_ppi[i],Temporary_res,Similarity_threshold); 
      
		for (int j = 0;j < Temporary_res.size();j++) {
            result.push_back(Temporary_res[j]);
        }
    }

    return result;
}

// This function writes predicted protein complexes in result
void write_proteins(vector <Result> Result_complex,string Result_file) {
    sort(Result_complex.begin(),Result_complex.end()); 
	ofstream fout(Result_file);
	int Tmpk = 0;
	string Tmps;
	    
	int Threshold_complex = min(int(Result_complex.size() * 0.05),int(Result_complex.size() - 1));
	double Cohesion_threshold = Result_complex[Threshold_complex].cohesion * 0.5;
		
    for (int i = 0;i < Result_complex.size();i++) {
        if (Result_complex[i].protein.size() <= 1) {
            continue;
        }
        if (i > Result_complex.size() / 2)
        	break;
        Tmpk += 1;
        Tmps = "";
        
        for (int j = 0;j < Result_complex[i].protein.size();j++) {
            Tmps = Tmps + " " + Id2protein[Result_complex[i].protein[j]];
        }

        fout << Tmps << endl;
    }
    fout.close();
    return;
}



void print_information(string Ppidata_file,string Result_file,double Balanced_index) {
	printf("_________________________________________\n");
	cout << "The PPI_file is "<< Ppidata_file << endl;
	cout << "The result_file is "<< Result_file << endl;
	printf("The balanced index is %.3lf\n",Balanced_index);
	printf("_________________________________________\n");
	printf("It will takes tens of minutes.\n");
}

// int main(int argc, char *argv[]) {
// 	double Balanced_index,Similarity_threshold = 0.45;
// 	string Ppidata_file,Result_file;
// 	PPI Original_ppi;
	
// 	Ppidata_file = argv[1];
// 	Result_file = argv[2];
// 	read_proteins(Original_ppi,Ppidata_file);
// 	if (argc >= 4) {
// 		sscanf(argv[3],"%lf",&Balanced_index);
// 	} else {
// 		Balanced_index = 1.5;
// 	}
	
// 	print_information(Ppidata_file,Result_file,Balanced_index);
//     get_balanced_interaction(Original_ppi,Balanced_index);
//     write_proteins(get_result(Original_ppi,Similarity_threshold),Result_file);
//     return 0;
// }

// int main(void) {
// 	double Balanced_index = 1.5,  Similarity_threshold = 0.45;
// 	string Ppidata_file,Result_file;
// 	PPI Original_ppi;
	
// 	// Ppidata_file = "./ppi.txt";
// 	Ppidata_file = "/home/jh/code/BOPS/dataset/collins.txt";
// 	Result_file = "./result_delete.txt";
// 	read_proteins_with_weight(Original_ppi,Ppidata_file);
    

//     std::cout << Id2protein.size() << std::endl;
// 	BioInformation bio;
//     std::cout << "before: " << Original_ppi.interaction.size() << std::endl;
//     weight_interaction(Original_ppi, bio);
//     // delete_interaction_by_pearson(Original_ppi, bio);
//     std::cout << "after delete: " << Original_ppi.interaction.size() << std::endl;

// 	map<int, set<int>> adjacency;
//     std:cout << "get_adiacency..." << std::endl;
// 	get_adjacency_list(Original_ppi,adjacency);
//     std::cout << "get_hocn..." << std::endl;
//     vector<vector<double>> weight(Id2protein.size()+1, vector<double>(Id2protein.size()+1, 0.0));
//     get_hocn(Original_ppi, adjacency, weight);
//     // get_quadrangle_weight(Original_ppi, adjacency, weight);
//     // max_sms(Original_ppi, adjacency, weight);
//     // double max = 0.0;
//     // for(int i = 0; i < Original_ppi.interaction.size(); ++i){
//     //     int protein1 = Original_ppi.interaction[i].proteina;
//     //     int protein2 = Original_ppi.interaction[i].proteinb;
//     //     // Original_ppi.interaction[i].interaction = Original_ppi.interaction[i].interaction + weight[protein1][protein2];
//     //     if(max < weight[protein1][protein2]) max = weight[protein1][protein2];
//     // }
    
//     for(int i = 0; i < Original_ppi.interaction.size(); ++i){
//         int protein1 = Original_ppi.interaction[i].proteina;
//         int protein2 = Original_ppi.interaction[i].proteinb;
//         // Original_ppi.interaction[i].interaction = (Original_ppi.interaction[i].interaction + weight[protein1][protein2]) * 0.5;
//         double score = (Original_ppi.interaction[i].interaction + weight[protein1][protein2]); 
//         Original_ppi.interaction[i].interaction = score;
//         std::cout <<Id2protein[protein1] << "\t" << Id2protein[protein2] << ":" << score << endl;
    
//         // Original_ppi.interaction[i].interaction = (Original_ppi.interaction[i].interaction + weight[protein1][protein2] / max) * 0.5;
//     }
//     //
//     print_information(Ppidata_file, Result_file, Balanced_index);
//     get_balanced_interaction(Original_ppi,Balanced_index);
//     write_proteins(get_result(Original_ppi,Similarity_threshold),Result_file);


//     return 0;
// }

