
CMC
Version 2.0
=============

This package contains the following programs:

1/ cd_ppi:  A method to assign reliability to edges in a PPI network.
         The reliability is assigned based on the concept of
         interated CD distance. See [2].

2/ filterNadd_ppi: A method for cleansing a PPI network. It uses
         cd_ppi to assign reliability to edges in a PPI network.
         Those below a specified threshold are removed as false positives.
         It also evaluates protein pairs not in the PPI network and
         decides if they are possible false negatives. Those pairs
         scoring above a specified threshold are added into the network.
         See [1,2]
         
3/ CMC:  A method to identify protein complexes from a PPI network.
         The method is based on merging of maximal clique. See [1].


A detailed set of instructions for these programs can be found after
the references below.



Credits:

This package was implemented by LIU Guimei and is supported in part
by URC grant "R-252-000-274-112: Graph-Based Protein Function Prediction"
and by NRF grant "R252-000-218-281: Informatics and Search Algorithms
for Lipidomics - Novel Tools and Applications".

If you use this program, please cite the following references:

[1] Guimei Liu, Limsoon Wong, Hon Nian Chua. 
    Complex Discovery from Weighted PPI Networks.
    Bioinformatics, 2009 (to appear).

[2] Guimei Liu, Jinyan Li, Limsoon Wong. 
    Assessing and Predicting Protein Interactions Using Both Local and
    Global Network Topological Metrics. 
    Proceedings of 19th International Conference on Genome Informatics (GIW),
    pages 138--149, Gold Coast, Australia, 3 December 2008.

[3] Hon Nian Chua, Limsoon Wong. 
    Increasing the Reliability of Protein Interactomes. 
    Drug Discovery Today, 13(15/16):652--658, August 2008


===cd_ppi===

cd_ppi:  A method to assign reliability to edges in a PPI network.
         The reliability is assigned based on the concept of
         interated CD distance. See [2].

Usage of cd_ppi:

  cd_ppi  train_ppi_file  test_ppi_file  method  #iterations  output_filename


Parameters: 

1. train_ppi_file: contains protein protein interactions. 
	Each line represents an interaction, and contains a pair of proteins. 
	The program uses this file to calculate the score of every protein pair. 
2. test_ppi_file: contains the set of interactions to be assessed.
	Its format is the same as that of "train_ppi_file". The program will 
	calculate the score of the proteins pairs in this file, and output 
	those protein pairs with non-zero score. If "train_ppi_file" and 
	"test_ppi_file" are the same file, then the program is assessing 
	the reliability of the interactions in "train_ppi_file". If the 
	value of "test_ppi_file" is "NULL", then the program will predict 
	new protein interactions that are not in "train_ppi_file", and
        output those proteins pairs that are not in "train_ppi_file" and 
	their interacting scores. 

3. method: takes the following values: 
                 CD: CD-distance            (See [2,3])
            AdjstCd: Adjusted CD-distance   (See [2])
                 FS: FS-weight              (See [3])
    
4. nmax_iterations: is the number of iterations. 

5. output_filename: Each line contains a pair of proteins and their
	 interacting score. 



====filterNadd_ppi=====

filterNadd_ppi: A method for cleansing a PPI network. It uses
         cd_ppi to assign reliability to edges in a PPI network.
         Those below a specified threshold are removed as false positives.
         It also evaluates protein pairs not in the PPI network and
         decides if they are possible false negatives. Those pairs
         scoring above a specified threshold are added into the network.
         See [1,2]
         
Usage of filterNadd_ppi:

filterNadd_ppi ppi_filename \
	#iterations  filter_method  filter_min_score \
	add_method  add_min_score \
	output_file

This program calls the "cd_ppi" program to calculate PPI scores. 

Parameters:

1. ppi_filename: contains the set of interactions. 
	Each line represents an interaction, and contains a pair of proteins. 
	The program uses this file to calculate the score of 
	every protein pair. 

2. #iterations: number of iterations for the chosen scoring methods. 
	Suggested value: 2. 

3. filter_method: the scoring method used to remove non-reliable interactions. 
	It takes the following values: 
                 CD: CD-distance             (See [3])
            AdjstCd: Adjusted CD-distance    (See [1,2])
                 FS: FS-weight               (SEe [3])

4. filter_min_score: the threshold used to filter interactions with low 
	score from "ppi_filename".  If its value is between 0 and 1, 
	then all the interactions with score lower than its value are removed. 
	If its value is an integer k that is larger than 1, then only the 
	top-k interactions are retained. 

5. add_method: the scoring method used to add new interactions. 
	It takes the following values: 
                 CD: CD-distance
            AdjstCd: Adjusted CD-distance 
                 FS: FS-weight              

6. add_min_score: the threshold used to add new predicted interactions. 
	If its value falls in (0, 1], then those proteins pairs that are 
	not in "ppi_filename" but have a score no less than the given 
	value are added to the final output file. If its value is
        an integer k that is larger than 1, then only the top-k new 
	proteins pairs will be added. If do not want to add new interactions,
	set its value to 0. 

7. output_file: contains the final proteins pairs. Each line contains 
	a pair of proteins and their score. 
                  

Example: 

1. The following command uses AdjstCD to calculate score, and removes 
interactions with 0 score, and add the top-1000 new interactions to the
PPI network.  

	filterNadd_ppi  dip.ppi.txt 1 AdjstCD 0  AdjstCD 1000  ppi.score.txt  

2. Same as the above command except that it uses iterated AdjstCD to
calculate score. The number of iterations is 2. 

	filterNadd_ppi  dip.ppi.txt 2 AdjstCD 0 AdjstCD  1000  ppi.score.txt





====CMC=====

CMC:  A method to identify protein complexes from a PPI network.
         The method is based on merging of maximal clique. See [1].

Usage of CMC:

CMC ppi_score_filename \
    min_deg_ratio \
    min_size \
    overlap_thres \
    merge_thres \
    output_file

This program calls the "quasiCliques" program to find maximal cliques 

Parameters:

1. ppi_score_filename: contains the set of interactions and their scores.
	 Each line represents an interaction, 
         and contains a pair of proteins and their score.

2. min_deg_ratio: set it to 1

3. min_size: the minimum size of the clusters generated

4. overlap_thres: the threshold used to remove or merge highly
	 overlapped clusters. Given two clusters C1 and C2, if
         the overlap between C1 and C2 is no less than filter_score*|C2|,
	 then C2 will either be removed or merged.

5. merge_thres: the threshold used to remove or merge highly 
	 overlapped clusters. Given two clusters C1 and C2, if
         the overlap between C1 and C2 is no less than filter_score*|C2|,
         and the inter-connectivity between C1 and C2 is no less than 
         merge_thres, then C2 is merged with C1, otherwise C2 is removed.

6. output_file: contains the list of clusters generated. 
	Each line represents a cluster. The string 
        before ':' is the identifier of the cluster, followed by the
	set of proteins in this cluster. 



Example:

   CMC  ppi.score.txt  1 4 0.5 0.25  clusters.txt



 
