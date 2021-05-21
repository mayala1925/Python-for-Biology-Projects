# Python-for-Biology-Projects
Projects that I finished for the class Practical Computation Biology for Biologists: Python

- Project 2 Instructions : We want to find out if adjacent genes change during ecdysone treatment. Gene expression changes upon addition of ecdysone is in the file "rpkm_ecdysone_to_ctl.txt", which has two columns separated by a tab. First column is Gene_name, second column named "Enrichment" is the log2 fold change of expression upon addition of ecdysone. The genes next to enhancers responsive to ecdysone are listed in "w_ecd_genes.list" and the control group is listed in "wo_ecd_genes.list". Make a boxplot of expression change comparing genes next to ecdysone enhancers and genes next to control enhancers.

--Output is the boxplot.

Use this to run project 2.

1. git clone https://github.com/mayala1925/Python-for-Biology-Projects.git
2. cd Python-for-Biology-Projects
3. python3 Project\ 2\ -\ Computational\ Python.py


Project 3 Instructions: We know that RNA-binding proteins like RBFOX2 often recognize their RNA targets through interaction with short RNA sequences (kmers).  Perhaps, then, RBFOX2 is recognizing a specific kmer that is enriched in the introns downstream of the exons it enhances.  Can we figure out what this kmer (or kmers) might be?

1) Go through each file and count the number of occurences of each 5-mer.  This might be best stored as a dictionary of the form {kmer : number of occurences}.
2) Calculate a log2 enrichment for each for each kmer in the affected relative to control by dividing the frequency of the kmer in the affected sequences to its frequency in control sequences and taking the log base 2 of the ratio.  Remember that the frequency of a kmer is that number of occurences of that kmer divided by the number of occurences of all kmers.  Frequencies must be between 0 and 1.  Enrichments will be positive if the kmer is enriched in the enhanced sequences and negative if enriched in the control sequences.
3) Calculate the statistical significance of that enrichment with a Fisher's exact test.  This can be done by importing the fisher_exact function from the scipy.stats module.  
4) Correct your Fisher's exact p-values for multiple hypothesis testing with a Bonferroni correction.  Very simply, multiply your original pvalue from fisher_exact by the number of kmers you are considering (there are 1024 possible 5mers).
5) Sort kmers by their corrected pvalues


Use this to run project 3. - BEWARE THIS WILL TAKE A FEW MINUTES TO RUN THE FISHER STATISTICS P-VALUES OF ENRICHMENT.
1. git clone https://github.com/mayala1925/Python-for-Biology-Projects.git
2. cd Python-for-Biology-Projects
3. python3 project3-fisher-stats.py 

-- Output is final dataframe sorted by kmer in ascending order and their corrected p values.
