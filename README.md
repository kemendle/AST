# AST
A challenge in phylogenetic inference of gene trees is how to properly sample a large pool of homologous sequences to derive a good representative subset of the sequences. Such a need may arise for various application scenarios such as (1) accuracy-oriented phylogeny reconstruction methods may not be able to deal with a large pool of sequences due to their required computing resources; (2) some applications may prefer using small gene trees as in detection of horizontal gene transfers; and (3) the pool of sequences may be biased towards heavily studied species. To deal with such an issue, manual selection has been usually used as there is no automated sequence sampling method that can solve the resampling issue.  

Here we present a new computational method AST to sample representative sequences, which maximizes the taxonomic diversity of the sampled sequences. AST can be applied to (1) infer the evolutionary history of gene (or gene family) and (2) detect horizontal gene transfer as long as preparing suitable input files.  

The basic algorithm and its applications are described in the following publication:Chan Zhou, Fenglou Mao, Yanbin Yin, Jinling Huang, Johann Peter Gogarten and Ying Xu, AST: An Automated Sequence-Sampling Method for Improving the Taxonomic Diversity of Gene Phylogenetic Trees, Plos One, 2014 Jun 03,DOI: 10.1371/journal.pone.0098844 (http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0098844 ) 

Software requirements: AST was developed in PERL language, so PERL version 5.8.1 or later must be pre-installed. It has been installed, tested and run successfully in Linux OS with more than 7G memory.  It may be run in MS Windows with installed PERL.

User instruction is provided in the "README" file included in that package.
