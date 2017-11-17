# Internship at the Jacques Monod Institute - Paris

Scripts written and used during a 2 months bioinformatics internship

Data integration and data analysis on the charaterization of replication timing across 6 human cell lines

### Matrix creation

-Create matrix of replication timing and its conservation.
-Use raw data from the lab
-Define segments on the genome
-Clustering of Replication timing (RT) with a AHC approach 

### Multiline viewer

Can be used to plot the replication timing profile of the 6 cell lines across all genome (lots of PDF). 
Smoothing the curve with a Loess method

### Matrix Complete

Query data from databases or parsing of files downloaded on UCSC genome Browser and adding it to the main data. 
SNP study in another script 

### Data analysis and visualization

Script *boxplotFunction.R* and *encadrements.R*.
Also *getGeneList.R* can be used to compute conservation groups with Gene Ontology/Panther 

