# Helper functions to make generating figures in R a little easier.
See [Sulfur B_MAGs](https://github.com/Silveira-Lab/sulfur_bmags) for a more complete description of integrating this into a larger metagenomics pipeline.  

Currently, this programs main function is to take a CoverM output, a bin directory, a GTDB-tk output directory, and a directory containing outputs from samtools flagstat and convert this information into dataframes that are human readable and easily used in generating statistics about the community composition. It caalculates fractional abundances of bins by computing the fractional abundances of contigs within a bin and summing them. The fractional abundance is calculated by the equation:  
$$f_i = {r_i \over T_j} * {L_m \over L_i}$$  
where:  
- $f_i$ is the fractional abundance of contig i  
- $r_i$ is the number of reads mapped to contig i  
- $T_j$ is the total number of reads in metagenome j  
- $L_m$ is mean length of all contigs  
- $L_i$ is the length of contig i   
### How to generate the CoverM file
Currently, the program requires a CoverM output with at least the length, counts, and covered bases of each of your contigs. This can be doe by mapping your reads from each sample to your contigs and then running the following command on the BAM/SAM files.  
```
coverm contig -b *.bam --methods length count covered_bases > coverm_output.tsv
```
### Bin directory
The bin directory should be a directory containing 1 file per bin with all of the contigs for each bin in the file corresponding to the bin. This is a pretty standard output for most binning softwares. Ensure that these are the only files as have coding sequences, genes, or some other file will likely result in problems.
### GTDB-tk directory
This should be the output directory of the gtdbtk classify_wf.
### Flagstat directory
This should include text files containing the outputs of `samtools flagstat` for each of your SAM/BAM files. This was done using samtools 1.21 and I do believe older versions have a different output so make sure you are using a recent version of samtools. If you aren't sure if your samtools will work try a flagstat command and make sure it outputs a line ending in "primary".
### Other Options
#### Cutoffs
If you want to be stringent in whether or not you consider a contig to be present you can set two types of cutoffs with this program:  
- Number of Reads (min_reads) - Set an integer value for the minimum number of reads required for the contig to be counted
- Coverage (cov_cut) - Set a minimum percentage (scale 0-1) of the contig that is covered by reads mapped  
#### Proportions of metagenomes
You can set whether or not you want the fractional abundances of a metagenome to be scaled to 1 i.e. a percentage of all of the fractional abundances in said metagenome. The default behavior is to scale them you can turn it off with prop=F.  
#### Adding bacterial counts
If you calculated the proportion you may want to assess the actual counts of these proportions by multiplying the proportions by your bacterial counts. You can do this by providing a 2 column dataframe with the sample names (col1) and the bacterial counts (col2).
#### Controling the taxonomic levels and names
We can set the taxonomic level at which you would like all rows to be aggregated. The defaults is to output a MAG per row. But if we want the fractional abundances or proportions at the class level then we can set tax_rank="class" and each row will then output the summed abundances for all MAGs in each class.  
We can also turn on and off the id strings before each taxonomic level i.e. "d__" for domain, "p__" for phylum... by setting rank_ids=F.
