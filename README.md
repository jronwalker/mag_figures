# Helper functions to make generating figures in R a little easier.
See [Sulfur B_MAGs](https://github.com/Silveira-Lab/sulfur_bmags) for a more complete description of integrating this into a larger metagenomics pipeline.  

Currently, this programs main function is to take a CoverM output, a bin directory, a GTDB-tk output directory, and a directory containing outputs from samtools flagstat and convert this information into dataframes that are human readable and easily used in generating statistics about the community composition. It caalculates fractional abundances of bins by computing the fractional abundances of contigs within a bin and summing them. The fractional abundance is calculated by the equation:
$$f_i = (r_i / T_j) * {(L_m \over L_i)}$$
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
If you want to be stringent in whether or not you consider a contig to be present you can 
