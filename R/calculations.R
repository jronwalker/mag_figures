#' Calculate reads mapped per contig
#'
#' The functions takes a directory containing the samtools idxstats output files for each sample and 
#' outputs a matrix usable for downstream analysis.
#' @param coverm_out Text file of the output of CoverM. CoverM command should use the methods length, count, and covered_bases.
#' @param bin_dir Path to directory containing nucleotide fasta files of bins.
#' @param gtdbtk_dir Path to GTDB-tk output directory
#' #' @param flagstat_dir Path to directory with files containing the raw output of samtools flagstat for each sample. 
#' File names should end in _flagstat followed by .csv or .tsv or .txt. Everything before "_flagstat" will be used as the sample names. 
#' Do not include another "_" after "_flagstat".
#' @param min_reads Minimum number of reads mapped to a contig to be considered present in a metagenome. Default is 10 reads.
#' @param cov_cut Minimum breadth of coverage for a contig to be considered there. The idea here is that X% of the contig
#' needs to be covered by reads within a sample to be considered present. Default is to not perform coverage filtering.
#' @param prop TRUE or FALSE indicating if you want the output to be on a 1-100 percent scale. FALSE will output the 
#' fractional abundance as defined by FRAP. The default is TRUE.
#' @param per_ml Data frame containing the sample name in the 1st column and the number of bacteria per mL of each sample. 
#' Make sure to include all samples if you have replicate for the metagenomes but are using one count value then you will need
#' multiple rows for each replicate with the appropriate names.
#' @param rank_ids Determines whether or not taxonomic information will retain the prefix indicating the taxonomic level. 
#' Default is TRUE which retains the prefixes. FALSE will remove the prefixes.
#' @param tax_rank The taxonomic rank at which you would like the fractional abundance to be calculated. 
#' Ranks that can be used are "bin", "domain", "phylum", "class", "order", "family", "genus", "species". 
#' The input should be a character string. Default is "bin".
#' @return The output will be a data frame containing the fractional abundance of each taxonomic rank chosen.
#' @importFrom dplyr group_by
#' @importFrom dplyr group_by_if
#' @importFrom dplyr select
#' @importFrom dplyr select_if
#' @importFrom dplyr summarise
#' @importFrom dplyr summarise_if
#' @importFrom magrittr %>%
#' @export
frac_abund <- function(coverm_out, bin_dir, gtdbtk_dir, flagstat_dir,
                       min_reads = 10, cov_cut, prop = T, per_ml,
                       rank_ids = TRUE, tax_rank = "bin"){
  #Get contigs in bins information
  c2b <- contigs2bins(bin_dir)
  #Get reads mapped to contigs
  covm <- parse_coverm(coverm_out)
  #Apply the minimum number of reads cutoff
  covm[["counts"]][covm[["counts"]] <= min_reads] <- 0
  #If statement to handle to coverage cutoff
  if (missing(cov_cut)){
  } else #Else is when it is present
  {
    covm$counts[covm$coverage < cov_cut] <- 0
  }
  #Merge the counts and contig/bin information
  b_counts <- merge(c2b, covm[["counts"]], by = "Contig")
  #Calculate the counts per bin
  b_counts <- b_counts %>% select(-Contig) %>% 
    dplyr::group_by(bin) %>% 
    dplyr::summarise(across(everything(), sum))
  #Get total read counts per metagenome
  totals <- reads_meta(flagstat_dir)
  #If then statement to make sure the sample names are the same. Output an error otherwise
  if (all(sort(colnames(b_counts[,-1])) == sort(totals$samples)) == TRUE){
    #Make the reads mapped per read in metagenome calculation. 
    per_meta <- sweep(b_counts[,totals$samples], 2, totals$totals, "/")
    #Give the row names the bin names used for later on to make sure they are in the same order.
    rownames(per_meta) <- b_counts$bin
  } else 
    {
    stop("The samples with mapped reads and the samples with total reads do not have matching names. Double check that 
         your flagstat and idxfiles are all there and name similarily. ")
      }
  #Merge the lengths and contig/bin information
  b_lengths <- merge(c2b, covm[["length"]], by = "Contig")
  #Calculate the counts per bin
  b_lengths <- b_lengths %>% select(-Contig) %>% 
    dplyr::group_by(bin) %>% 
    dplyr::summarise(across(everything(), sum))
  #Make the length coefficient table
  b_lengths$len_coeff <- mean(b_lengths$length)/b_lengths$length
  #Multiple rows by the length coefficient
  frac_abund <- sweep(per_meta[b_lengths$bin,], 1, b_lengths$len_coeff, "*")
  #Move the row names to a new column
  frac_abund <- cbind(row.names(frac_abund), frac_abund)
  #Rename the new column
  colnames(frac_abund)[1] <- "bin"
  #Remove the rownames 
  rownames(frac_abund) <- seq(1, length(rownames(frac_abund)))
  #If statement to determine how to handle taxonomic ranks
  if (tax_rank == "bin"){
    #Get taxanomic information
    taxa <- read_gtdbtk(gtdbtk_dir, rank_ids = rank_ids)
    #Merge taxa information with fractional abundance
    frac_taxa <- merge(taxa, frac_abund, by = "bin") %>% as.data.frame()
  } else 
    {
      #Get taxanomic information
      taxa <- read_gtdbtk(gtdbtk_dir, rank_ids = rank_ids)
      #Determine the column of the taxonomic rank chosen
      tax_pos <- which(colnames(taxa) == tax_rank)
      #Subset taxa to just the ranks of interest
      taxa <- taxa[,c(1:tax_pos)] 
      #Merge taxa information with fractional abundance
      frac_taxa <- merge(taxa, frac_abund, by = "bin")
      #Add up fractional abundances across taxonomic groups
      frac_taxa <- frac_taxa %>% select(-bin) %>% #Remove the bin column
        dplyr::group_by_if(is.character) %>%  #group by the character columns
        dplyr::summarise_if(is.numeric, sum) %>% as.data.frame() #add up each of the numeric columns
    }
  #If statemnt to determine if proportions should be calculation or left at fractional abundances
  if (prop == T){
    prop_df <- frac_taxa %>% dplyr::select_if(is.numeric) %>% 
      as.matrix() %>%
      prop.table(2)
    prop_taxa <- frac_taxa %>% dplyr::select_if(is.character) %>%
      cbind(prop_df)
  } else 
    {
      prop_taxa <- frac_taxa
    }
  #If statement on how to handle per_ml calculations
  if (missing(per_ml)) { #If matrix is missing no calculation
    return(prop_taxa)
    } else if (all(dim(per_ml) == c((ncol(covm[["counts"]])-1), 2), is.data.frame(per_ml))) 
      { #Now check to make sure sample lengths are equal
      if (prop == T){ #What to do when its in proportions were calculated
        ml_vec <- per_ml[,2]
        names(ml_vec) <- per_ml[,1]
        ml_vec <- ml_vec[colnames(prop_df)]
        bin_ml <- sweep(prop_df, 2, ml_vec, "*")
        ml_taxa <- frac_taxa %>% dplyr::select_if(is.character) %>%
          cbind(bin_ml)
        return(ml_taxa)
        } else if (prop == F) #What to do if the proportions were not calculated
          {
          warning("You multiplied your fractional abundance by bacterial counts without converting to proportions!
            This will likely result in some really weird numbers that don't make a whole lot of sense but you make the rules.")
          ml_vec <- per_ml[,2]
          names(ml_vec) <- per_ml[,1]
          prop_df <- frac_taxa %>% dplyr::select_if(is.numeric)
          ml_vec <- ml_vec[colnames(prop_df)]
          bin_ml <- sweep(prop_df, 2, ml_vec, "*")
          ml_taxa <- frac_taxa %>% dplyr::select_if(is.character) %>%
            cbind(bin_ml)
          return(ml_taxa)
          }
      } else if (all(dim(per_ml) == c((ncol(covm[["counts"]])-1), 2), is.data.frame(per_ml)) == FALSE) #Report error if samples are not equal length
        {
        stop("Something is wrong with your input files. Either your per_ml dataframe is not 2 columns or your number of 
           samples is not consistent between your stats and per_mL")
        } 
}

#' Generate table ready for ggplot2.
#'
#' The functions takes a community data frame where columns are taxonomic information and abundances
#' and rows are unique taxa (MAGs, ASVs, OTUs) and returns a long formatted table for ggplot2. 
#' 
#' @param abund_df Data frame object containing taxonomic information and proportional/relative abundance information. 
#' Works best with the frac_abund output but any data frame object should work. Taxonomic column names should 
#' be one or multiple of the following: "bin", "domain", "phylum", "class", "order", "family", "genus", "species".
#' Sample column names can be anything but the columns should be numeric.
#' @param rank Taxonomic rank at which you would like to make the graph. If none chosen the lowest will be used.
#' @param min_cutoff The cutoff to aggregate a group into the "Others". The cutoff is compared with the maximum
#' relative abundance of that group in your each of your samples. If that group does not surpass the cutoff in
#' any sample it is added to "Others". Default is 0.01 (i.e. 1% in proportional data).
#' @return The output will be a data frame containing the fractional abundance of each taxonomic rank chosen.
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by_at
#' @importFrom dplyr summarise_if
#' @export
comm_long <- function(abund_df, rank, min_cutoff = 0.01) {
  trs <- c("domain", "phylum", "class", "order", "family", "genus", "species", "bin")
  if (missing(rank)) {
    #Return the matches of colnames and all ranks then find the highest value and return the string at that location
    rank <- trs[max(match(names(abund_df), trs, nomatch = 0))]
  } else {
    rank <- rank
  }
  #Find all of the ranks lower than the rank being used
  rm_rank <- trs[match(rank, trs)+1:length(trs)]
  #Find all the ranks including and above the rank used
  k_rank <- trs[1:match(rank, trs)]
  #Drop column that are in the rm_rank variable
  abund_df <- abund_df[,!colnames(abund_df) %in% rm_rank]
  #Group by the rank chosen and those above it and sum by groups
  abund_df <- abund_df %>% dplyr::group_by_at(k_rank) %>% 
    dplyr::summarise_if(is.numeric, sum) %>% 
    as.data.frame()
  #Determine the max in each row
  maxes <- abund_df %>% select_if(is.numeric) %>% apply(1,max)
  #Determine if the max meets the cutoff or not
  maxes <- maxes >= min_cutoff
  #Retain those that meet the cutoff
  keep <- abund_df[maxes, ]
  #Get all the groups that didn't make the cut
  left <- abund_df[!maxes, ]
  #Add up the unwanted relative abundances
  others <- left %>% select_if(is.numeric) %>% 
    apply(2, sum) %>% t() %>% as.data.frame()
  #Give them their taxonomic ranks and define as Other
  others[k_rank] <- "Others"
  #Re-order the columns to the original order
  others <- others[, colnames(keep)]
  #Bind the keepers with the others
  abund_df <- rbind(keep, others)
  #Make table long
  abund_long <- tidyr::pivot_longer(abund_df, names(select_if(abund_df, is.numeric)), 
                                    names_to = "Sample", values_to = "Rel. Abund")
  return(abund_long)
}
