dir.create("Individual_genes")
dir.create("Individual_genes_aligned")
dir.create("Individual_genes_aligned_trimmed")
dir.create("Individual_genes_gaps_filled")
dir.create("Individual_genes_gaps_filled/_Genes_need_more_check")
dir.create("Individual_genes_translated")

library(ape)
library(seqinr)
library(DECIPHER)
library(Biostrings)
library(stringr)

# obtain a list of gene names
The_reference <- read.fasta(file = "Data/Scyliorhinus_1266_gene_list_revised.txt", seqtype = "DNA", as.string = TRUE,
                               forceDNAtolower = FALSE,set.attributes = FALSE)
Gene_names_1266 <- names(The_reference)

# import the fasta file that contains all gene-captured sequences for multiple samples
Seq_list <- read.fasta(file = "Data/Poroderma_8_species_all_genes.txt",
           seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE)

Gene_to_list <- lapply(1:length(Gene_names_1266), function (x ) grep(Gene_names_1266[x], names(Seq_list)))


# to make a matrix showing which sample miss which gene
Sample_list <- unique(gsub("_.*","",names(Seq_list)))  # to get the list of sample names (GN numbers)
Gene_sample_matrix <- lapply(1:length(Gene_names_1266), function (x ) Sample_list %in% gsub("_.*","",names(Seq_list[Gene_to_list[[x]]])))

Gene_sample_matrix2 <- data.frame(lapply(as.data.frame(Gene_sample_matrix), function(x) {gsub("TRUE","Captured", gsub("FALSE","_", x))}))
colnames(Gene_sample_matrix2) <- Gene_names_1266
rownames(Gene_sample_matrix2) <- Sample_list                         

write.csv(t(Gene_sample_matrix2), "The_gene_sample_matrix.csv", row.names = TRUE)  # must use write.csv; table is transposed


# make a list to show the genes not captured in any sample (only in the reference)
Gene_in_one_sample <- list()
for (i in 1:length(Gene_names_1266))
  {
  if (length(Gene_to_list[[i]])==1)
    {
    Gene_in_one_sample<- append(Gene_in_one_sample,Gene_names_1266[i])
    }
  }
write.table(Gene_in_one_sample, "List_of_genes_not_captured_in_any_sample.csv",
            append = FALSE, quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)

# write in the folder "Individual_genes" DNA sequences sorted by gene
lapply(1:length(Gene_names_1266), function (x) write.fasta(sequences = Seq_list[Gene_to_list[[x]]],
                names = names(Seq_list[Gene_to_list[[x]]]), 
                file.out = paste0("Individual_genes/",names(Seq_list)[x],".fasta")))


# write in the folder "Individual_genes_aligned" the aligned DNA sequences sorted by gene
Indiv_genes<- list.files(path="Individual_genes/.", pattern=".fasta")

for (i in 1:length(Gene_names_1266))
  {
  seqs <- readDNAStringSet(paste0("Individual_genes/", Indiv_genes[i]))
  names(seqs) <- gsub("_.*","",names(seqs))                              # keep only GN numbers in the sample name
  if (length(seqs)<3) {next}                                            # skip the gene that can only be found in a small number of samples. Important parameter !!!!!!!
  seqs_gapless <- RemoveGaps(seqs, removeGaps = "all")
  writeXStringSet(AlignSeqs(seqs_gapless), file=paste0("Individual_genes_aligned/Aligned_", Indiv_genes[i]))
}


# trim sample sequences to the same length as the reference sequence and write results to the folder "Individual_genes_aligned_trimmed"
Indiv_align<- list.files(path="Individual_genes_aligned/.", pattern=".fasta")

Percent_ambiguity <- list()
Seq_name <- list()
for (i in 1:length(Indiv_align))
  {
  aligned_seq <- readDNAStringSet(paste0("Individual_genes_aligned/", Indiv_align[i]))

  first_nt <- min(start(matchPattern('a', aligned_seq[[1]])[1]),start(matchPattern('t', aligned_seq[[1]])[1]),
                  start(matchPattern('c', aligned_seq[[1]])[1]),start(matchPattern('g', aligned_seq[[1]])[1]))
  
  last_nt <- max(start(matchPattern('a', aligned_seq[[1]])[-1]),start(matchPattern('t', aligned_seq[[1]])[-1]),
                start(matchPattern('c', aligned_seq[[1]])[-1]),start(matchPattern('g', aligned_seq[[1]])[-1]))
  
  trimmed_seq <- substring(as.character(aligned_seq), first_nt-0, last_nt+0) # change the first "0" and the second "0" to keep more nucleotides. 
                                                                             # If trim too much, you will get an "out-of-bounds" error.
  
  writeXStringSet(DNAStringSet(trimmed_seq), file=paste0("Individual_genes_aligned_trimmed/Trimmed_", Indiv_align[i]))
  
  seq_chr <- as.character(trimmed_seq)
  Number_ambiguity <- sum(str_count(seq_chr, "-")+str_count(seq_chr, "N")+str_count(seq_chr, "R")+str_count(seq_chr, "Y")
                          +str_count(seq_chr, "Y")+str_count(seq_chr, "M")+str_count(seq_chr, "K")+str_count(seq_chr, "W")
                          +str_count(seq_chr, "S")+str_count(seq_chr, "H")+str_count(seq_chr, "V")+str_count(seq_chr, "D")
                          +str_count(seq_chr, "B"))
  
  Char_in_gene <- sum(nchar(seq_chr))

  Percent_ambiguity <- append(Percent_ambiguity,round(Number_ambiguity*100/Char_in_gene,digits = 1))
  Gene_name_short <- gsub("[:.:].*","", Indiv_align[i])
  Seq_name<- append(Seq_name, Gene_name_short)
}

ambiguity_list <- as.data.frame(cbind(Seq_name,Percent_ambiguity))
List_ambiguity_sorted <- ambiguity_list[order(as.numeric(Percent_ambiguity),decreasing = TRUE),]  # must add "as.numeric"
List_ambiguity_sorted_chr <- apply(List_ambiguity_sorted,2,as.character)  # make dataframe two lists of characters

write.table(List_ambiguity_sorted_chr, "Percentage_of_ambiguous_nucleotides_found_in_each_gene.csv",append = FALSE, quote = FALSE, sep = ",", 
            row.names = FALSE, col.names = TRUE)


# translate each individual gene alignment into aa sequences
Gene_align<- list.files(path="Individual_genes_aligned_trimmed/.", pattern=".fasta")

Number_stop <- list()
Seq_name_list <- list()
for (i in 1:length(Gene_align))
  {
  Single_gene_align <- readDNAStringSet(file=paste0("Individual_genes_aligned_trimmed/", Gene_align[i]),format="fasta")
  
  First_in_align_noGap1 <- RemoveGaps(Single_gene_align[1])   # remove gaps; lost names
  First_in_align_noGap2 <- DNAStringSet(paste0("g",First_in_align_noGap1))   # shift the reading frame one nucleotide
  First_in_align_noGap3 <- DNAStringSet(paste0("gg",First_in_align_noGap1))   # shift the reading frame two nucleotides
 
  aa_seq1 <- translate(First_in_align_noGap1, genetic.code=GENETIC_CODE, no.init.codon=FALSE,if.fuzzy.codon="error")
  aa_seq2 <- translate(First_in_align_noGap2, genetic.code=GENETIC_CODE, no.init.codon=FALSE,if.fuzzy.codon="error")
  aa_seq3 <- translate(First_in_align_noGap3, genetic.code=GENETIC_CODE, no.init.codon=FALSE,if.fuzzy.codon="error") 
  
  ORF1_stop <- str_count(str_sub(as.character(aa_seq1), 1,), "\\*") 
  ORF2_stop <- str_count(str_sub(as.character(aa_seq2), 1,), "\\*") 
  ORF3_stop <- str_count(str_sub(as.character(aa_seq3), 1,), "\\*")
  
  Best_ORF <- which.min(c(ORF1_stop,ORF2_stop,ORF3_stop)) # find the ORF that gives the fewest stop codons
  N_to_add <- DNAStringSet(str_sub("NN",0,Best_ORF-1))  # add 0-2 "N" to the beginning of every sequence to shift the ORF
  
  #print(paste0("The best ORF for ",Gene_align[i], " is ", Best_ORF))
  
  Gene_noGap0 <- DNAStringSet(gsub("-", "N", Single_gene_align))           # replace gaps with N
  
  Gene_noGap1 <- xscat(N_to_add,Gene_noGap0)           # shift the reading frame 0-2 nucleotides; lost names
  #Gene_noGap1 <- DNAStringSet(paste0(N_to_add,Gene_noGap0))           # also works; lost names
  
  # make all ORF starts from 1
  if (Best_ORF==1) {
    Gene_noGap2 <- Gene_noGap1
    }
  
  if (Best_ORF==2 | Best_ORF==3 ) {
    Gene_noGap2 <- subseq(Gene_noGap1, start = 4)
  }
  
  # make the length of each alignment dividable by 3
  if (width(Gene_noGap2[1])%%3==0) {
    Gene_noGap2 <- Gene_noGap2
  }
  
  if (width(Gene_noGap2[1])%%3==1) {
    Gene_noGap2 <- subseq(Gene_noGap2, end=width(Gene_noGap2[1])-1)
  }
  
  if (width(Gene_noGap2[1])%%3==2) {
    Gene_noGap2 <- subseq(Gene_noGap2, end=width(Gene_noGap2[1])-2)
  }
  
  names(Gene_noGap2) <- names(Gene_noGap0)           # assign names again
  
  writeXStringSet(Gene_noGap2, file=paste0("Individual_genes_gaps_filled/Gapless_", Gene_align[i]))
  
  aa_sequence <- translate(Gene_noGap2, genetic.code=GENETIC_CODE, no.init.codon=FALSE,if.fuzzy.codon="solve")
  
  writeXStringSet(aa_sequence,file=paste0("Individual_genes_translated/","AA_", Gene_align[i]))
  
  Number_stop_gene <- sum(str_count(str_sub(as.character(aa_sequence), 1,), "\\*"))  # count the total number of stop codons in an alignment
  Number_stop <- append(Number_stop,Number_stop_gene)
  
  #Gene_name_short <- gsub("Trimmed_Aligned_List_of_", "", (gsub("[:.:].*","", Gene_align[i])))
  Gene_name_gapless <- paste0("Gapless_",Gene_align[i])
  Seq_name_list<- append(Seq_name_list, Gene_name_gapless)
}

stop_list <- as.data.frame(cbind(Seq_name_list,Number_stop))
List_stop_sorted <- stop_list[order(as.numeric(Number_stop),decreasing = TRUE),]  # must add "as.numeric"
List_stop_sorted_chr <- apply(List_stop_sorted,2,as.character)  # make a dataframe with two lists of characters

write.table(List_stop_sorted_chr, "List_of_stop_codons_found_in_each_gene.csv",append = FALSE, quote = FALSE, sep = ",", 
            row.names = FALSE, col.names = TRUE)

# Move genes with stop codon detected to the folder "_Genes_need_more_check"
Gene_recheck <- stop_list$Seq_name_list[which(stop_list$Number_stop>0)]

for (i in 1:length(Gene_recheck)){
  file.rename(from =paste0("Individual_genes_gaps_filled/",Gene_recheck[i]),
              to=paste0("Individual_genes_gaps_filled/_Genes_need_more_check/",Gene_recheck[i]))
}

# Move genes with too many ambiguous sites (including "-" and "N") to the folder "_Genes_need_more_check"
High_ambiguity <- ambiguity_list$Seq_name[which(ambiguity_list$Percent_ambiguity>=10)]
High_ambiguity_unique <- setdiff(as.list(paste0("Gapless_Trimmed_",High_ambiguity,".fasta")),Gene_recheck) # exclude those already moved to the folder due to stop codons

for (i in 1:length(High_ambiguity_unique)){
  file.rename(from =paste0("Individual_genes_gaps_filled/",High_ambiguity_unique[i]),
              to=paste0("Individual_genes_gaps_filled/_Genes_need_more_check/",High_ambiguity_unique[i]))
}




######################################################################################################################################
# The following codes are for concatenating genes

#dir.create("Genes_to_concatenate")

Genes_to_concat <- list.files(path="Genes_to_concatenate/.", pattern=".fasta")      # place genes you want to concatenate in the folder

#sample_to_concat <- c("Amblyraja","GN22813","GN22829","GN22833","GN22834")              # this is the final list of samples
sample_to_concat <- Sample_list

# some sample missing some genes. Need to add those genes in and their sequences contain only "N"
for (i in 1:length(Genes_to_concat))
  {
  Gene_seq <- read.fasta(file = paste0("Genes_to_concatenate/",Genes_to_concat[i]), seqtype = "DNA", 
           as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE)
  Missed_sample <- setdiff(sample_to_concat, names(Gene_seq))
  for (j in 1:length(Missed_sample))
    {
    new_seq <- paste0(rep("N",nchar(Gene_seq[[1]][1])),collapse = "")
    names(new_seq) <- Missed_sample[j]
    Gene_seq <- append(Gene_seq, new_seq)
  }
  write.fasta(sequences = Gene_seq, names = names(Gene_seq), 
              file.out = paste0("Genes_to_concatenate/__Revised_",Genes_to_concat[i]))
}


# To double check if each of the gene to be concatenated is dividable by 3
Check_ORF <- list()
for (i in 1:length(Genes_to_concat))
{
  Gene_seq <- read.fasta(file = paste0("Genes_to_concatenate/",Genes_to_concat[i]), seqtype = "DNA", 
                         as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE)
  if (nchar(Gene_seq$Reference)%%3!=0)
    {
     Check_ORF <- append(Check_ORF, paste0(nchar(Gene_seq$Reference)/3, "_", Genes_to_concat[i]))
    }
}
write.table(Check_ORF,file = "Check_ORF_results.txt", append = FALSE, quote = FALSE, row.names = FALSE, col.names = FALSE)



# put all revised sequences into one file and then concatenate genes by sample names
Gene_concat_all <- lapply(1:length(Genes_to_concat), function (x) read.fasta(file = paste0("Genes_to_concatenate/__Revised_",Genes_to_concat[x]), 
                            seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE))

All_gene_concat <- list()
for (i in 1:length(sample_to_concat))
  {
  One_gene_concat <- c()
  for (j in 1:length(Genes_to_concat))
    {
    sample_ind <- which(names(Gene_concat_all[[j]])==sample_to_concat[i])
    One_gene_concat <- paste0(One_gene_concat, as.character(Gene_concat_all[[j]][sample_ind]))
  }
  All_gene_concat <- append(All_gene_concat, One_gene_concat)
}

# Replace GN numbers with real species names
Species_fullname <- read.table("Data/List_of_species_names.csv")
#Species_fullnameGN <- lapply(1:nrow(Species_fullname), function (x) gsub(".*_GN","GN",Species_fullname$V1[x]))     # if GN number in the back
Species_fullnameGN <- lapply(1:nrow(Species_fullname), function (x) gsub("_.*","",Species_fullname$V1[x]))      # if GN number in the front

matched_full <- lapply(2:(nrow(Species_fullname)+1), function (x) Species_fullname$V1[which(Species_fullnameGN==sample_to_concat[x])])

write.fasta(sequences = All_gene_concat, names = c("Reference", unlist(matched_full)), file.out = paste0("Concatenated_gene_sequences.txt"))

unlink("Genes_to_concatenate/__Revised_*")







