# USING CRISPRVARIANT TO STYLE MULTIPLE ALIGNMENT

# Need Ref as a DNAString
  #datos$Ref already as DNAString
amp <- DNAString(datos$Ref$`TYR MOLDE`)
# PAM region & co
pam_region = 132
snip_loc = 136
guide <- IRanges::IRanges(pam_region,pam_region+23)
focus_interest_range <- IRanges(116,156)

crispr_data <- list(pam = pam_region, target_loc = snip_loc, guide_range = guide)

# All fastq data as a DNAStringSet
  # Use datos$Table to identify printing order
######### IF ACCEPTED METHOD #############
# Probably should use indel_final xlsx data
######### window of interest fix
# Another way to trim data is to trim it early and align later
#########
# Tabla <- read.xlsx(file = '', col = TRUE)
Tabla <- Fuse_entries(datos$Tabla) # Use fuse entries to ignore SNP and adapt Frequency

fasta_top10 <- inner_join(datos$fasta, 
                          slice(Tabla, 2:11) %>% select(ID, Freq, insert_range, deletion_range), by = 'ID') %>% 
                          arrange(desc((Freq)))

# indel_data <- c("-4ins","-3del","-1_-2del", "-3_-4del","-4_-5ins","-1_-3del",
                # "-3_-8del","-5_-6del","+6_-9del","-1_3del")
# allelic_id <- str_c('Indel: ',indel_data, ' Freq: ',fasta_top10$Freq)

############# BETA ##############
# Notation change: good for puntual changes, not valid for ranges.
# notat_ins <- ifelse(!is.na(fasta_top10$insert_range),paste0('+',as.numeric(unlist(lapply(str_split(fasta_top10$insert_range, pattern = ':'),'[[',1)))-pam_region-2,'Ins'),NA)
# notat_dels <- ifelse(!is.na(fasta_top10$deletion_range),paste0('+',as.numeric(unlist(lapply(str_split(fasta_top10$deletion_range, pattern = ':'),'[[',1)))-pam_region-2,'Del'),NA)
############# BETA ##############

reads <- Biostrings::DNAStringSet(fasta_top10$Sequence)
allelic_frequencies <- str_c('Freq', fasta_top10$Freq, sep = ': ')
names(reads) <- allelic_frequencies

# NeedleMan Wunsh pairwise alignment

pa <- Biostrings::pairwiseAlignment(DNAStringSet(replicate(length(reads),amp)),reads)
rc <- Biostrings::pairwiseAlignment(DNAStringSet(replicate(length(reads),amp)),reverseComplement(reads))

results <- getReferenceGaps(pa, reverse.alns = rc)

# Are we restrincting the visualization to a target focus area?
plotting_data <- plotAlignments_data(amp = amp, 
                                     results = results, 
                                     crispr_data = crispr_data, 
                                     focus_interest_range = focus_interest_range)
######
names(plotting_data$results$seqs) <- names(reads)
reads <- as.character(plotting_data$results$seqs)
# reads <- as.character(spam)
id_to_nm <- names(reads)
names(id_to_nm) <- seq_along(reads)
plotting_data$results$insertions$cigar <- id_to_nm[plotting_data$results$insertions$id]
plotting_data$results$insertions$count <- 1

xtick.breaks <- seq(from = 1, to = nchar(plotting_data$amp), by = 5)
xtick.labs <- xtick.breaks - plotting_data$crispr_data$target_loc
xtick.labs <- ifelse(xtick.labs <= 0, xtick.labs-1,xtick.labs)


############# PLOT data ################

p <- plotAlignments(plotting_data$amp, alns = reads, ins.sites = plotting_data$results$insertions,
                    pam.start = plotting_data$crispr_data$pam, pam.end = plotting_data$crispr_data$pam+2, 
                    target.loc = plotting_data$crispr_data$target_loc, 
                    guide.loc = plotting_data$crispr_data$guide_range,
                    legend.cols = 8, plot.text.size = 1.5, axis.text.size = 6,
                    legend.symbol.size = 3, legend.text.size = 5, line.weight = 0.5, 
                    highlight.guide = TRUE,
                    y.size = 50, tile.height = 0.40,
                    xtick.breaks = xtick.breaks, xtick.labs = xtick.labs)

p <- p + theme(plot.margin = grid::unit(c(0.3,0.3,0,0.3), "lines"),
               legend.margin = margin(0,unit = "lines"),
               axis.ticks = element_blank())

# PLOT only the reference

plotAlignments(amp, alns = NULL, target.loc = 215, ins.sites = data.frame(), renumbered = FALSE)


getReferenceGaps <- function(paired.alns, reverse.alns = NULL){
  elementLengths=elementNROWS
  browser()
  ########################
  # getReferenceGaps converts pairwise alignments into alignments with respect to a 
  # reference, and returns a list of data that can be input to CrispRVariants::plotVariants.matrix
  #
  # This function takes paired.alns, a Biostrings::PairwiseAlignments object of a set of reads 
  # to a reference sequence, and optionally reverse.alns, a matching PairwiseAlignments 
  # object of the reverse complement reads.
  ##########################
  
  refs <- as.character(pattern(paired.alns))
  qrys <- DNAStringSet(as.character(subject(paired.alns))) 
  
  # Use the reverse complement alignment when it is better
  if (! is.null(reverse.alns)){
    rc_better <- score(reverse.alns) > score(paired.alns)
    refs[rc_better] <- as.character(pattern(reverse.alns))[rc_better]
    tmp <- DNAStringSet(as.character(subject(reverse.alns)[rc_better]))
    qrys[rc_better] <- tmp
  }
  
  # Split into nucleotides
  ref_chrs <- strsplit(refs, "")
  
  # Get locations of gaps in reference sequence
  # make these into a list of insertions in the reads
  irl <- lapply(seq_along(ref_chrs), function(i){
    x <- ref_chrs[[i]]
    rls <- rle(x == "-")
    ir <- as(PartitioningByEnd(cumsum(rls$lengths)), "IRanges")
    gap_rngs <- ir[rls$values]
    gap_seqs <- Views(qrys[[i]], ir[rls$values])
    ungap_qry <- paste0(Views(qrys[[i]], ir[!rls$values]), collapse = "")     
    list(insertions = gap_rngs, ungapped = ungap_qry,
         gap_seqs = as.character(gap_seqs))
  })
  qrys <- DNAStringSet(unlist(lapply(irl, "[[", "ungapped")))
  ins <- lapply(irl, "[[", "insertions")
  #########
  #
  # Do gaps not match because of leading deletions not being printed?
  #
  #########
  # Offset the starts for the previous insertions
  temp <- IRangesList(ins)
  tt <- unlist(cumsum(width(temp))) # count inserted sequences
  tt <- c(0,tt[1:(length(tt)-1)]) # offset by one to consider only previous ins
  ee <- unlist(elementLengths(ins))
  ee <- cumsum(ee[ee>0])
  tt[ee[1:(length(ee)-1)] + 1] <- 0 # Zeros at the start of new reads
  starts <- unlist(start(temp)) - tt
  
  # Make a table of insertions
  insertions_t <- data.frame(start = starts, 
                             seq = unlist(lapply(irl, "[[", "gap_seqs")),
                             id = rep(seq_along(qrys), elementLengths(ins)))
  
  list(insertions = insertions_t, seqs = qrys)
}

#' @title Creates a text alignment from a set of cigar strings
#'@description Creates a one-to-one text alignment of a set of cigar strings with respect
#'to the reference sequence by collapsing insertions and introducing gaps
#'across deletions.
#'
#'When genomic coordinates for the alignment start and the target region
#'are provided, aligned sequences are cropped to the target region
#'@author Helen Lindsay
#'@param cigar A list of cigar strings to align
#'@param dnaseq The set of sequences corresponding to the cigars, as Biostrings::DNAStrings
#'@param target The target region to return, as GRanges.  Sequences overlapping
#'the target region are trimmed to exactly match it.
#'@param del_char The character to represent deleted bases. Default "-"
#'@param aln_start Genomic start locations of aligned sequences. Should be
#'used in conjunction with target_start and target_end.
#'@return The sequences with insertions collapsed and deletions padded
#'@rdname seqsToAln
seqsToAln <- function(cigar, dnaseq, target, del_char = "-", aln_start = NULL){
  # Additional trimming is necessary because deletion operations may overhang
  # boundaries of target region
  target_start <- start(target)
  target_end <- end(target)
  strand <- as.character(strand(target))
  if (strand == "*") strand <- "+"
  
  result <- GenomicAlignments::sequenceLayer(dnaseq, cigar, D.letter = del_char,
                                             N.letter = "N")
  sq_len <- GenomicAlignments::width(result)
  
  if (! is.null(aln_start) & ! is.null(target_start) & ! is.null(target_end) ){
    trim_start <- target_start - (aln_start - 1)
    
    if (any(trim_start < 0)){
      stop("dnaseq to be trimmed must start before the target location")
    }
    
    trim_end <- trim_start + (target_end - target_start)
    if (any(trim_end > sq_len)){
      stop("dnaseq is not long enough to trim to the target region")
    }
    result <- subseq(result, trim_start,trim_end)
  }
  if (strand == "-"){
    result <- reverseComplement(result)
  }
  result <- as.character(result)
  result
}

Fuse_entries <- function(x){# Should we fuse both indels ignoring SNP sheets?
  
  y <- x %>% select(-score) %>%
    group_by(Deletions, Deleted_bps, Insertions, Inserted_bps, insert_range, deletion_range) %>% 
    summarise(Abundance = sum(Abundance), Freq = sum(Freq),ID,SEQ) %>% 
    arrange(desc(Abundance)) %>% 
    distinct(Deletions, Deleted_bps, Insertions, Inserted_bps, insert_range, deletion_range, .keep_all = TRUE) %>% 
    relocate(Abundance, Freq, ID, .before = Deletions) %>% 
    relocate(ID, .before = Abundance) %>% 
    ungroup()
  return(y)
  
}

plotAlignments_data <- function(amp, results, crispr_data, focus_interest_range=IRanges(1,nchar(amp))){
  # Trims sequences and displaces insert data, pam and guide position
  
  if(width(focus_interest_range) == nchar(amp)){
    # No triming required
    return(list(amp = amp, 
                results = results, 
                crispr_data = crispr_data, 
                focus_interest_range = focus_interest_range))
    
  }else if(width(focus_interest_range) < nchar(amp)){
    # focus shorter than seq length: trimming is a go
    
    amp <- subseq(amp, start(focus_interest_range), end(focus_interest_range))
    
    results$seqs <- subseq(results$seqs, start(focus_interest_range), end(focus_interest_range))
    results$insertions$start <- results$insertions$start-start(focus_interest_range)
    
    crispr_data$pam <- crispr_data$pam - start(focus_interest_range)
    crispr_data$target_loc <- crispr_data$target_loc - start(focus_interest_range)
    
    start(crispr_data$guide_range) <- start(crispr_data$guide_range) - start(focus_interest_range)
    end(crispr_data$guide_range) <- end(crispr_data$guide_range) - start(focus_interest_range)
    
    return(list(amp = amp, results = results, crispr_data = crispr_data))
  }else if(width(focus_interest_range) > nchar(amp)){
    
    return(stop('target range larger than sequence length'))
    
  }
}
