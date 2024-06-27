# library(dplyr)
# library(GenomicRanges)
# library(genomation)
# library(tidyr)

#' CIRCUS (CIrcular RNA Coordinate Update System)
#' @import dplyr
#' @import GenomicRanges
#' @import genomation
#' @import tidyr
#' @importFrom methods as
#' @importFrom utils read.csv read.table
#' @importFrom S4Vectors subjectHits DataFrame

#' @name CIRCUS
#' @description Converts BED12 or BED16 coordinates of circular RNAs into Universal Identifiers (UID)
#' according to nomenclature proposed by Chen et al. (PMID: 36658223).
#'
#' @param ref_gpf_path A txt file of two columns with no headers: first column with gene names,
#' second column with transcript IDs. It should match all transcripts in ref_path.
#' See 'transcript2gene.txt' in test data.
#' @param ref_path BED 12 file with no headers of reference transcripts for annotation. Must
#' contain all 12 columns indicating exon composition of each transcript. 4th column should
#' contain transcript IDs. see 'annotation.bed' in test data.
#' @param bed_path BED file with circRNAs to convert into circRNA UID. 4th column should
#' contain an original identifier of choice (e.g., circRNA_1). Indicate whether this file is
#' BED12 or BED6 with the bed6 parameter.
#' @param bed6 Boolean (T/F) indicator of whether circRNA BED file is BED12 or BED6. Note:
#' When BED6 input is used, output circRNA UIDs will only contain information about the 5'
#' and 3' terminal exon blocks (e.g., circZFAND6(3,5)) without intermediate exons (e.g., a
#' BED12 input would result in circZFAND6(3,4L,5)).
#'
#' @return Data frame with circRNA UID and annotation information of input circRNAs.
#' "isoform_index" indicates alphabetized indices in the case of duplicated UIDs of circRNA
#' isoforms that share same exon/intron anatomy but with different splice sites.
#' @examples
#' library(CIRCUS)
#' test_bed12 <- CIRCUS(ref_gpf_path = system.file("extdata",'transcript2gene.txt',
#' package="CIRCUS"),ref_path = system.file("extdata",'annotation.bed',package="CIRCUS"),
#' bed_path = system.file("extdata",'bed12_coords.bed',package="CIRCUS"),bed6=FALSE)
#' test_bed6 <- CIRCUS(ref_gpf_path = system.file("extdata",'transcript2gene.txt',
#' package="CIRCUS"),ref_path = system.file("extdata",'annotation.bed',package="CIRCUS"),
#' bed_path = system.file("extdata",'bed6_coords.bed',package="CIRCUS"),bed6=TRUE)
#' @export

CIRCUS <- function(ref_gpf_path,ref_path,bed_path,bed6=FALSE){

  #######List of gene to transcript names as input No. 1#############
  ref_gpf <- read.csv(ref_gpf_path, header = FALSE, stringsAsFactors = FALSE, sep = '\t')
  names(ref_gpf) <- c('gene','name')

  #######BED12 annotation as input No. 2#############
  ref_bed <- read.csv(ref_path, header = FALSE, stringsAsFactors = FALSE, sep = '\t')
  names(ref_bed) <- c('chrom', 'start', 'end','name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts')
  ref_bed <- ref_bed %>% select(chrom, start, end, name, strand, blockCount, blockSizes, blockStarts)

  if (bed6 == F){

    #######Case 1: BED12 circRNA coordinates as input No. 3#############
    message('BED12 circRNA coordinates as input')
    ##Load BED12 file and keep relevant columns for annotation##
    bed12 <- read.table(bed_path, header = F, sep = "\t")
    names(bed12) <- c('chrom', 'start', 'end','name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts')
    bed12$n_row <- seq_len(nrow(bed12))
    bed<-bed12[,c(1:3,6,10:12)]
    names(bed) <- c('chrom','start','end','strand','blockCount','blockSizes','blockStarts')

    ##Genomation annotation##
    gene.obj = readTranscriptFeatures(ref_path,
                                      up.flank=1000,
                                      down.flank=1000,
                                      remove.unusual = TRUE)
    sig.ann.obj = annotateWithGeneParts(as(bed, "GRanges"), gene.obj,strand = TRUE)
    tss = sig.ann.obj@dist.to.TSS
    bed$n_row <- seq_len(nrow(bed))
    m = as(merge(bed,tss,by.x = c("n_row","strand"), by.y = c("target.row","feature.strand"),all.x=T),"GRanges")
    m_df <- data.frame(m)
    m_df$seqnames <- as.character(m_df$seqnames)

    #Identify the indices of the exons in the annotation object that overlap the circRNAs
    hits = findOverlaps(m,gene.obj$exons,select="all")
    #Extract the indices
    idx <- unique(subjectHits(hits))
    #Get the feature names and exon numbers for the annotations that match the indices
    values <- DataFrame(gene.obj$exons[idx])
    values_df <- data.frame(values)

    bed_df <- merge(m_df,values_df,by.x='feature.name',by.y='name',all.x=TRUE)

    ###For + strand###
    message('Processing + strand')
    bed_df_positive <- bed_df %>%
      filter(strand == "+")%>%
      arrange(n_row,score)

    #Identify novel exons that do not match any annotated exons
    bed_df_positive$blockCount <- as.numeric(bed_df_positive$blockCount)
    bed_df_positive$exon_name <- ifelse(
      is.na(bed_df_positive$X.seqnames) & !is.na(bed_df_positive$feature.name),
      "NE",
      NA
    )

    get_exons_sites_positive <- function(bed, skip=0) {
      exon_start_size <- data.frame(
        start <- as.integer(unlist(strsplit(as.character(bed$blockStarts), ","))),
        size <- as.integer(unlist(strsplit(as.character(bed$blockSizes), ",")))
      )
      df <- bed[rep(1:nrow(bed), bed$blockCount), ] # replicate rows occurs as many as its exon number
      df$end <- df$start + exon_start_size$start + exon_start_size$size
      df$start <- df$start + exon_start_size$start
      return(df)
    }

    bed_exon_positive <- get_exons_sites_positive(bed_df_positive)
    bed_exon_positive <- bed_exon_positive %>%
      group_by(n_row) %>%
      mutate(sub_row = paste0(n_row,'-',dense_rank(paste(seqnames, start, end, strand)))) %>%
      ungroup() %>%
      arrange(n_row,sub_row)

    #Discuss all cases of annotated exons
    bed_exon_positive <- bed_exon_positive %>%
      mutate(
        exon_name =
          case_when(
            start > X.start & (end == X.end + 1 | end == X.end - 1 | end == X.end) ~ paste0('S', score),
            (start == X.start - 1 | start == X.start + 1 | start == X.start) & end < X.end ~ paste0(score, 'S'),
            start > X.start & end < X.end ~ paste0('S', score, 'S'),
            start < X.start & (end == X.end + 1 | end == X.end - 1 | end == X.end) &
              (sub_row != lag(sub_row) | (sub_row == lag(sub_row) & start > lag(X.end)) | is.na(lag(sub_row))) ~ paste0('L', score),
            start < X.start & (end == X.end + 1 | end == X.end - 1 | end == X.end) &
              sub_row == lag(sub_row) & (start == lag(X.end) | start == lag(X.end)+1 | start == lag(X.end)-1) ~ paste0('RI',',', score),
            (start == X.start - 1 | start == X.start + 1 | start == X.start) & end > X.end &
              (sub_row != lead(sub_row) | (sub_row == lead(sub_row) & end < lead(X.start)) | is.na(lead(sub_row))) ~ paste0(score, 'L'),
            (start == X.start - 1 | start == X.start + 1 | start == X.start) & end > X.end &
              sub_row == lead(sub_row) & (end == lead(X.start) | end == lead(X.start) +1 | end == lead(X.start) -1) ~ paste0(score, ',', 'RI'),
            start < X.start & end > X.end & (sub_row != lead(sub_row) | (sub_row == lead(sub_row) & end < lead(X.start)) | is.na(lead(sub_row))) &
              (sub_row != lag(sub_row) | (sub_row == lag(sub_row) & start > lag(X.end)) | is.na(lag(sub_row))) ~ paste0('L', score, 'L'),
            start < X.start & end > X.end & sub_row == lag(sub_row) & (start == lag(X.end)| start == lag(X.end) + 1 | start == lag(X.end) -1 ) &
              (sub_row != lead(sub_row) | (sub_row == lead(sub_row) & end < lead(X.start)) | is.na(lead(sub_row)))  ~ paste0('RI', ',', score, 'L'),
            start < X.start & end > X.end & sub_row == lead(sub_row) & end < lead(X.start) &
              (sub_row != lag(sub_row) | (sub_row == lag(sub_row) & start > lag(X.end)) | is.na(lag(sub_row))) ~ paste0('L', score, ',', 'RI'),
            start < X.start & end > X.end & sub_row == lead(sub_row) & end < lead(X.start) & sub_row == lag(sub_row) &
              (start == lag(X.end)| start == lag(X.end) + 1 | start == lag(X.end) -1 ) & (end == lead(X.start) | end == lead(X.start) +1 | end == lead(X.start) -1) ~ paste0('RI', ',', score, ',', 'RI'),
            start < X.start & end > X.start & end < X.end & (sub_row != lag(sub_row) | (sub_row == lag(sub_row) & start > lag(X.end)) | is.na(lag(sub_row))) ~ paste0('L', score, 'S'),
            start < X.start & end > X.start & end < X.end & sub_row == lag(sub_row) & (start == lag(X.end)| start == lag(X.end) + 1 | start == lag(X.end) -1 ) ~ paste0('RI', ',', score, 'S'),
            start > X.start & start < X.end & end > X.end & (sub_row != lead(sub_row) | (sub_row == lead(sub_row) & end < lead(X.start)) | is.na(lead(sub_row))) ~ paste0('S', score, 'L'),
            start > X.start & start < X.end & end > X.end & sub_row == lead(sub_row) & (end == lead(X.start) | end == lead(X.start) +1 | end == lead(X.start) -1) ~ paste0('S', score, ',','RI'),
            start >= X.end & (sub_row != lead(sub_row) | (sub_row == lead(sub_row) & end <= lead(X.start)) | is.na(lead(sub_row))) ~ 'NE',
            end <= X.start  & (sub_row != lag(sub_row) | (sub_row == lag(sub_row) & start >= lag(X.end) & lag(exon_name) != 'NE') | is.na(lag(sub_row))) ~ 'NE',
            TRUE ~ exon_name),
        exon_name
      )

    #Identify ciRNA cases
    bed_exon_positive <- bed_exon_positive %>%
      mutate(
        exon_name =
          case_when(
            (start == X.end | start == X.end + 1 | start == X.end - 1) & n_row == lead(n_row) &
              (end == lead(X.start) | end == lead(X.start) + 1 | end == lead(X.start) - 1) ~ paste0('I_',score),
            TRUE ~ exon_name),
        exon_name
      )

    #Label all the rest of the exons
    bed_exon_positive <- bed_exon_positive %>%
      mutate(
        exon_name =
          case_when (
            (start == X.start - 1 | start == X.start + 1 | start == X.start) &
              (end == X.end + 1 | end == X.end - 1 | end == X.end) ~ as.character(score),
            TRUE ~ exon_name),
        exon_name
      )

    #Merge transcripts with their gene names
    bed_exon_positive <- merge(bed_exon_positive,ref_gpf,by.x='feature.name',by.y='name',all.x=TRUE)

    #Rank all exons by exon index
    bed_exon_positive <- bed_exon_positive[order(bed_exon_positive$X.score),]

    #Generate circRNA UID
    df_positive <- bed_exon_positive %>%
      filter(!is.na(exon_name)) %>%
      group_by(n_row) %>%
      mutate(new_exon_name = paste(exon_name, collapse = ",")) %>%
      ungroup() %>%
      distinct(n_row, .keep_all = TRUE)
    df_positive1 <- df_positive %>% select(feature.name,gene,n_row,new_exon_name)
    df_positive1$circRNA_name <- NA
    df_positive1 <- df_positive1 %>%
      mutate(circRNA_name =
               case_when(
                 !grepl('I_',new_exon_name)  ~ paste0('circ',gene,'(',new_exon_name,')'),
                 TRUE ~ gsub('I_','',paste0('ci',gene,'(',new_exon_name,')'))))

    ###For - strand###
    message('Processing - strand')
    bed_df_negative <- bed_df %>%
      filter(strand == "-")%>%
      arrange(n_row,score)

    bed_df_negative$blockCount <- as.numeric(bed_df_negative$blockCount)

    #Identify novel exons that do not match any annotated exons
    bed_df_negative$exon_name <- ifelse(
      is.na(bed_df_negative$X.seqnames) & !is.na(bed_df_negative$feature.name),
      "NE",
      NA
    )

    get_exons_sites_negative <- function(bed, skip=0) {
      exon_start_size <- data.frame(
        start <- as.integer(unlist(strsplit(as.character(bed$blockStarts), ","))),
        size <- as.integer(unlist(strsplit(as.character(bed$blockSizes), ",")))
      )
      df <- bed[rep(1:nrow(bed), bed$blockCount), ] # replicate rows occurs as many as its exon number
      df$start <- df$end - exon_start_size$start - exon_start_size$size
      df$end <- df$end - exon_start_size$start
      return(df)
    }

    bed_exon_negative <- get_exons_sites_negative(bed_df_negative)
    bed_exon_negative <- bed_exon_negative %>%
      mutate(start = as.numeric(start)) %>%
      group_by(n_row) %>%
      arrange(desc(start), .by_group = TRUE) %>%
      mutate(sub_row = paste0(n_row, '-', dense_rank(start))) %>%
      ungroup()%>%
      arrange(n_row,sub_row,X.start)

    #Discuss all cases of annotated exons
    bed_exon_negative <- bed_exon_negative %>%
      mutate(
        exon_name =
          case_when(
            start > X.start & (end == X.end + 1 | end == X.end - 1 | end == X.end) ~ paste0('S', score),
            (start == X.start - 1 | start == X.start + 1 | start == X.start) & end < X.end ~ paste0(score, 'S'),
            start > X.start & end < X.end ~ paste0('S', score, 'S'),
            start < X.start & (end == X.end + 1 | end == X.end - 1 | end == X.end) &
              (sub_row != lag(sub_row) | (sub_row == lag(sub_row) & start > lag(X.end)) | is.na(lag(sub_row))) ~ paste0('L', score),
            start < X.start & (end == X.end + 1 | end == X.end - 1 | end == X.end) &
              sub_row == lag(sub_row) & (start == lag(X.end) | start == lag(X.end)+1 | start == lag(X.end)-1) ~ paste0('RI',',', score),
            (start == X.start - 1 | start == X.start + 1 | start == X.start) & end > X.end &
              (sub_row != lead(sub_row) | (sub_row == lead(sub_row) & end < lead(X.start)) | is.na(lead(sub_row))) ~ paste0(score, 'L'),
            (start == X.start - 1 | start == X.start + 1 | start == X.start) & end > X.end &
              sub_row == lead(sub_row) & (end == lead(X.start) | end == lead(X.start) +1 | end == lead(X.start) -1) ~ paste0(score, ',', 'RI'),
            start < X.start & end > X.end & (sub_row != lead(sub_row) | (sub_row == lead(sub_row) & end < lead(X.start)) | is.na(lead(sub_row))) &
              (sub_row != lag(sub_row) | (sub_row == lag(sub_row) & start > lag(X.end)) | is.na(lag(sub_row))) ~ paste0('L', score, 'L'),
            start < X.start & end > X.end & sub_row == lag(sub_row) & (start == lag(X.end)| start == lag(X.end) + 1 | start == lag(X.end) -1 ) &
              (sub_row != lead(sub_row) | (sub_row == lead(sub_row) & end < lead(X.start)) | is.na(lead(sub_row)))  ~ paste0('RI', ',', score, 'L'),
            start < X.start & end > X.end & sub_row == lead(sub_row) & end < lead(X.start) &
              (sub_row != lag(sub_row) | (sub_row == lag(sub_row) & start > lag(X.end)) | is.na(lag(sub_row))) ~ paste0('L', score, ',', 'RI'),
            start < X.start & end > X.end & sub_row == lead(sub_row) & end < lead(X.start) & sub_row == lag(sub_row) &
              (start == lag(X.end)| start == lag(X.end) + 1 | start == lag(X.end) -1 ) & (end == lead(X.start) | end == lead(X.start) +1 | end == lead(X.start) -1) ~ paste0('RI', ',', score, ',', 'RI'),
            start < X.start & end > X.start & end < X.end & (sub_row != lag(sub_row) | (sub_row == lag(sub_row) & start > lag(X.end)) | is.na(lag(sub_row))) ~ paste0('L', score, 'S'),
            start < X.start & end > X.start & end < X.end & sub_row == lag(sub_row) & (start == lag(X.end)| start == lag(X.end) + 1 | start == lag(X.end) -1 ) ~ paste0('RI', ',', score, 'S'),
            start > X.start & start < X.end & end > X.end & (sub_row != lead(sub_row) | (sub_row == lead(sub_row) & end < lead(X.start)) | is.na(lead(sub_row))) ~ paste0('S', score, 'L'),
            start > X.start & start < X.end & end > X.end & sub_row == lead(sub_row) & (end == lead(X.start) | end == lead(X.start) +1 | end == lead(X.start) -1) ~ paste0('S', score, ',','RI'),
            start >= X.end & (sub_row != lead(sub_row) | (sub_row == lead(sub_row) & end <= lead(X.start)) | is.na(lead(sub_row))) ~ 'NE',
            end <= X.start  & (sub_row != lag(sub_row) | (sub_row == lag(sub_row) & start >= lag(X.end) & lag(exon_name) != 'NE') | is.na(lag(sub_row))) ~ 'NE',
            TRUE ~ exon_name),
        exon_name
      )

    #Identify ciRNA cases
    bed_exon_negative <- bed_exon_negative %>%
      mutate(
        exon_name =
          case_when(
            (start == X.end | start == X.end + 1 | start == X.end - 1) & n_row == lead(n_row) &
              (end == lead(X.start) | end == lead(X.start) + 1 | end == lead(X.start) - 1) ~ paste0('I_',score),
            TRUE ~ exon_name),
        exon_name
      )

    #Label all the rest of the exons
    bed_exon_negative <- bed_exon_negative %>%
      mutate(
        exon_name =
          case_when (
            (start == X.start - 1 | start == X.start + 1 | start == X.start) &
              (end == X.end + 1 | end == X.end - 1 | end == X.end) ~ as.character(score),
            TRUE ~ exon_name),
        exon_name
      )

    #Merge transcripts with their gene names
    bed_exon_negative <- merge(bed_exon_negative,ref_gpf,by.x='feature.name',by.y='name',all.x=TRUE)
    #Rank all exons by exon index
    bed_exon_negative <- bed_exon_negative[order(bed_exon_negative$X.score),]

    #Generate circRNA UID
    df_negative <- bed_exon_negative %>%
      filter(!is.na(exon_name)) %>%
      group_by(n_row) %>%
      arrange(desc(sub_row), .by_group = TRUE) %>%
      mutate(new_exon_name = paste(exon_name, collapse = ",")) %>%
      ungroup() %>%
      distinct(n_row, .keep_all = TRUE)
    df_negative1 <- df_negative %>% select(feature.name,gene,n_row,new_exon_name)
    df_negative1$circRNA_name <- NA
    df_negative1 <- df_negative1 %>%
      mutate(circRNA_name =
               case_when(
                 !grepl('I_',new_exon_name)  ~ paste0('circ',gene,'(',new_exon_name,')'),
                 TRUE ~ gsub('I_','',paste0('ci',gene,'(',new_exon_name,')'))))


    ###Combine + & - strand###
    message('Combining strands and generating UIDs')
    df <- rbind(df_positive1,df_negative1)
    df <- df %>%
      arrange(n_row)
    annotated_bed<- merge(bed,df, by = 'n_row', all.x = TRUE)

    #Add alphabets to isoforms sharing the same UID
    message('Indexing ambiguous isoforms')
    annotated_bed$make_unique <- make.unique(annotated_bed$circRNA_name,sep = ";")
    annotated_bed$make_unique[is.na(annotated_bed$circRNA_name)] <- NA
    annotated_bed$make_unique[annotated_bed$circRNA_name==annotated_bed$make_unique] <- NA

    annotated_bed<- separate(annotated_bed,
                             col=make_unique,
                             into=c("UID","index"),
                             sep = ";",
                             remove = FALSE)
    annotated_bed$isoform_index <- LETTERS[as.numeric(annotated_bed$index)]
    annotated_bed$circRNA_UID[!is.na(annotated_bed$make_unique)] <- paste0(annotated_bed$circRNA_name[!is.na(annotated_bed$make_unique)],
                                                                           ".",
                                                                           annotated_bed$isoform_index[!is.na(annotated_bed$make_unique)])
    annotated_bed$circRNA_UID[is.na(annotated_bed$make_unique)] <- annotated_bed$circRNA_name[is.na(annotated_bed$make_unique)]


    ##Final annotated output##
    annotated_bed <- merge(annotated_bed,bed12[,c("name","n_row")],by="n_row",all.x=T)
    annotated_bed <- annotated_bed[c("name","n_row","chrom","start","end","strand","blockCount","blockSizes","blockStarts",'feature.name',"gene","new_exon_name","isoform_index","circRNA_UID")]
    names(annotated_bed) <- c("original_circRNA_name","n_row","chrom","start","end","strand","blockCount","blockSizes","blockStarts",'transcript',"gene","exon_names","isoform_index","circRNA_UID")
    message('DONE')
    annotated_bed
  } else {

    #######Case 2: BED6 circRNA coordinates as input No. 3#############
    message('BED6 circRNA coordinates as input')
    ##Load BED6 file and keep relevant columns for annotation##
    bed6 <- read.table(bed_path, header = F, sep = "\t")
    names(bed6) <- c('chrom', 'start', 'end','name', 'score', 'strand')
    bed6$n_row <- seq_len(nrow(bed6))
    bed<-bed6[,c(1:3,6)]
    names(bed) <- c('chrom','start','end','strand')

    ##Genomation annotation##
    gene.obj = readTranscriptFeatures(ref_path,
                                      up.flank=1000,
                                      down.flank=1000,
                                      remove.unusual = TRUE)
    sig.ann.obj = annotateWithGeneParts(as(bed, "GRanges"), gene.obj,strand = TRUE)
    tss = sig.ann.obj@dist.to.TSS
    bed$n_row <- seq_len(nrow(bed))
    m = as(merge(bed,tss,by.x = c("n_row","strand"), by.y = c("target.row","feature.strand"),all.x=T),"GRanges")
    m_df <- data.frame(m)
    m_df$seqnames <- as.character(m_df$seqnames)

    #Identify the indices of the exons in the annotation object that overlap the circRNAs
    hits = findOverlaps(m,gene.obj$exons,select="all")
    #Extract the indices
    idx <- unique(subjectHits(hits))
    #Get the feature names and exon numbers for the annotations that match the indices
    values <- DataFrame(gene.obj$exons[idx])
    values_df <- data.frame(values)

    bed_df <- merge(m_df,values_df,by.x='feature.name',by.y='name',all.x=TRUE)
    bed_df$abs_start_diff <- abs(bed_df$X.start - bed_df$start)
    bed_df$abs_end_diff <- abs(bed_df$X.end - bed_df$end)

    ###For + strand###
    message('Processing + strand')
    bed_df_positive <- bed_df %>%
      filter(strand == "+")%>%
      arrange(n_row,score)

    #Identify novel exons that do not match any annotated exons
    bed_df_positive$exon_name <- ifelse(
      is.na(bed_df_positive$X.seqnames) & !is.na(bed_df_positive$feature.name),
      "NE",
      NA
    )

    #Identify ciRNA cases
    bed_df_positive <- bed_df_positive %>%
      mutate(
        exon_name =
          case_when(
            (start == X.end | start == X.end + 1 | start == X.end - 1) & n_row == lead(n_row) &
              (end == lead(X.start) | end == lead(X.start) + 1 | end == lead(X.start) - 1) ~ paste0('I_',score),
            TRUE ~ exon_name),
        exon_name
      )

    bed_df_positive$start_exon_name <- NA
    bed_df_positive$end_exon_name <- NA

    #Discuss all cases of annotated exons
    bed_df_positive <- bed_df_positive %>%
      mutate(
        start_exon_name =
          case_when(
            start == X.start + 1 | start == X.start - 1 | start == X.start ~ paste0(score),
            start > X.start & start < X.end ~ paste0('S', score),
            start < X.start ~ paste0('L',score),
            TRUE ~ start_exon_name),
        start_exon_name
      )

    bed_df_positive <- bed_df_positive %>%
      mutate(
        end_exon_name =
          case_when(
            end == X.end + 1 | end == X.end - 1 | end == X.end ~ paste0(score),
            end > X.start & end < X.end ~ paste0(score,'S'),
            end > X.end ~ paste0(score,'L'),
            end < X.start & n_row == lag(n_row) & end > lag(X.end) ~ paste0(lag(score),'L'),
            TRUE ~ end_exon_name),
        end_exon_name
      )

    #Label all the rest of the exons
    bed_df_positive <- bed_df_positive %>%
      mutate(
        exon_name =
          case_when (
            (start == X.start + 1 | start == X.start - 1 | start == X.start) &
              (end == X.end + 1 | end == X.end - 1 | end == X.end) ~ paste0(score),
            TRUE ~ exon_name),
        exon_name
      )

    #Merge transcripts with their gene names
    bed_df_positive <- merge(bed_df_positive,ref_gpf,by.x='feature.name',by.y='name',all.x=TRUE)

    #Rank all exons by exon index
    bed_df_positive <- bed_df_positive[order(bed_df_positive$X.score),]


    df_positive1 <- bed_df_positive %>%
      filter(!is.na(feature.name)) %>%
      group_by(n_row) %>%
      filter(all(is.na(exon_name))) %>%
      ungroup()

    df_positive1<- df_positive1%>%
      group_by(n_row)%>%
      mutate(
        exon_name = case_when(
          !is.na(start_exon_name[which.min(abs_start_diff)]) & !is.na(end_exon_name[which.min(abs_end_diff)]) ~ {
            min_index_start <- which.min(abs_start_diff)
            min_index_end <- which.min(abs_end_diff)

            start_exon <- start_exon_name[min_index_start]
            end_exon <- end_exon_name[min_index_end]

            paste0(start_exon, ',', end_exon)
          },
          TRUE ~ exon_name
        )
      ) %>%
      ungroup() %>%
      distinct(n_row, .keep_all = TRUE)


    combine_exon_names <- function(exon_name) {
      parts<-strsplit(exon_name,',')[[1]]
      if (length(parts) ==2){
        start_exon <- parts[1]
        end_exon <- parts[2]
        start_num <- regmatches(start_exon, regexpr('\\d+', start_exon))
        end_num <- regmatches(end_exon, regexpr('\\d+',end_exon))
        if (length(start_num) > 0 && length(end_num) > 0 && start_num == end_num) {
          combined_exon <- paste0(start_exon, gsub(start_num, '', end_exon))
        } else {
          combined_exon <- paste0(exon_name)
        }
        return(combined_exon)
      } else
        return(exon_name)
    }

    df_positive1 <- df_positive1 %>%
      mutate(exon_name = sapply(exon_name,combine_exon_names))

    df_positive2 <- bed_df_positive %>%
      filter(!is.na(feature.name)) %>%
      filter(!is.na(exon_name))

    df_positive <- rbind(df_positive1, df_positive2)
    df_positive <- df_positive %>%
      mutate(exon_name = case_when(is.na(exon_name) ~ 'NE', TRUE ~ exon_name))

    df_positive <- df_positive %>%
      select(feature.name, gene,seqnames, start, end, width, strand,n_row, dist.to.feature,exon_name)

    ###For - strand###
    message('Processing - strand')
    bed_df_negative <- bed_df %>%
      filter(strand == "-")%>%
      arrange(n_row,desc(score))

    #Identify novel exons that do not match any annotated exons
    bed_df_negative$exon_name <- ifelse(
      is.na(bed_df_negative$X.seqnames) & !is.na(bed_df_negative$feature.name),
      "NE",
      NA
    )

    #Identify ciRNA cases
    bed_df_negative <- bed_df_negative %>%
      mutate(
        exon_name =
          case_when(
            (start == X.end | start == X.end + 1 | start == X.end - 1) & n_row == lead(n_row) &
              (end == lead(X.start) | end == lead(X.start) + 1 | end == lead(X.start) - 1) ~ paste0('I_',score),
            TRUE ~ exon_name),
        exon_name
      )

    bed_df_negative$start_exon_name <- NA
    bed_df_negative$end_exon_name <- NA

    #Discuss all cases of annotated exons
    bed_df_negative <- bed_df_negative %>%
      mutate(
        end_exon_name =
          case_when(
            start == X.start + 1 | start == X.start - 1 | start == X.start ~ paste0(score),
            start > X.start & start < X.end ~ paste0(score,'S'),
            start < X.start ~ paste0(score, 'L'),
            start > X.end & n_row == lead(n_row) & start < lead(X.start) ~ paste0(lead(score),'L'),
            TRUE ~ end_exon_name),
        end_exon_name
      )

    bed_df_negative <- bed_df_negative %>%
      mutate(
        start_exon_name =
          case_when(
            end == X.end + 1 | end == X.end - 1 | end == X.end ~ paste0(score),
            end > X.start & end < X.end ~ paste0('S',score),
            end > X.end ~ paste0('L',score),
            TRUE ~ start_exon_name),
        start_exon_name
      )

    #Label all the rest of the exons
    bed_df_negative <- bed_df_negative %>%
      mutate(
        exon_name =
          case_when (
            (start == X.start + 1 | start == X.start - 1 | start == X.start) &
              (end == X.end + 1 | end == X.end - 1 | end == X.end) ~ paste0(score),
            TRUE ~ exon_name),
        exon_name
      )

    #Merge transcripts with their gene names
    bed_df_negative <- merge(bed_df_negative,ref_gpf,by.x='feature.name',by.y='name',all.x=TRUE)

    #Rank all exons by exon index
    bed_df_negative <- bed_df_negative[order(bed_df_negative$X.score),]

    df_negative1 <- bed_df_negative %>%
      filter(!is.na(feature.name)) %>%
      group_by(n_row) %>%
      filter(all(is.na(exon_name))) %>%
      ungroup()

    df_negative1<- df_negative1%>%
      group_by(n_row)%>%
      mutate(
        exon_name = case_when(
          !is.na(start_exon_name[which.min(abs_end_diff)]) & !is.na(end_exon_name[which.min(abs_start_diff)]) ~ {
            min_index_start <- which.min(abs_end_diff)
            min_index_end <- which.min(abs_start_diff)

            start_exon <- start_exon_name[min_index_start]
            end_exon <- end_exon_name[min_index_end]

            paste0(start_exon, ',', end_exon)
          },
          TRUE ~ exon_name
        )
      ) %>%
      ungroup() %>%
      distinct(n_row, .keep_all = TRUE)


    combine_exon_names <- function(exon_name) {
      parts<-strsplit(exon_name,',')[[1]]
      if (length(parts) ==2){
        start_exon <- parts[1]
        end_exon <- parts[2]
        start_num <- regmatches(start_exon, regexpr('\\d+', start_exon))
        end_num <- regmatches(end_exon, regexpr('\\d+',end_exon))
        if (length(start_num) > 0 && length(end_num) > 0 && start_num == end_num) {
          combined_exon <- paste0(start_exon, gsub(start_num, '', end_exon))
        } else {
          combined_exon <- paste0(exon_name)
        }
        return(combined_exon)
      } else
        return(exon_name)
    }

    df_negative1 <- df_negative1 %>%
      mutate(exon_name = sapply(exon_name,combine_exon_names))

    df_negative2 <- bed_df_negative %>%
      filter(!is.na(feature.name)) %>%
      filter(!is.na(exon_name))

    df_negative <- rbind(df_negative1, df_negative2)
    df_negative <- df_negative %>%
      mutate(exon_name = case_when(is.na(exon_name) ~ 'NE', TRUE ~ exon_name))

    df_negative <- df_negative %>%
      select(feature.name, gene,seqnames, start, end, width, strand,n_row, dist.to.feature,exon_name)

    ###Combine + & - strand###
    message('Combining strands and generating UIDs')
    df <- rbind(df_positive, df_negative)
    df <- df %>%
      select(n_row, width, feature.name, gene,dist.to.feature, exon_name) %>%
      arrange(n_row)
    df$circRNA_name <- NA
    df <- df %>%
      mutate(circRNA_name =
               case_when(
                 !grepl('I_',exon_name)  ~ paste0('circ',gene,'(',exon_name,')'),
                 TRUE ~ gsub('I_','',paste0('ci',gene,'(',exon_name,')'))))

    annotated_bed<- merge(bed,df, by = 'n_row', all.x = TRUE)

    #Add alphabets to isoforms sharing the same UID
    message('Indexing ambiguous isoforms')
    annotated_bed$make_unique <- make.unique(annotated_bed$circRNA_name,sep = ";")
    annotated_bed$make_unique[is.na(annotated_bed$circRNA_name)] <- NA
    annotated_bed$make_unique[annotated_bed$circRNA_name==annotated_bed$make_unique] <- NA

    annotated_bed<- separate(annotated_bed,
                             col=make_unique,
                             into=c("UID","index"),
                             sep = ";",
                             remove = FALSE)
    annotated_bed$isoform_index <- LETTERS[as.numeric(annotated_bed$index)]
    annotated_bed$circRNA_UID[!is.na(annotated_bed$make_unique)] <- paste0(annotated_bed$circRNA_name[!is.na(annotated_bed$make_unique)],
                                                                           ".",
                                                                           annotated_bed$isoform_index[!is.na(annotated_bed$make_unique)])
    annotated_bed$circRNA_UID[is.na(annotated_bed$make_unique)] <- annotated_bed$circRNA_name[is.na(annotated_bed$make_unique)]

    ##Final annotated output##
    annotated_bed <- merge(annotated_bed,bed6[,c("name","n_row")],by="n_row",all.x=T)
    annotated_bed <- annotated_bed[c("name","n_row","chrom","start","end","strand",'feature.name',"gene","exon_name","isoform_index","circRNA_UID")]
    names(annotated_bed) <- c("original_circRNA_name","n_row","chrom","start","end","strand",'transcript',"gene","exon_names","isoform_index","circRNA_UID")
    message('DONE')
    annotated_bed
  }
}

