
# ==========================================
# 1. DNA Encoding
# ==========================================

#' DNA Sequence Encoding
#'
#' This function converts DNA sequences into
#' individual nucleotide position features.
#'
#' @param dna_strings A character vector containing DNA sequences.
#'
#' @return A data frame containing nucleotide positions
#' as categorical variables.

dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")),
                  ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}


# ==========================================
# 2. Multiple Sample Prediction
# ==========================================

#' Multiple Sample Prediction for m6A Sites
#'
#' This function predicts m6A probabilities and
#' classification status for multiple samples.
#'
#' @param ml_fit A trained random forest model.
#' @param feature_df A data frame containing input features.
#' @param positive_threshold A probability threshold (default 0.5).
#'
#' @return A data frame containing the original features,
#' predicted m6A probabilities and predicted m6A status.
#'
#' @examples
#' ml_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds", package = "m6APrediction"
#' ))
#' example_df <- read.csv(system.file(
#'   "extdata", "m6A_input_example.csv",
#'   package = "m6APrediction"
#' ))
#' prediction_multiple(ml_fit, example_df, positive_threshold = 0.6)
#'
#' @export
#' @import randomForest

prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df)))

  feature_df$RNA_type <- factor(as.character(feature_df$RNA_type),
                                levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))

  feature_df$RNA_region <- factor(as.character(feature_df$RNA_region),
                                  levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  encoded_dna <- dna_encoding(as.character(feature_df$DNA_5mer))

  model_df <- data.frame(
    gc_content = feature_df$gc_content,
    RNA_type = feature_df$RNA_type,
    RNA_region = feature_df$RNA_region,
    exon_length = feature_df$exon_length,
    distance_to_junction = feature_df$distance_to_junction,
    evolutionary_conservation = feature_df$evolutionary_conservation,
    encoded_dna,
    check.names = FALSE
  )

  probability <- predict(ml_fit, newdata = model_df, type = "prob")[, "Positive"]

  feature_df$predicted_m6A_prob <- probability

  feature_df$predicted_m6A_status <- ifelse(
    probability > positive_threshold,
    "Positive", "Negative"
  )

  return(feature_df)
}


# ==========================================
# 3. Single Sample Prediction
# ==========================================

#' Single Sample Prediction for m6A Sites
#'
#' This function predicts the m6A probability and
#' classification status for an individual sample.
#'
#' @param ml_fit A trained random forest model.
#' @param gc_content A numeric value for GC content.
#' @param RNA_type A character string specifying RNA type.
#' @param RNA_region A character string specifying RNA region.
#' @param exon_length A numeric value for exon length.
#' @param distance_to_junction A numeric value for distance to junction.
#' @param evolutionary_conservation A numeric conservation value.
#' @param DNA_5mer A character string representing a DNA 5-mer.
#' @param positive_threshold A probability threshold (default 0.5).
#'
#' @return A named vector containing predicted m6A probability
#' and predicted m6A status.
#'
#' @examples
#' ml_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds", package = "m6APrediction"
#' ))
#' prediction_single(
#'   ml_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "GGACA",
#'   positive_threshold = 0.5
#' )
#'
#' @export

prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){

  feature_df <- data.frame(
    gc_content = gc_content,
    RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
    RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer
  )

  prediction <- prediction_multiple(ml_fit, feature_df, positive_threshold)

  returned_vector <- c(
    predicted_m6A_prob = as.character(prediction$predicted_m6A_prob[1]),
    predicted_m6A_status = as.character(prediction$predicted_m6A_status[1])
  )

  return(returned_vector)
}

