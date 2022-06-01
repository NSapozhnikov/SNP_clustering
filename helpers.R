## Functions
suppressPackageStartupMessages(require(dplyr))

ReduceHaps <- function(d4) {
  # Drop LOCUS column
  d4$LOCUS <- NULL
  # Create combinations of indexes
  n <- nrow(d4)
  ind1 <- combn(n, 2, simplify = F)
  # Iterate over the indexes and setup the rows to be removed
  ind2 <- sapply(ind1, function(x) {
    h1 <- d4 %>% slice(x[1]) %>% unlist()
    h2 <- d4 %>% slice(x[2]) %>% unlist()
    if(h1["CHR"] == h2["CHR"] & h1["START"] == h2["START"]) {
      ifelse(h1["SIZE"] < h2["SIZE"], return(x[1]), return(x[2]))
    } 
  }) %>% unlist() %>% unique()
  # Remove haplotypes of less size and order by genomic coordinates
  d5 <- d4 %>% slice(-ind2) %>% arrange(CHR, START)
  d5
}

GetNumClusters <- function(x) {
  # Group by CHR, START 
  dt <- x %>% count(CHR, START)
  # Get the number of clusters
  n <- nrow(dt)
  n
}

SubsetHaplotypesByPvalue <- function(x) {
  # Group by CHR, START 
  d1 <- x %>% count(CHR, START)
  # Get the number of clusters
  n <- nrow(d1)
  # Get threshold value of p-value
  p <- 0.05/n
  # Subset haplotypes with p-value less then threshold one
  x %>% filter(P < p)
}

Consolidate <- function(files) {
  ## Consolidate *sign.haps.csv files
  ## Args:
  ##  files: a vector with paths to *.csv files created with Plink
  ## Value:
  ##  A data frame with all haplotypes.
  # Load data
  d1 <- lapply(files, read.csv)
  d2 <- do.call("rbind", d1)
  # Extract the genomic coordinates
  d3 <- data.frame(
    CHR = stringr::str_match(d2$LOCUS, "([0-9]+):")[, 2] %>% as.integer(),
    START = stringr::str_match(d2$LOCUS, ":([0-9]+)-")[, 2] %>% as.integer(),
    END = stringr::str_match(d2$LOCUS, "-([0-9]+)")[, 2] %>% as.integer())
  # Create data frame with genomic coordinates
  d2$LOCUS <- NULL
  d4 <- cbind(d3, d2)
  # Create combinations of indexes
  n <- nrow(d4)
  ind1 <- combn(n, 2, simplify = F)
  # Iterate over the indexes and setup the rows to be removed
  ind2 <- sapply(ind1, function(x) {
    h1 <- d4 %>% slice(x[1]) %>% unlist()
    h2 <- d4 %>% slice(x[2]) %>% unlist()
    if(h1["CHR"] == h2["CHR"] & h1["START"] == h2["START"]) {
      ifelse(h1["SIZE"] < h2["SIZE"], return(x[1]), return(x[2]))
    } 
  }) %>% unlist() %>% unique()
  # Remove haplotypes of less size and order by genomic coordinates
  d5 <- d4 %>% slice(-ind2) %>% arrange(CHR, START)
  d5
}

