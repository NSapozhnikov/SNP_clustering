## Analysis of haplotype test results for dataset A
## Clustering: DBSCAN
suppressPackageStartupMessages(require(dplyr))
source("helpers.R")

# Initiate
files <- c("/home/nsap/hap_assoc/out_merged_sorted_nonOmnibus.xlsx",
           "/home/nsap/hap_assoc/out_merged_sorted_omnibus.xlsx")

# Check files exists
plyr::l_ply(files, function(x) if(!file.exists(x)) stop(x, "doesn't exist", call. = F))

## Analyse the results of individual haplotype testing

# Load data
d1 <- readxl::read_xlsx(files[1])
# Add explicit genomic coordinates and cluster properties
d2 <- d1 %>% 
  mutate(CHR = as.integer(stringr::str_match(LOCUS, "([0-9]+):")[, 2]),
         START = as.integer(stringr::str_match(LOCUS, ":([0-9]+)-")[, 2]),
         END = as.integer(stringr::str_match(LOCUS, "-([0-9]+)$")[, 2])) %>%
  mutate(SIZE = nchar(HAPLOTYPE), LENGTH = END - START)

# Get the number of clusters per chromosomes 
nc <- d2 %>% group_by(CHR) %>% group_map(~ GetNumClusters(.x), .keep = T) %>% unlist()
# Plot the number of clusters per chromosome
nr <- table(d2$CHR)
barplot(nc, names.arg = names(nr), xlab = "хромосомы",
        ylab = "Число кластеров")

# Count p-value threshold
padj <- 0.05/sum(nc)

# Subset haplotypes
h1 <- d2 %>% filter(P <= padj)

# Reduce haplotypes
h2 <- ReduceHaps(h1)

# Save
write.csv(h2, "out/significant-haplotypes-03.csv", row.names = F)


