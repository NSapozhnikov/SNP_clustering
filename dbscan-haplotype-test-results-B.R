## Process the results of haplotype test for dataset B to subset the significant regions 
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(ggplot2))
source("helpers.R")


# Initiate
folder <- "data2"
files <- list.files(folder, "chr[0-9]+_nonOmnibus.csv", full.names = T)

# Load data
d1 <- plyr::ldply(files, data.table::fread)
# Add genomic coordinates, haplotype size and length
d2 <- d1 %>%  
  mutate(CHR = as.integer(stringr::str_match(LOCUS, "([0-9]+):")[, 2]),
         START = as.integer(stringr::str_match(LOCUS, ":([0-9]+)-")[, 2]),
         END = as.integer(stringr::str_match(LOCUS, "-([0-9]+)$")[, 2])) %>%
  mutate(SIZE = nchar(HAPLOTYPE), LENGTH = END - START)
# Clean
rm(d1)

# Get the number of clusters per chromosomes 
nc <- d2 %>% group_by(CHR) %>% group_map(~ GetNumClusters(.x), .keep = T) %>% unlist()
# Plot the number of clusters per chromosome
nr <- table(d2$CHR)
#barplot(nc, names.arg = names(nr), xlab = "хромосомы", ylab = "Число кластеров")

# Get the number of clusters
n <- sum(nc)
# Get threshold value of p-value
p <- 0.05/n
# Subset haplotypes with p-value less then threshold one 
h1 <- d2 %>% filter(P <= p)
rm(d2)
# Reduce haplotypes
h2 <- ReduceHaps(h1)

# Compare frequencies
df1 <- reshape2::melt(h2[, c("F_A", "F_U")])
ggplot(df1) + geom_boxplot(aes(y = value, x = variable)) + 
  labs(y = "Частота", x = "") +
  theme(axis.text=element_text(size = 12),
        axis.title=element_text(size = 14))

# Identify outlier frequecies
boxplot.stats(h2$F_A)
boxplot.stats(h2$F_U)

# Subset haplotypes by frequency
h3 <- h2 %>% filter(F_A >= 0.107 | F_U >= 0.107)

# Plot the distribution of selected haplotypes 
barplot(table(h3$CHR), xlab = "Хромосома", ylab = "Число гаплотипов")

plot(h3$CHR, h3$F_A, xaxt = "n", xlab = "Хромосома", ylab = "Частота", )
axis(1, at = 1:22, labels = 1:22)

# Save
write.csv(h3, "out/significant-haplotypes-04.csv", row.names = F)


