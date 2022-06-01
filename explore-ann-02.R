suppressPackageStartupMessages(require(dplyr))

# Initiate
wd <- "/home/gennady/data/stroke/nsap"
files <- sprintf("%s/regions/significant-haplotypes-0%s.eff.tab", wd, 3:4)
names(files) <- c("A", "B")
# Load 
l1 <- sapply(files, function(x) read.table(x, stringsAsFactors = F), simplify = F)
# Remove duplicated 
l2 <- sapply(l1, function(x) x %>% filter(!duplicated(x)), simplify = F)

# Count unique genes
sapply(l2, function(x) length(unique(x$V4)))
# Save the lists of unique genes
mapply(function(x, y) {
  f <- sprintf("%s/regions/genes-0%s.txt", wd, y)
  write(unique(l2[[x]]$V4), f)
}, 1:2, 3:4)
# Intersect the lists from two datasets
ggvenn::ggvenn(list(A = unique(l2[[1]]$V4), B = unique(l2[[2]]$V4)))









