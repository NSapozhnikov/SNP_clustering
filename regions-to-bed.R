# Initiate
files <- sprintf("/home/gennady/data/stroke/nsap/regions/significant-haplotypes-0%s.csv", 1:4)

# Load 
l <- sapply(files, read.csv)

# Save
plyr::l_ply(1:4, function(x) {
  write.table(l[[x]][, c("CHR", "START", "END")], 
              sprintf("/home/gennady/data/stroke/nsap/regions/significant-haplotypes-0%s.bed", x),
              sep = "\t", row.names = F, col.names = F)
  })
