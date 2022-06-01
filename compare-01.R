## Compare the significant locuses obtained for datasets A, B, and GWAS made by S.Nikolich

# Initiate
files <- sprintf("/home/gennady/proj/nsap/out/significant-haplotypes-0%s.csv", 1:2)

# Load data
l1 <- plyr::llply(files, read.csv)


# Convert into GRanges
gr1 <- plyr::llply(l1, function(x) {
  GenomicRanges::GRanges(seqnames = x$CHR,
                         ranges = IRanges::IRanges(start = x$START, end = x$END),
                         seqinfo = paste(1:22))
  })  

# Overlap GRanges
ov <- GenomicRanges::countOverlaps(gr1[[1]], gr1[[2]], type = "any")
