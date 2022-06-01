## Consolidate the significant haplotypes
source("helpers.R")

# Get the list of files from dataset B
files <- list.files("data2", "*sign.haps.csv", full.names = T)
out <- Consolidate(files)
# Save 
write.csv(out, "out/significant-haplotypes-02.csv", row.names = F)
# Explore the results
summary(out)

# Visulize
old <- par(mfrow = c(1,2))
hist(out$SIZE, breaks = 20, main = "", xlab = "Размер гаплотипа",
     ylab = "Количество")
hist(out$LENGTH, breaks = 20, main = "", xlab = "Длина гаплотипа, п.н.",
     ylab = "Количество")
par(old)
