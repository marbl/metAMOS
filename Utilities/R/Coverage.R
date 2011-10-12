# compute a few values for:
# all data
# 90% 
# 80%
# 70%
# ...
#
# The values computed
# number contigs > 3X coverage
# total size of contigs > 3X coverage
# max contig size
# min contig size
# mean contig size
# median contig size
# size @ 1Mbp
# size @ 2Mbp
# size @ 4Mbp
# size @ 10Mbp

data <- read.table("SRS048791_LANL.cov", sep='\t', header=F)

out <- array(0,c(20,7))
xval <- c(1:20)
xval <- xval * 100/20
for (i in 1:20) {
  out[i,] <- asmStats(data$V2[data$V4*i/20 > 3])
#  colnames(out) <- colnames(asmStats(data$V2[data$V4*i/10 > 3]))
}

par(mfrow=c(2,2), mar=c(4,5,2,1))
plot(xval,out[,2]/1000000, type="l", xlab="Percent reads", ylab="Total Size (Mbp)", ylim=c(0,100), cex.lab=2, cex.axis=2)
plot(xval,out[,3], type="l", xlab="Percent reads", ylab="Max contig (bp)", cex.axis=2, cex.lab=2)
plot(xval,out[,6], type="l", xlab="Percent reads", ylab="Median contig (bp)", cex.axis=2, cex.lab=2)
plot(xval,out[,7], type="l", xlab="Percent reads", ylab="Size @ 10Mbp (bp)", cex.axis=2, cex.lab=2)


asmStats <- function(arr)
{

  outA <- array(0, c(0,7)) 
  colnames(outA) <- list("Number", "TotSize", "MaxContig", "MinContig", "MeanContig", "MedianContig", "Size10M")
  outA["Number"] <- length(arr)
  outA["TotSize"] <- sum(arr)
  outA["MaxContig"] <- max(arr)
  outA["MinContig"] <- min(arr)
  outA["MeanContig"] <- mean(arr)
  outA["MedianContig"] <- median(arr)
  outA["Size10M"] <- nSize(arr, 10000000)
 outA
}

nSize <- function(arr, size)
{
   sa <- sort(arr, dec=T)
   sums <- cumsum(sa)
   small <- sums[sums<size]
   sa[length(small)+1]
}