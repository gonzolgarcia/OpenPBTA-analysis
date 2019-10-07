## SOMATIC FILTERING FOR LUMPY
# assuming the null p(ref) = p(alt) = 0.5
# 
# minimum tumor alt count (given normal alt = 0) that satisfies p binom < 0.01; see expr below
# => binom.test(0,7,0.5,alternative="less")$p.value  

# in case normal alt > 0 we calculate the p-val for the binomial test 

# in case normal alt >= 5 we assume is not somatic


fileName <- "data/release-v5-20190924/pbta-sv-lumpy.tsv.gz"
conn <-  gzfile(fileName,open="r")
linn <- readLines(conn)

b <- sapply(1:length(linn), function(i)  strsplit(linn[i],"\t")[[1]][14:15] )	# get tumor and normal read counts
normal <- as.numeric(sapply(b[1,], function(x) strsplit(x,":")[[1]][2]))
tumor <- as.numeric(sapply(b[2,], function(x) strsplit(x,":")[[1]][2]))

# prefilter based on minimun alt tumor (7) and maximum alt normal (5)
somatic_filter_1 <- intersect(which(normal < 5),which(tumor >= 7))

# For the remaining cases do binom test
tumor_alt_count <- tumor[somatic_filter_1]	# tumor alt count 
normal_alt_count <- normal[somatic_filter_1]	# normal alt count
alt_count <- tumor_alt_count + normal_alt_count
binom <- sapply(1:length(tumor_alt_count), function(i) binom.test(normal_alt_count[i],alt_count[i],0.5,alternative="less")$p.value)

# subselect cases that pass binom test
somatic <- somatic_filter_1[which(binom < 0.01)]

# write file including header 
write.table(linn[c(1,somatic)],file="data/release-v5-20190924/pbta-sv-lumpy-somatic.tsv",col.names=F,quote=F,row.names=F)


