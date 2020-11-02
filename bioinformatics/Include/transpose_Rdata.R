library(data.table)

args <- commandArgs(TRUE)
folder <- args[1]
folder
input <- paste(args[1],"/",args[2],".PCA_normalized.filtered.sample_zscores.RD.txt",sep="")
input
samples <- c(strsplit(args[3],","))
samples <- as.list(samples)[[1]]
samples
output <- paste(args[1],"/",args[2],".logratio.txt",sep="")
output
data <- fread(input, header=F)

data_transpose <- transpose(data)
names(data_transpose) <- as.character(unlist(data_transpose[1,]))
data_transpose <- data_transpose[-1,]
data_starseq <- subset(data_transpose,select=data_transpose[,which(colnames(data_transpose) %in% samples)])
write.table(data_starseq, output, quote=F, row.names=F, sep="\t")
print("general logratio file ok")

### for each sample of the run, export data : chr:start-end\tlog-ratio => sample.txt
for (x in samples){
    if(x %in% colnames(data_transpose) && x != "Matrix"){
	out <- file.path(folder, paste0("logratio_perSample/",x,".logratio.txt"))
	print(out)
	data_write <- data_transpose[,c("Matrix",x),with=FALSE]
    	write.table(data_write,file=as.character(out),quote=FALSE,row.names=FALSE,sep="\t")
	}
}
