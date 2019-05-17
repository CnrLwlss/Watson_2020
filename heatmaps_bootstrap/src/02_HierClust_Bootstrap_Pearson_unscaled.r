library(pvclust)

path = "../../data/"
for (i in list.files(path)) {
  dataset <- read.csv(paste(path,i, sep=""))
  colnames(dataset)[1]="X"
  dataset <- dataset[!duplicated(dataset$X),]
  row.names(dataset) <- dataset$X
  dataset$X <- NULL
  dataset <- as.matrix(dataset)
  len = length(unlist(strsplit(i,"\\_")))
  result <- pvclust(dataset, method.dist=function(x)as.dist(1-cor(x, use="pairwise", method="pearson")), method.hclust="ward.D2", nboot=10000, parallel=TRUE)
  svg(file=paste0("../Report",i,"_Bootstrap_Pearson_unscaled_10000.svg"),width=ifelse(dim(dataset)[2] > 300, 28.5 , ifelse(dim(dataset)[2] > 200, 18.5 , 9.25 )),height=6.5)
  plot(result, cex = 0.4, cex.pv=0.4)
  pvrect(result, alpha=0.95)
  dev.off()
}

