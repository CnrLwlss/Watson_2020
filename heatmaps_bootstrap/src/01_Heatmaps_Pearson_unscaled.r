library(ComplexHeatmap)
library(circlize)

path = "../../data/"
for (i in list.files(path)) {
dataset <- read.csv(paste(path,i, sep=""))
colnames(dataset)[1]="X"
dataset <- dataset[!duplicated(dataset$X),]
row.names(dataset) <- dataset$X
dataset$X <- NULL
dataset <- as.matrix(dataset)

if (grepl("DSF", i) == TRUE){
  color = colorRamp2(c(-20, 0, 20), c("red", "white", "blue"))
}else {
  color = colorRamp2(c(-100, 0, 100), c("red", "white", "blue"))  
}
len = length(unlist(strsplit(i,"\\_")))
svg(file=paste0("../Report",i,"_Pearson_unscaled.svg"),
    width=ifelse(dim(dataset)[2] > 300, 26.5 , ifelse(dim(dataset)[2] > 200, 16.5 , 7.25 )),
    height=ifelse(dim(dataset)[2] > 100, 8, 5))
print(Heatmap(dataset, name="Value",
        col = color,
        na_col = "grey",
        clustering_distance_rows = function(x)as.dist(1-cor(t(x), use="pairwise", method="pearson")), 
        clustering_distance_columns = function(x)as.dist(1-cor(t(x), use="pairwise", method="pearson")), 
        clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2",
        show_row_names = FALSE,
        column_title = paste0(unlist(strsplit(i,"\\_"))[1:len-1], collapse=" "),
        column_title_gp = gpar(fontsize=10, fontface="bold"),
        column_names_gp = gpar(fontsize=5),
        column_dend_height = unit(30, "mm")
))

dev.off()
}

