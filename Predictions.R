library(parallel)

ctest = function(A,B){
  res = try(cor.test(as.numeric(A),as.numeric(B),use="pairwise.complete.obs"))
  if(class(res)=="try-error"){
   res = cor.test(1:5,6:10)
   res$p.value = 1
   res$estimate[["cor"]]=cor(A,B,use="pairwise.complete.obs")
  }
  return(res)
}

correlations = function(dat){
 getcorrs = lapply(genes,function(x) ctest(dat$UNKNOWN,dat[[x]]))
 cors = NULL
 pvals = NULL
 for(i in seq_along(getcorrs)){
  cors[i]=getcorrs[[i]]$estimate[[1]]
  pvals[i]=getcorrs[[i]]$p.value
 }
 df = data.frame(genes=genes,cors=cors,pvals=pvals,stringsAsFactors=FALSE)
 df = df[order(df$cors,decreasing=TRUE),]
 return(df)
}

ranks = function(dat) correlations(dat)$genes

getfrac = function(dat,hits,nsamp,nreps=1000,ntop=10){
  NKinase = length(dat$INHIBITOR)
  corrmat = replicate(nreps, ranks(dat[sample(1:NKinase,nsamp,replace=FALSE),]))
  res = try(mean(sapply(1:nreps,function(x) sum(hits%in%corrmat[1:ntop,x])==length(hits))))
  if(class(res) == "try-error") {res = -99}
  return(res)
}

readDat = function(fpath){
 dat = read.delim(fpath, stringsAsFactors=FALSE, sep=ifelse(grepl(".csv",fpath),",","\t"))
 dat$X = NULL
 colnames(dat)[1] = "INHIBITOR"
 colnames(dat)[2] = "UNKNOWN"
 dat = dat[!is.na(dat$UNKNOWN),]
 return(dat)
}

getnticks = function(vals,nmax=20){
 x = vals-1
 y = seq_len(x)
 divs = y[ x%%y == 0 ]
 nticks = x/divs

 nticks[which(nticks<nmax, arr.ind=TRUE)[1]]
}

fpaths = c(
"data\\INCENP S446 Nanosyn.txt",
"data\\Integrin pY Nanosyn.txt",
"data\\EGFR Nanosyn.txt",
"data\\H3T3ph DSF.txt",
"data\\H3T3ph Harvey 1uM.txt",
"data\\H3T3ph Peterson.txt",
"data\\INCENP S446 DSF.txt"
)

fpaths = c(
"data\\INCENP S446 Nanosyn degapped.csv",
"data\\Integrin pY Nanosyn degapped.csv",
"data\\EGFR Nanosyn degapped.csv",
"data\\H3T3ph DSF degapped.csv",
"data\\H3T3ph Harvey 1uM degapped.csv",
"data\\H3T3ph Peterson degapped.csv",
"data\\H3T3ph Nanosyn degapped.csv",
"data\\INCENP S446 DSF degapped.csv"
)

fpaths = c(
"data\\H3T3phPeterson.csv",
"data\\H3T3phHarvey1.csv",
"data\\H3T3phDSF.csv",
#"data\\H3S28phZarrinkar.csv",
"data\\H3S28phZarrinkarpt5.csv",
"data\\H3S28phPeterson.csv",
"data\\H3S28phNanosynHi.csv",
"data\\EGFRNanosynLo.csv",
"data\\EGFRNanosynHi.csv"
)

# Had to delete spurious column at end of "data\\H3S28phNanosynHi.csv"

#truehits = c(7,13,3,1,1,1,1,2)
truehits = c(1,1,1,1,1,1,1,1)
names(truehits) = fpaths

# First, just check hits
for(fpath in fpaths){
 ftitle = strsplit(basename(fpath),"\\.")[[1]][1]   
 dat = readDat(fpath)
 genes = colnames(dat)[3:length(colnames(dat))]
 corrs = correlations(dat)
 print(ftitle)
 print(head(corrs,20))
}

no_cores <- detectCores()/2 -1
clust <- makeCluster(no_cores)

#cutoff = 0.4
#nhits = 3
nreps = 10000
tol = 0.99

fres = list()

for(fpath in fpaths){
 ftitle = strsplit(basename(fpath),"\\.")[[1]][1]   
 dat = readDat(fpath)
 genes = colnames(dat)[3:length(colnames(dat))]
 corrs = correlations(dat)
 tops = ceiling(length(corrs$genes)*c(0.01,0.05,0.1))
 percents = tops/length(corrs$genes)
 print(fpath)
 print(tops)
 print(percents)
 

 #params = data.frame(hits = c(3,3,3,3,3),tops = c(3,4,5,10,15))
 #params = data.frame(hits = c(1,3,3,3),tops = c(10,3,5,10))
 #params = data.frame(hits = rep(truehits[fpath],4), tops = truehits[fpath]+(0:3)*2)
 #params = data.frame(hits = truehits[fpath], tops = truehits[fpath])
 params = data.frame(hits = c(1,1,1), tops = c(1,5,10)) 
 #params = data.frame(hits = c(1,1,1), tops = tops)

 result = list()

 clusterExport(clust,list("getfrac", "ranks", "correlations","corrs","genes","nreps","dat","params","genes","ctest"))

 for(i in 1:length(params[,1])){
   newn = sprintf("fracs_%02d_%02d",params[i,"hits"],params[i,"tops"])
   hits = corrs$genes[1:params$hits[i]]
   #sampsizes = max(length(hits),3):length(dat$INHIBITOR)
   sampsizes = round(seq(max(length(hits),3),length(dat$INHIBITOR),length=50))
   clusterExport(clust,list("hits","i","sampsizes"))
   result[[newn]] = parSapply(clust, sampsizes, function(x) getfrac(dat,hits,nsamp=x,nreps=nreps,ntop=params$tops[i]))
   #result[[newn]] = sapply(sampsizes, function(x) getfrac(dat,hits,nsamp=x,nreps=nreps,ntop=params$tops[i]))
 }
 fres[[ftitle]][["result"]] = result
 fres[[ftitle]][["corrs"]] = corrs
 fres[[ftitle]][["dat"]] = dat
 fres[[ftitle]][["sampsizes"]] = sampsizes
 fres[[ftitle]][["truehits"]] = truehits[fpath]
 fres[[ftitle]][["params"]] = params

}

stopCluster(clust)

pdf("Sampling_fixed.pdf", width=11.69,height=8.27)
for(ftitle in names(fres)){

 result = fres[[ftitle]][["result"]]
 corrs = fres[[ftitle]][["corrs"]]$cors
 pvals = fres[[ftitle]][["corrs"]]$pvals
 dat = fres[[ftitle]][["dat"]]
 sampsizes = fres[[ftitle]][["sampsizes"]]
 truehitz = fres[[ftitle]][["truehits"]]

 #tops = ceiling(length(fres[[ftitle]][["corrs"]]$genes)*c(0.01,0.05,0.1))
 tops = fres[[ftitle]][["params"]]$tops
 percents = tops/length(fres[[ftitle]][["corrs"]]$genes)
 params = data.frame(hits = c(1,1,1), tops = tops)

 # Piecewise linear nature of plot suggests that correlations above 0.25 look interesting
 #hits = names(corrs[corrs > cutoff])
 nhits = truehitz
 cutoff = (corrs[nhits] + corrs[nhits+1])/2

 op=par(mfrow=c(1,2))
 mlab = paste(ftitle,paste(length(dat$INHIBITOR),"inhibitors",dim(dat)[2],"genes"),sep="\n")

 plot(corrs,type="n",xlab="Query kinase rank by correlation",ylab="Correlation with unknown kinase",main=mlab,ylim=range(corrs))
 abline(h = c(cutoff,0,-cutoff), col="red",lwd=3,lty=c(2,1,2))
 points(corrs,type="b")
 text(1:nhits,corrs[1:nhits],fres[[ftitle]][["corrs"]]$genes[1:nhits],pos=4,cex=0.5)

 plot(NULL,xlab="Size of sample",ylab=paste("Fraction of hits detected in",nreps,"samples"),type="l", main=paste("Random samples from",ftitle),
 xlim=c(0,round(ceiling(length(dat$INHIBITOR)+5),5)),
 ylim=c(0,1.0)
 #xaxp  = c(1, length(dat$INHIBITOR), getnticks(length(dat$INHIBITOR),nmax=30))
 )
 print(ftitle)
 print(paste("Size of full library:",max(sampsizes)))
 for(i in seq_along(names(result))) {
   n = names(result)[i]
   points(sampsizes,result[[n]],type="l",col=i,lwd=2)
   lastbelow = sampsizes[length(result[[n]]) - min(which(rev(result[[n]])<tol)) + 1]
   #AUC = sum(result[[n]][result[[n]]>0])/max(sampsizes)
   AUC = sum(filter(result[[n]][result[[n]]>0],c(0.5,0.5),sides=2)*filter(sampsizes[result[[n]]>0],c(1,-1),sides=2),na.rm=TRUE)/max(sampsizes)
   print(n)
   print(paste("Last library size below tolerance:",lastbelow))
   print(paste("AUC:",formatC(AUC,4)))
 }
 print("")

 legend("bottomright",legend=paste(rep(truehitz,length(params$tops)),params$tops,sep="/"),col=1:length(params$hits),lwd=2)
 par(op)
}


dev.off()

clist=list()
for(ftitle in sort(names(fres))){
  print(ftitle)
  print(fres[[ftitle]][["corrs"]]$genes[1:20])
  cvals = fres[[ftitle]][["corrs"]]$cors
  zvals = (cvals-(mean(cvals)))/sd(cvals)

  clist[[ftitle]] = data.frame(corr=cvals,label= sub(" ","\n",gsub(" degapped","",ftitle)),hit=FALSE,z=FALSE,stringsAsFactors=FALSE)
  clist[[ftitle]]$hit[1:fres[[ftitle]][["truehits"]]]=TRUE
  clist[[ftitle]]$z=zvals
  clist[[ftitle]]$p=fres[[ftitle]][["corrs"]]$pvals


}
corrdf = do.call("rbind",clist)
liblabs = sort(unique(corrdf$label))
nums = 1:length(liblabs)
names(nums) =liblabs
corrdf$num = nums[corrdf$label]
corrdf$jit = runif(length(corrdf$num),0,0.3)

op=par(mfrow=c(1,2),mar=c(10, 4, 4, 2) + 0.1)

pcut=0.001
plot(corrdf$num+corrdf$jit,corrdf$corr,pch=16,col=rgb(0,0,0,0.2),ylab="Correlation",xlab="",axes=FALSE)
points((corrdf$num+corrdf$jit)[corrdf$hit],corrdf$corr[corrdf$hit],pch=16,col="red",cex=0.75)
#points((corrdf$num+corrdf$jit)[abs(corrdf$z)>zcutoff],corrdf$corr[abs(corrdf$z)>zcutoff],pch=16,col="blue",cex=0.25)
points((corrdf$num+corrdf$jit)[corrdf$p<pcut],corrdf$corr[corrdf$p<pcut],pch=16,col="blue",cex=0.25)
#legend("bottomleft",c("
axis(2)
axis(1,at=1:length(liblabs),labels=names(nums),las=2)

plot(corrdf$num+corrdf$jit,corrdf$z,pch=16,col=rgb(0,0,0,0.2),ylab="z-score",xlab="",axes=FALSE)
points((corrdf$num+corrdf$jit)[corrdf$hit],corrdf$z[corrdf$hit],pch=16,col="red",cex=0.75)
points((corrdf$num+corrdf$jit)[abs(corrdf$z)>zcutoff],corrdf$z[abs(corrdf$z)>zcutoff],pch=16,col="blue",cex=0.25)
axis(2)
axis(1,at=1:length(liblabs),labels=names(nums),las=2)

par(op)
