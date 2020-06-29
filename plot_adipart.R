plot_adipart<-function(adipart.obj,label="model",terms){
  nlevels<-attr(adipart.obj,"n.levels")
  indices<-c(1,(1:(nlevels-1))+nlevels)
  hierarchies<-rev(names(adipart.obj$statistic)[indices])
  if(label=="model"){
    labels<-attr(adipart.obj,"terms")[1:(nlevels-1)]
  } else {
    labels<-terms
  }
  hierarchies<-paste(hierarchies,rev(c("within",rep("between",nlevels-1))),sep=".")
  hierarchies<-sapply(strsplit(hierarchies,".",fixed=T),function(x){
    paste(x[1],x[3],labels[as.numeric(x[2])])
  })
  barwidth=0.6
  lines<-data.frame(x=(1-barwidth)*rep(1,nlevels)+1.5*barwidth, xend=(1-barwidth)*rep(2,nlevels)+1.5*barwidth,
                    y=cumsum(adipart.obj$statistic[indices]),yend=cumsum(adipart.obj$oecosimu$means[indices]),
                    labx=(1-barwidth)*1.5+1.5*barwidth)
  lines$laby<-(lines$y+lines$yend)/2
  lines$lab<-paste(signif(adipart.obj$oecosimu$z[indices],4),
                   symnum(adipart.obj$oecosimu$pval[indices], corr = FALSE, na = FALSE, 
                          cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                          symbols = c("***", "**", "*", ".", " ")),
                   sep="\n")
  plotdf<-data.frame(partition = factor(rep(names(adipart.obj$statistic)[indices],2),levels=rev(names(adipart.obj$statistic)[indices])),
                     source = factor(rep(c("Observed","Expected"),each=nlevels),levels=c("Observed","Expected")),
                     value = c(adipart.obj$statistic[indices],adipart.obj$oecosimu$means[indices]))
  ggplot(data=plotdf,aes(x=source,y=value))+
    geom_col(aes(fill=partition),width=barwidth)+
    geom_segment(data=lines,aes(x=x,xend=xend,y=y,yend=yend),linetype=3)+
    geom_label(data=lines,aes(x=labx,y=laby,label=lab),label.size=0,alpha=0.5)+
    labs(x="",y=paste0("OTU ",attr(adipart.obj$oecosimu$simulated,"index")),fill="")+
    scale_fill_discrete(labels=hierarchies)+
    theme(legend.position = "top")
}