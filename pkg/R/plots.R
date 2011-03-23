#' @usage
#' plotLoadingsHeat <- function(RscApost, xlabel=NA, ylabel=NA, 
#' color=c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", 
#'      "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"), 
#' sorting=NA)
plotLoadingsHeat <- function(RscApost, xlabel=NA, ylabel=NA, 
                    color=c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", 
                          "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"), 
                    sorting=NA) {
  #ldf$Group = factor(ldf$Group, levels=rownames(scApl)[order(scApl[,1])])
	ldf = melt(RscApost)
	colnames(ldf) = c("Group", "Factor", "value")
  if (!any(is.na(sorting))) ldf$Group = factor(ldf$Group, levels=rownames(RscApost)[sorting])
  ldf$Factor = as.factor(ldf$Factor)
	breaks = seq(-1,1,by = 0.2)
	cl = colorRampPalette(color)(21)
	lim = c(-1.0, 1)
	p = ggplot(ldf, aes(x=Factor, y=Group, fill=value))
	p = p + scale_fill_gradientn(colour=cl, limits = lim, breaks=breaks) + geom_tile()
	if (!is.na(xlabel)) p = p+xlab(xlabel)
	if (!is.na(ylabel)) p = p+ylab(ylabel) 
  return(p)
}

plotLoadingsTxt <- function(RscApost) {
	ldf = melt(RscApost)
	ldf = cbind(ldf, as.numeric(ldf$value>0), as.factor(1:dim(RscApost)[1]))
	colnames(ldf) = c("Group", "Factor", "value", "sign", "y")
	ldf$Factor = as.factor(ldf$Factor)
	ldf$value = abs(ldf$value)
	p <- ggplot(data=ldf, aes(x=Factor, y=y, label=Group, color=sign, size=value)) 
	p <- p + geom_text(position="dodge") + xlab("") +scale_area(to = c(1, 5))
	print(p)
}