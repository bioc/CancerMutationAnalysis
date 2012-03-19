ebfdr <- function(obs.score,
                  null.score,
                  mass,
                  estimate.p0=FALSE,
                  p0.step=1,
                  p0=1,
                  showFigure,
                  cutoffFdr,
                  score.name){
  
  obs.score <- obs.score[!is.na(obs.score)]
  null.score <- null.score[!is.na(null.score)]
  xxx <- sort(obs.score)
  nx <- length(xxx)
  obs.tail <- null.tail <- rep(0,nx)
  
  ##calculate null.tail and observed.tail
  calc.tail <- function(xxx.i, vect1, vect2)
    {
      return(c(sum(vect1 >= xxx.i), sum(vect2 >= xxx.i)))
    }
  tails <- sapply(xxx, calc.tail, null.score, obs.score)
  colnames(tails) <- NULL
  null.tail <- unlist(tails[1,])/mass
  obs.tail <- unlist(tails[2,])/1
  
  Fdr <- p0 * null.tail / obs.tail
  dens <- density(xxx)
  steps0 <- c(null.tail[1:(nx-1)]-null.tail[2:nx],0)
  if (sum(steps0) == 0 ) {
    fdr <- rep(0,nx)
  }
  else{
    dens0 <- density(xxx,weights=steps0/sum(steps0))
    fdr.grid <- p0 * ( max(null.tail) * ( dens0$y - min(dens0$y) ) ) /
      ( nx * ( dens$y - min(dens$y) + 1.0e-4) )
    ngrid <- length(dens$x)
    fff <- ff0 <- fdr <- rep(0,nx)
    
    ##calculate fff, ff0, and fdr at each observed score
    ##first get the closest point on the density grid
    get.iclose <- function(xx, dens.x, n)
      {
        max ( (1:n) [ abs(dens.x-xx) == min(abs(dens.x-xx)) ] )
      }
    iclose <- sapply(xxx, get.iclose, dens$x, ngrid)
    fff <-  dens$y[iclose]
    ff0 <-  dens0$y[iclose]
    fdr <- fdr.grid[iclose]
  }
  
  # monotonization (to the right of the max Fdr)
  
  maxFdrID <- max( (1:nx)[ Fdr==max(Fdr) ] )
  if (maxFdrID<nx) {
    for (i in (maxFdrID+1):nx ){
      Fdr[i] <- min(Fdr[i],Fdr[i-1])
      fdr[i] <- min(fdr[i],fdr[i-1])
    }
  }
  
  Fdr[Fdr>1] <- 1
  fdr[fdr>1] <- 1

  if(showFigure)
    {
      ##make a smoother density plot
      dens0 <- density(xxx,weights=steps0/sum(steps0), bw=0.2)
      ##get the largest Fdr smaller than cutoffFdr
      cutoffFdr.index <- sum(sort(Fdr) < cutoffFdr)
      cutoffFdr <- sort(Fdr)[cutoffFdr.index]
      cutoffFdr.index <- which(Fdr == cutoffFdr)
      ##now get the observed score corresponding to this cutoff
      cutoff.obs <- min(sort(obs.score)[cutoffFdr.index])
      ##how many real genes have scores above this cutoff?
      genes.above <- sum(obs.score >= cutoff.obs)
      ##how many null genes on average have scores above this cutoff?
      null.genes.above <- round(sum(null.score >= cutoff.obs)/mass, 2)
      ##get cutoff for double that number of genes
      xmin <- sort(obs.score, decreasing=TRUE)[2*genes.above]
      ##xmax <- sort(obs.score, decreasing=TRUE)[min(3, genes.above)]
      xmax <- round(min(cutoff.obs+3, sort(obs.score, decreasing=TRUE)[2]))+1
      ##get ymax
      ymax <- dens0$y[sum(dens0$x <= xmin)]
      ##turn off warnings (to not have warning for rug values being clipped)
      options(warn=-1)
      plot(dens0, xlim = c(xmin, xmax), ylim = c(0, ymax), lwd = 2,
           xlab="score", type="l",
           ylab="Density of null scores",
           main=paste(score.name, "score", sep=" "))
      xlabText <- max(0.25*(xmax+xmin), round(xmin)+0.5)
      text(c(xlabText,xlabText),
           c(ymax*0.55,ymax*0.45),
           labels=c(paste("Real genes to the right of cutoff:",
             genes.above,sep=" "),
             paste("Null genes to the right of cutoff:",
                   null.genes.above,sep=" ")),pos=4,
           col=c("blue","black"), cex=1.2)
      abline(v=cutoff.obs)
      rug(obs.score, col="blue")
      options(warn=0)
    }
    
  # estimate p0

  if (estimate.p0) {
    scores <- c(obs.score, null.score)
    maximum <- ceiling(max(scores))
    minimum <- floor(min(scores))
    
    hh <- hist(obs.score,breaks=seq(minimum,maximum,p0.step),
               plot=FALSE)
    h0 <- hist(null.score,breaks=seq(minimum,maximum,p0.step),
               plot=FALSE)
    
    nonzero <- hh$density > 0

    p0 = 1 / max(h0$density[nonzero]/hh$density[nonzero], na.rm = TRUE)

  }
  
  out <- as.data.frame(cbind(xxx,obs.tail,null.tail,Fdr,fdr,p0))
  rownames(out) <- names(xxx)
  colnames(out) <- c("Score","F","F0","Fdr","fdr","p0")

  return(out)

}
