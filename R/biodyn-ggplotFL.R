#' plot(biodyn)
#'
#' Basic plot for biodyn.
#'
#' @aliases plot,biodyn,missing-method
#' @docType methods
#' @rdname plot
#' @examples
#' \dontrun{
#'   bd=simBiodyn()
#'   bd=fwd(bd,harvest=rlnorm(100,harvest(bd)[,-1]))
#'   plot(bd)
#'   }

setMethod("plot", signature(x="biodyn", y="missing"),
      function(x, xlab="Year", ylab="", 
               fn=list(Stock  =stock,
                       Catch  =catch,
                       Harvest=harvest),
               probs=c(0.25,0.5,0.75),na.rm=TRUE,type=7,
               rel=FALSE,...) {
            
          # extract info to plot
          fqs = FLQuants(mlply(fn, function(fn) fn(x)))
          
          p <- plot(fqs,na.rm=na.rm,probs=probs,type=type)+
            xlab(xlab)+ylab(ylab)
            
          return(p)})

#' plot(biodyns) 
#' 
#' plot() method for biodyns
#'
#' 
#' @aliases plot,biodyns,missing-method
#' @rdname plot
#' @examples
#' \dontrun{
#'   data(bd)
#'   bds <- biodyns(runA=bd, runB=bd)
#'   plot(bds)
#'   }

setMethod("plot", signature(x="biodyns", y="missing"),
          function(x, xlab="Year", ylab="", 
                   fn=list(Stock  =stock,
                           Catch  =catch,
                           Harvest=harvest),
                   probs=c(0.25,0.5,0.75),na.rm=TRUE,type=7,
                   rel=FALSE,...) {
            
            cos <- c('red', 'blue')
            
            # extract slots by stock
            fqs <- lapply(x, function(y) FLQuants(Stock=stock(y),
                                                  Catch=catch(y), Harvest=fbar(y)))
            
            # get median & 85% quantiles if iters
            its <- unlist(lapply(x, function(x) dims(x)$iter))
            if(any(its > 1))
            {
              # quantiles
              fqs <- lapply(fqs, function(y) as.data.frame(lapply(y, quantile,
                                                                  c(0.10, 0.50, 0.90))))
            } else {
              fqs <- lapply(fqs, as.data.frame)
              fqs <- lapply(fqs, function(x) {x$iter <- "50%"; return(x)})
            }
            
            # stock names
            stk <- rep.int(names(fqs), unlist(lapply(fqs, nrow)))
            # rbind dfs
            fqs <- do.call(rbind, fqs)
            rownames(fqs) <- NULL
            # add stock names
            fqs <- transform(fqs, stock=stk)
            
            # cast with quantiles in columns
            df <- dcast(fqs, age+year+unit+season+area+qname+stock~iter, value.var="data")
            
            # plot data vs. year + facet on qname +
            p <- ggplot(data=df, aes(x=year, y=`50%`, group=stock)) +
              facet_grid(qname~., scales="free") +
              # line + xlab + ylab + limits to include 0 +
              geom_line(aes(colour=stock)) + xlab(xlab) + ylab(ylab) + expand_limits(y=0) +
              # no legend
              theme(legend.title = element_blank())
            
            # object w/ iters?
            if(any(unlist(lapply(x, function(y) dims(y)$iter)) > 1)) {
              p <- p +
                # 75% quantile ribbon in red, alpha=0.25
                geom_ribbon(aes(x=year, ymin = `10%`, ymax = `90%`, group=stock,
                                colour=stock, fill=stock), alpha = .20, linetype = 0)
              # 90% quantile ribbon in red, aplha=0.10
            }
            return(p)}) 
