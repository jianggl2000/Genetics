#density plot
#data ploted by columns

myplotDensity = function(object, maintitle="", colors=NULL, xlab=NULL, ...){
    xrange = c(min(object),max(object))*1.1
    nARR = ncol(object)
    if(is.null("colors"))
        colors <- rainbow(nARR, s = 1, v = 1, start = 0, end = max(1, nARR - 1)/nARR)
    yrange=c(0,0)
    for(i in 1:nARR){
        if(yrange[2] < max(density(object[,i])$y))
            yrange[2] = max(density(object[,i])$y)
    }
    yrange[2] = yrange[2]*1.1
    for (n in 1:nARR) {
        if (n == 1) 
            plot(density(object[, n], na.rm = TRUE), col = colors[n], xlim=xrange, ylim = yrange, main = maintitle, xlab=xlab, ...)
        else lines(density(object[, n], na.rm = TRUE), col = colors[n], ...)
    }
}

cat("Please run myplotDensity(mydata)\n")