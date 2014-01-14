#### internals ####
no.dimnames <- function(a){
  d <- list()
  l <- 0
  for (i in dim(a)){
    d[[l<-l+1]]<-rep("",i)
  }
  dimnames(a)<-d
  a
}


#### Benford ####

DF <- function(data){
  collapsed <- 10*(data)/10^trunc(log10(data))
  DF<-(mean(collapsed)-39.0865)/(39.0685)
  DF*100
}

mantissa.arc.test <- function(data){ #data must be the mantissa of the log10
  x.coord <- cos(2*pi*data)
  y.coord <- sin(2*pi*data)
  L2 <- (mean(x.coord))^2 + (mean(y.coord))^2
  names(L2) <- "L2"
  p.value <- exp(-(L2)*length(data))
  return(list(L2=L2, p.value=p.value))
}

extract.digits <- function (data, number.of.digits = 2, sign="positive", second.order = FALSE, discrete=TRUE, round=3) {
  
  if(class(data)!="numeric"){stop("Data must be a numeric vector")}
  
  ## cleaning data for analysis - only > 0 and either only positive or only negative
  if (sign == "positive") {positives <- data[data>0 & !is.na(data)]} 
  if (sign == "negative") {positives <- data[data<0 & !is.na(data)]*(-1)}
  if (sign == "both")     {positives <- abs(data[data!=0 & !is.na(data)])}
  if (second.order){
    if (number.of.digits>4){warning("There might be some floating point precision issues on the Second Order distribution")}
    n <- length(positives)
    first <- sort(positives)[1:(n-1)]
    second <- sort(positives)[2:n]
    positives <-  if(discrete){round(second - first, number.of.digits+round)}else{second - first} 
    positives <- positives[positives>0]
  }
  return(data.frame(
    data = positives,
    data.digits = trunc((10^((floor(log10(positives))*-1)+number.of.digits-1))*positives)))
}

p.these.digits <- function(d){
  
  if (class(d)!="numeric" & class(d)!="integer") stop ("d must be numeric or integer")
  d <- trunc(d)
  if (d<0) d <- d*(-1)
  prob <- log10(1+1/d)
  return(prob)
}

p.this.digit.at.n <- function(d,n){
  if (d<0) d <- d*(-1)
  n1 <- strsplit(as.character(d), "")[[1]]
  n1 <- length(n1)
  if (n1>1) stop("d must have only 1 digit. This function evaluates 1 digit at position n")
  if (class(d)!="numeric") stop ("d must be numeric")
  if (class(n)!="numeric") stop ("d must be numeric")
  if (n<1) stop ("n must be greater than 1")
  if (n==1) return(log10(1+1/d))
  sum = 0
  k <- 10^(n-2)
  j <- (10^(n-1))-1
  soma = 0
  for (i in k:j){
    soma <- soma + log10(1 + 1/(10*i + d))
  }
  soma
}

generate.benford.digits <- function (number.of.digits) {
  number.of.digits <- as.integer(number.of.digits)
  begin <- 10^(number.of.digits-1); ends <- 10^(number.of.digits)-1;
  benford.digits <- begin:ends
  return(benford.digits)
}

generate.benford.distribution <- function (benford.digits) {
  benford.dist <- sapply(benford.digits, p.these.digits)
  return(benford.dist)
}

generate.empirical.distribution <- function (data, number.of.digits,sign, second.order=FALSE, benford.digits, discrete=TRUE, round=3){
   x <- NULL
   v <- NULL
   data.frame <- extract.digits(data, number.of.digits, sign, second.order, discrete=discrete, round=round)
   n <- length(data.frame$data.digits)
#    dist.freq <- table(c(data.frame$data.digits, benford.digits))-1
#   require(data.table)
   DF<- data.table(x=c(data.frame$data.digits, benford.digits),
                   v=c(data.frame$data.digits, benford.digits) )
   DFcount <- DF[,length(x)-1, by=v][order(v)]
   dist.freq <- DFcount$V1
   dist <- dist.freq/n
   return(list(
     data= data.frame$data, 
     data.digits= data.frame$data.digits, 
     dist = dist, 
     dist.freq=dist.freq))
}

extract.mantissa <- function (positives) {
  log <- log10(positives)
  log[log<0] <- log[log<0]+as.integer(log[log<0])*(-1)+1
  mantissa <- log - trunc(log)
}

generate.summation <- function (benford.digits, data, data.digits) {
x<- NULL
v <- NULL
#   table <-aggregate(data,by=list(data.digits), sum) 
#   names(table)<- c("digits", "value")
#   require(data.table)
  table <- data.table(x=data.digits, v=data)
  table <- table[, sum(v), by=x][order(x)]
  setnames(table,c("x", "V1"), c("digits", "value"))
  
  if(length(which(!benford.digits %in% table$digits))!=0){
  add <- data.frame(digits=which(!benford.digits %in% table$digits), value=0)
  add
  table<-rbind(table, add)
  table<-table[order(table$digits),]
  }
  
  summation <- table$value
  return(summation)
}



#### Basic Calculations ####

excess.kurtosis <- function (x) {
  
  
  (mean((x - mean(x))^4)/(mean((x - mean(x))^2)^2))-3}

skewness <- function (x) {
  
  
  (mean((x - mean(x))^3)/(mean((x - mean(x))^2)^(3/2)))}

#### plot ####

plotting.data.vs.benford <- function (x, ...) {
  
  xmarks <- barplot(x[["bfd"]]$data.dist.freq, col="lightblue", 
                    main= "Digits Distribution",
                    xlab = "Digits", ylab="Freq",
                    ylim = c(0,max(c(x[["bfd"]]$data.dist.freq, x[["bfd"]]$benford.dist.freq))*1.1))
  
  axis(1, at=xmarks, labels=x[["bfd"]]$digits)
  
  lines(xmarks, x[["bfd"]]$benford.dist.freq, lty=2, lwd=2, col="red")
  
#   legend("topright", lty=c(1,2), lwd=2, col=c("lightblue", "red"), 
#          legend=c("Data", "Benford"))
}

plotting.second.order <- function (x, ...) {
  y<-x[["bfd"]]$benford.so.dist.freq
  xmarks <- barplot(x[["bfd"]]$data.second.order.dist.freq, col="lightblue", 
                    main= "Digits Distribution \nSecond Order Test",
                    xlab = "Digits", ylab="Freq",
                    ylim = c(0,max(c(x[["bfd"]]$data.second.order.dist.freq, y))*1.1)) 
  
  axis(1, at=xmarks, labels=x[["bfd"]]$digits)
  
  lines(xmarks, y, lty=2, lwd=2, col="red")
  
#   legend("topright", lty=c(1,2), lwd=2, col=c("lightblue", "red"), 
#          legend=c("Data", "Benford"))
}

plotting.summation <- function (x, ...) {
  xmarks<-barplot(x[["bfd"]]$data.summation, col="lightblue", 
                  main="Summation Distribution by digits",
                  xlab="Digits", ylab="Summation",
                  ylim = c(0,max(x[["bfd"]]$data.summation))*1.1)
  
  axis(1, at=xmarks, labels=x[["bfd"]]$digits)
  
  lines(x=xmarks, y=rep(mean(x[["bfd"]]$data.summation), length(xmarks)), col="red", lty=2)
}

plotting.ordered.mantissa <- function (x, ...) {
  
  plot(sort(x[["data"]]$data.mantissa), pch=".",col="blue", main ="Ordered Mantissa",
       xlab="Ordered Observation",
       ylab= "Mantissa")
  
  abline(a=0, b=1/length(x[["data"]]$data.mantissa), col="red", lty=2)
}

plotting.chi_squared <- function (x, ...) {
  
  plot(x[["bfd"]]$digits, x[["bfd"]]$squared.diff, pch="x", col="blue", 
       xlab="Digits",
       ylab="Chi-squared", main= "Chi-Squared Difference",xaxt="n")
  
  axis(1, at= x[["bfd"]]$digits)
}

plotting.abs.diff <- function (x, ...) {
  plot(x[["bfd"]]$digits, x[["bfd"]]$absolute.diff, pch="x", col="blue", xlab="Digits",
       ylab="Absolute Difference", 
       main= "Absolute Difference",xaxt="n")
  
  axis(1, at= x[["bfd"]]$digits)
}

plotting.ex.summation <- function (x, ...) {
  plot(x[["bfd"]]$digits, x[["bfd"]]$abs.excess.summation, pch="x", col="blue", 
       xlab="Digits",
       ylab="Absolute Excess Summation", 
       main= "Summation Difference",xaxt="n")
  axis(1, at= x[["bfd"]]$digits)
}

plotting.legend <- function (x) {
  plot(1, type = "n", axes=FALSE, xlab="", ylab="", main= paste("Legend \n Dataset:", x[["info"]]$data.name))
  plot_colors <- c("lightblue","blue","red")
  legend(x = "top",inset = 0,
         legend = c("data", "data", "benford"), 
         col=plot_colors, lwd=2, lty=c(1,1,2))
}
