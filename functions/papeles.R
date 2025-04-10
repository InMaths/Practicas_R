# papeles.R
#
# Copyright (C) 2011  Jesus Fernandez (fernandej IN THE DOMAIN unican DOT es)
# Universidad de Cantabria, Spain
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# Ultima modificacion:
#   2017-11-03 Jesus Fernandez
#     Incluyo papel de minimos y nuevo argumento ppos para elegir punteo
#   2021-05-06 Maria Dolores Frias
#     Incluyo la funcion papel.minimal.lines
#     papel.exponencial devuelve los parametros en una lista
#     papel.maximal y papel.minimal toman ahora por defecto Hazen
#library(Hmisc) # minor.tick
if (! exists("minor.tick")) {minor.tick <- function(...){TRUE}}

eta.norm <- function(p){
  return(qnorm(p))
}

eta.exp <- function(p){
  return(-log(1-p))
}

eta.max <- function(p) -log(-log(p))
eta.max.inv <- function(x) exp(-exp(-x))

eta.min <- function(p) log(-log(1-p))
eta.min.inv <- function(x) 1-exp(-exp(x))

xi.exp <- function(x){
  return(x)
}

xi.gumbel.max <- function(x, lamb="unused!"){
  return(x)
}

xi.weibull.max <- function(x, lamb){
  return(-log(lamb-x))
}

xi.frechet.max <- function(x, lamb){
  return(log(x-lamb))
}

xi.gumbel.min <- function(x, lamb="unused!") x
xi.weibull.min <- function(x, lamb) log(x-lamb)
xi.frechet.min <- function(x, lamb) -log(lamb-x)

punteo.basic <- function(n){
  return((1:n)/(n+1))
}
punteo.basic.inv <- function(p,n) p*(n+1)

punteo.hazen <- function(n) (1:n-0.5)/n
punteo.hazen.inv <- function(p,n) n*p+0.5

punteo.blom <- function(n){
  return((1:n-3/8)/(n+0.25))
}

punteo.gringorten <- function(n){
  return((1:n-0.44)/(n+0.12))
}

std.pvals<-c(seq(0.001,0.01,0.001), seq(0.01,0.1,0.01), seq(0.1,0.9,0.1), seq(0.91,0.99,0.01))

papel.normal <- function(xs, trim=0, show.fit=F, x.transform="none", ...){
  x.i <- sort(xs)
  p.i <- punteo.hazen(length(x.i))
  plevels <- std.pvals
  labellevels <- std.pvals[-c(1,2,3,4,6,7,8,9,14,16,17,18,19)]
  xi.i <- x.i
  par(mar=c(4,4,2,4))
  if (x.transform == "none"){
    xlab <- substitute(xi==x)
  }
  else if (x.transform == "log"){
    # Opcion no implementada
    xlab <- substitute(xi==ln(x))
  }
  ylab <- substitute(eta)
  plot(xi.i, eta.norm(p.i), xlab=xlab, ylab=ylab, ...)
  minor.tick(nx=10, tick.ratio=0.5)
  minor.tick(nx=5, tick.ratio=1)
  abline(h=eta.norm(plevels)) #, col="gray")
  mtext(sprintf("%5.3f", labellevels), 4, 0.2,
    at=eta.norm(labellevels),
    las=1, cex=0.6)
  if (show.fit){
    abline(h=0,lwd=3)
    abline(h=1,lty=2)
    these <- (1+trim):(length(p.i)-trim)
    points(xi.i[these],eta.norm(p.i[these]), pch=19)
    mylm <- lm(eta.norm(p.i[these])~xi.i[these])
    abline(mylm, lwd=2, col="red")
    a<-mylm$coefficients[1]
    b<-mylm$coefficients[2]
    eta0 <- -a/b
    abline(v=eta0); mtext(round(eta0,2), 3, 0.7, at=eta0, cex=0.8)
    eta1 <- (1-a)/b
    abline(v=eta1); mtext(round(eta1,2), 3, 0.7, at=eta1, cex=0.8)
    if (x.transform == "none"){
      return(list(
        mean=as.numeric(eta0), sd=as.numeric(eta1-eta0)
      ))
    }
  }
}

papel.exponencial <- function(xs, trim=0, show.fit=F, allow.shift=TRUE, ...){
  x.i <- sort(xs)
  p.i <- punteo.hazen(length(x.i))
  plevels <- std.pvals
  labellevels <- std.pvals[-c(1,2,3,4,6,7,8,9,14,16,17,18,19)]
  xi.i <- xi.exp(x.i)
  par(mar=c(4,4,2,4))
  xlab <- substitute(xi==x)
  ylab <- substitute(eta)
  plot(xi.i, eta.exp(p.i), xlab=xlab, ylab=ylab, ...)
  minor.tick(nx=10, tick.ratio=0.5)
  minor.tick(nx=5, tick.ratio=1)
  abline(h=eta.exp(plevels)) #, col="gray")
  mtext(sprintf("%5.3f", labellevels), 4, 0.2,
    at=eta.exp(labellevels),
    las=1, cex=0.6)
  if (show.fit){
    abline(h=0,lwd=3)
    abline(h=1,lty=2)
    these <- 1:(length(p.i)-trim)
    points(xi.i[these],eta.exp(p.i[these]), pch=19)
    if (allow.shift){
      mylm <- lm(eta.exp(p.i[these])~xi.i[these])
      a<-mylm$coefficients[1]
      b<-mylm$coefficients[2]
    }
    else {
      mylm <- lm(eta.exp(p.i[these])~ -1+xi.i[these])
      a<-0
      b<-mylm$coefficients[1]
    }
    abline(mylm, lwd=2, col="red")
    eta0 <- -a/b
    abline(v=eta0); mtext(round(eta0,2), 3, 0.7, at=eta0, cex=0.8)
    eta1 <- (1-a)/b
    abline(v=eta1); mtext(round(eta1,2), 3, 0.7, at=eta1, cex=0.8)
    return(list(shift=as.numeric(eta0), rate=as.numeric(1/(eta1-eta0))))
  }
}

papel.maximal <- function(xs, fit.only=0, show.fit=F, x.transform="gumbel", lambda=0, ppos="hazen", ...){
  n <- length(xs)
  x.i <- sort(xs)
  p.i <- do.call(sprintf("punteo.%s", ppos), list(n))
  plevels <- std.pvals
  labellevels <- std.pvals[-c(1,2,3,4,6,7,8,9,14,16,17,18,19)]
  xi.i <- do.call(paste("xi.",x.transform,".max", sep=''), list(x=x.i, lamb=lambda))
  par(mar=c(4,4,2,4))
  if (x.transform == "gumbel"){
    xlab <- substitute(xi==x)
  }
  else if (x.transform == "weibull"){
    lamb <- round(lambda,2)
    xlab <- substitute(xi==-ln(lamb-x))
  }
  else if (x.transform == "frechet"){
    lamb <- round(lambda,2)
    xlab <- substitute(xi==ln(x-lamb))
  }
  ylab <- substitute(eta)
  plot(xi.i, eta.max(p.i), xlab=xlab, ylab=ylab, ...)
  minor.tick(nx=10, tick.ratio=0.5)
  minor.tick(nx=5, tick.ratio=1)
  abline(h=eta.max(plevels)) #, col="gray")
  mtext(sprintf("%5.3f", labellevels), 4, 0.2,
    at=eta.max(labellevels),
    las=1, cex=0.6)
  mtext(sprintf("%6.2f",(1/(1-labellevels))), 4, 3.5,
    at=eta.max(labellevels),
    las=1, cex=0.6, adj=1)
  if (show.fit){
    abline(h=0,lwd=3)
    abline(h=1,lty=2)
    these <- 1:n
    if (fit.only != 0){
      cutfrom <- n-fit.only
      these <- cutfrom:n
      points(xi.i[these],eta.max(p.i[these]), pch=19)
    }
    mylm <- lm(eta.max(p.i[these])~xi.i[these])
    abline(mylm, lwd=2, col="red")
    a<-mylm$coefficients[1]
    b<-mylm$coefficients[2]
    # Add confidence bounds
    xx <- seq(min(xi.i[these]), max(xi.i[these]), length.out = 100)
    rr <- do.call(sprintf("punteo.%s.inv", ppos), list(p=eta.max.inv(a+b*xx), n=n)) # Invert basic plotting pos to get order (r)
    lines(xx, eta.max(qbeta(0.975, rr, n-rr+1)), lty="dotted", lwd=3, col="red")
    lines(xx, eta.max(qbeta(0.025, rr, n-rr+1)), lty="dotted", lwd=3, col="red")
    # Find the parameters
    eta0 <- -a/b
    abline(v=eta0); mtext(round(eta0,2), 3, 0.7, at=eta0, cex=0.8)
    eta1 <- (1-a)/b
    abline(v=eta1); mtext(round(eta1,2), 3, 0.7, at=eta1, cex=0.8)
    if (x.transform == "gumbel"){
      return(list(loc=as.numeric(eta0), scale=as.numeric(eta1-eta0)))
    }
    else if (x.transform == "weibull"){
      return(list(loc=as.numeric(lambda), scale=as.numeric(exp(-eta0)), shape=as.numeric(1/(eta1-eta0))))
    }
    else if (x.transform == "frechet"){
      return(list(loc=as.numeric(lambda), scale=as.numeric(exp(eta0)), shape=as.numeric(1/(eta1-eta0))))
    }
  }
}

papel.maximal.lines <- function(xs, ppos="hazen", ...){
  x.i <- sort(xs)
  p.i <- do.call(sprintf("punteo.%s", ppos), list(length(x.i)))
  lines(x.i, eta.max(p.i), ...)
}

papel.minimal <- function(xs, fit.only=0, show.fit=F, x.transform="gumbel", lambda=0, ppos="hazen", ...){
  n <- length(xs)
  x.i <- sort(xs)
  p.i <- do.call(sprintf("punteo.%s", ppos), list(n))
  plevels <- std.pvals
  labellevels <- std.pvals[-c(1,2,3,4,6,7,8,9,14,16,17,18,19)]
  xi.i <- do.call(paste("xi.",x.transform,".min", sep=''), list(x=x.i, lamb=lambda))
  par(mar=c(4,4,2,4))
  if (x.transform == "gumbel"){
    xlab <- substitute(xi==x)
  }
  else if (x.transform == "weibull"){
    lamb <- round(lambda,2)
    xlab <- substitute(xi==ln(x-lamb))
  }
  else if (x.transform == "frechet"){
    lamb <- round(lambda,2)
    xlab <- substitute(xi==-ln(lamb-x))
  }
  ylab <- substitute(eta)
  plot(xi.i, eta.min(p.i), xlab=xlab, ylab=ylab, ...)
  minor.tick(nx=10, tick.ratio=0.5)
  minor.tick(nx=5, tick.ratio=1)
  abline(h=eta.min(plevels)) #, col="gray")
  mtext(sprintf("%5.3f", labellevels), 4, 0.2,
        at=eta.min(labellevels),
        las=1, cex=0.6)
  mtext(sprintf("%6.2f",(1/labellevels)), 4, 3.5,
        at=eta.min(labellevels),
        las=1, cex=0.6, adj=1)
  if (show.fit){
    abline(h=0,lwd=3)
    abline(h=1,lty=2)
    these <- 1:n
    if (fit.only != 0){
      these <- 1:fit.only
      points(xi.i[these],eta.min(p.i[these]), pch=19)
    }
    mylm <- lm(eta.min(p.i[these])~xi.i[these])
    abline(mylm, lwd=2, col="red")
    a<-mylm$coefficients[1]
    b<-mylm$coefficients[2]
    # Add confidence bounds
    xx <- seq(min(xi.i[these]), max(xi.i[these]), length.out = 100)
    rr <- do.call(sprintf("punteo.%s.inv", ppos), list(p=eta.min.inv(a+b*xx), n=n)) # Invert basic plotting pos to get order (r)
    lines(xx, eta.min(qbeta(0.975, rr, n-rr+1)), lty="dotted", lwd=3, col="red")
    lines(xx, eta.min(qbeta(0.025, rr, n-rr+1)), lty="dotted", lwd=3, col="red")
    # Find the parameters
    eta0 <- -a/b
    abline(v=eta0); mtext(round(eta0,2), 3, 0.7, at=eta0, cex=0.8)
    eta1 <- (1-a)/b
    abline(v=eta1); mtext(round(eta1,2), 3, 0.7, at=eta1, cex=0.8)
    if (x.transform == "gumbel"){
      print(sprintf("loc: %f scale: %f", eta0, eta1-eta0))
    }
    else if (x.transform == "weibull"){
      print(sprintf("loc: %f scale: %f shape: %f", lambda, exp(eta0), 1/(eta1-eta0)))
    }
    else if (x.transform == "frechet"){
      print(sprintf("loc: %f scale: %f shape: %f", lambda, exp(-eta0), -1/(eta1-eta0)))
    }
  }
}

papel.minimal.lines <- function(xs, ppos="hazen", ...){
  x.i <- sort(xs)
  p.i <- do.call(sprintf("punteo.%s", ppos), list(length(x.i)))
  lines(x.i, eta.min(p.i), ...)
}
