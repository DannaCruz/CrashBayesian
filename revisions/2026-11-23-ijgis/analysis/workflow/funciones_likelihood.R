# Vectorized, numerically stable NB2-log likelihood integrated against a normal iid effect.
log_nb_quadrature <- function(eta, phi, y, sigma, quadrature) {
  S<-length(eta);Q<-length(quadrature$nodes)
  log_ratio<-matrix(quadrature$nodes,S,Q,byrow=TRUE)*sigma+eta-log(phi)
  log1pexp<-pmax(log_ratio,0)+log1p(exp(-abs(log_ratio)))
  constants<-lgamma(y+phi)-lgamma(phi)-lgamma(y+1)
  vals<-constants+y*log_ratio-(y+phi)*log1pexp+
        matrix(log(quadrature$weights),S,Q,byrow=TRUE)
  mx<-matrixStats::rowMaxs(vals)
  mx+log(rowSums(exp(vals-mx)))
}
