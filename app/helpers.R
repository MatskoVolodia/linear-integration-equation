# public functions
k0      <- function(z) { -log(z/2 + C)*i0(z) + sigma0(z) }
i0      <- function(z) { sum( ( rep(1,N+1)/( factorial(0:N)^2 ) ) * ( (z^2)/4 ) ^ (0:N) ) }
sigma0  <- function(z) { sum( ( ksi(0:N)/( factorial(0:N)^2 ) ) * ( (z^2)/4 ) ^ (0:N) ) }
ksi     <- function(k) { unlist(lapply(k, function(k) { sum( ( rep(1,k)/(1:k) )) })) }

Rj      <- function(t, tj) { (-1/N) * (Rjtemp1(t,tj) + Rjtemp2(t,tj)) }

k1      <- function(z) { 1/z + (log(z/2)+C)*i1(z) + sigma1(z) }
i1      <- function(z) { sum( rep(1, N+1)/i1temp(0:N) * (z/2) ^ ((0:N)*2 + 1) ) }
sigma1  <- function(z) { (-1/2) * sum(ksi((0:N)+1) + ksi(0:N) / i1temp(0:N) * (z/2) ^ ((0:N)*2 + 1)) }
h       <- function(t, tau) { (htemp(t,tau,1,2) - htemp(t,tau,2,1)) / moddiff(t, tau)  }

H1prob  <- function(t,tau) { -kappa * k0(kappa*moddiff(t,tau)) }
H11     <- function(t,tau) { (-1/2) * i0(kappa(moddiff(t,tau))) }
H12     <- function(t,tau) { H1prob(t,tau) - H11(t,tau)*log(4*(sin((t-tau)/2)^2)) }
H1      <- function(t,tau) { ifelse(t==tau, xmod(t), H11(t,tau) + H12(t,tau)) }

H2prob  <- function(t,tau) { -kappa * k1(kappa*moddiff(t,tau)) * h(t,tau) }
H21     <- function(t,tau) { -kappa * i1(kappa(moddiff(t,tau))) * h(t,tau) }
H22     <- function(t,tau) { H2prob(t,tau) - H11(t,tau)*log(4*(sin((t-tau)/2)^2)) }
H2      <- function(t,tau) { ifelse(t==tau, H2lim(t), H21(t,tau) + H22(t,tau)) }
H2lim   <- function(t) { (derivX[[2]](t)*deriv2X[[1]](t) - derivX[[1]](t) * deriv2X[[2]](t)) / (xmod(t) ^ 2) }

# not for direct use
xmod    <- function(val) { sqrt(x[[1]](val)^2 + x[[2]](val)^2)}
Rjtemp1 <- function(t, tj) { sum( (rep(1,N-1) / (1:(N-1)) * cos((1:(N-1)) * (t-tj)) ) ) }
Rjtemp2 <- function(t, tj) { sum( (rep(1,N-1) / rep(1/2/N) * cos(N * (t-tj)) ) ) }
htemp   <- function(t, tau, i, j) { (x[[i]](tau) - x[[i]](t)) * derivX[[j]](tau) }
moddiff <- function(t, tau) { sqrt(difft(t,tau,1)^2+difft(t,tau,2)^2) }
difft   <- function(t,tau,i) { x[[i]](t) - x[[i]](tau) }
i1temp  <- function(n) { unlist(lapply(n, function(n) { factorial(n) * factorial(n+1) }))}