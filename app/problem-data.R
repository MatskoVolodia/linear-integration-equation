x = c(
  function(t) { 3*cos(t) },
  function(t) { 3*sin(t) }
)

derivX = c(
  Deriv(x[[1]], 't'),
  Deriv(x[[2]], 't')
)

deriv2X = c(
  Deriv(Deriv(x[[1]], 't'), 't'),
  Deriv(Deriv(x[[2]], 't'), 't')
)