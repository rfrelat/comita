#' Phase-randomised surrogate of a time series (base R)
#'
#' Keeps the amplitude spectrum (and hence the autocorrelation structure and
#' mean) of `x` and randomises its Fourier phases, as in
#' `fractal::surrogate(x, method = "phase")`.
#'
#' @param x numeric vector, at least 4 values, no NA
#' @return numeric vector of the same length as `x`
#' @keywords internal
phase_surrogate <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (anyNA(x)) {
    stop("x contains NA values", call. = FALSE)
  }
  if (n < 4) {
    stop("need at least 4 values", call. = FALSE)
  }

  f <- stats::fft(x)

  # Positive frequencies excluding DC (index 1) and, for even n, Nyquist.
  m <- (n - 1) %/% 2
  idx <- 2:(m + 1)

  # Random phases for the positive frequencies, original amplitudes kept.
  phi <- stats::runif(m, 0, 2 * pi)
  f[idx] <- Mod(f[idx]) * exp(1i * phi)

  # Enforce conjugate symmetry so the inverse transform is real.
  # Frequency k (1-based) mirrors to n - k + 2.
  f[n - idx + 2] <- Conj(f[idx])

  # Nyquist term (even n) must stay real; DC is left unchanged so the mean is kept.
  if (n %% 2 == 0) {
    f[n / 2 + 1] <- Re(f[n / 2 + 1])
  }

  # R's inverse fft is unnormalised, so divide by n.
  return(Re(stats::fft(f, inverse = TRUE)) / n)
}
