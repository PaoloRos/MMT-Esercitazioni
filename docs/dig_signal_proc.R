# This file contains some main function related to the digital signal processing

# --- Signal operations ---

# Scalar product
`%ps%` <- function(a, b) mean(a*b)

# Norm of a vector
norm_ps <- function(a) sqrt(a %ps% a)

# Scalar normalized product between n decimated vectors: i.e. cosine similarity
psDnN <- function(a, b, n = 1) {
  stopifnot(n >= 1, n %% 1 == 0)
  
  idx <-  seq(1, min(length(a), length(b)), by = n)
  
  a_d <- a[idx]
  b_d <- b[idx]
  
  mean(a_d * b_d) / (norm_ps(a_d)*norm_ps(b_d))
}

# --- Fourier Transform ---

# Harmonic components
cos_k <- function(t, f, k=1) cos(2*pi*f*k*t)
sin_k <- function(t, f, k=1) sin(2*pi*f*k*t)

compute_fft <- function(df, Ta, name_sig, monolateral = TRUE) {
  
  # Compute the FFT of a signal
  # df: data frame
  # Ta: acquisition time
  # name_sig: name of the signal contained in `df`
  
  N <- nrow(df)
  
  df <- df %>% mutate(
    n = row_number(),
    f = n/Ta,
    fft = fft({{ name_sig }}),
    mod = ifelse(monolateral, (2-(f==0)), 1) * Mod(fft) / length(n),
    phase = Arg(fft)/pi*180
  )
  
  if(monolateral)
    return(df %>% slice_head(n = floor(N/2)))
  
  df
}

# --- Filters ---

# Gaussian mask: law-pass filter, FIR and linear
gaussian_kernel <- function(
    x, sigma, radius = ceiling(3 * sigma), pad = c("none", "replicate", "reflect")
    ) {
  
  stopifnot(is.numeric(x), length(x) > 0, is.numeric(sigma), sigma > 0)
  pad <- match.arg(pad)
  
  # Kernels calculation
  k <- seq(-radius, radius)
  h <- exp(-(k^2) / (2 * sigma^2))
  h <- h / sum(h)  
  
  if (pad == "none") {
    y <- stats::filter(x, h, sides = 2) # convolution
    return(as.numeric(y))
  }
  
  # Padding of signal tails
  if (pad == "replicate") {
    x_pad <- c(rep(x[1], radius), x, rep(x[length(x)], radius))
  } else { # reflect
    
    left  <- rev(x[2:(radius + 1)])
    right <- rev(x[(length(x) - radius):(length(x) - 1)])
    x_pad <- c(left, x, right)
  }
  
  # Convolution and truncation of the excess tails
  y_pad <- stats::filter(x_pad, h, sides = 2)
  y <- y_pad[(radius + 1):(radius + length(x))]
  
  as.numeric(y)
}
