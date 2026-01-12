# This file contains some main function related to the digital signal processing

# Scalar product
`%ps%` <- function(a, b) mean(a*b)

# Norm of a vector
norm_ps <- function(a) sqrt(a %ps% a)

# Harmonic components

cos_k <- function(t, f, k=1) cos(2*pi*f*k*t)
sin_k <- function(t, f, k=1) sin(2*pi*f*k*t)