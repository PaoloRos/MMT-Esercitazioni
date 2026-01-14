# This file contains some functions for graphics formatting

ggbodeplot <- function(tf, fmin=1, fmax=1e4, df=0.01) {
  
  # vector of points for each order of magnitude (OOM):
  pts <- 10^seq(0, 1, df) %>% tail(-1)
  # vector of OOMs:
  ooms <- 10^(floor(log10(fmin)):ceiling(log10(fmax)-1))
  # combine pts and ooms:
  freqs <- as.vector(pts %o% ooms)
  
  # warning: bode wants pulsation!
  bode(tf, freqs*2*pi) %>% {
    tibble(
      f=.$w/(2*pi), 
      `Magnitude (dB)`=.$mag, 
      `Phase (deg)`=.$phase)} %>%
    pivot_longer(-f) %>% 
    ggplot(aes(x=f, y=value)) +
      geom_line() +
      scale_x_log10(
        minor_breaks=scales::minor_breaks_n(10), 
        labels= ~ latex2exp::TeX(paste0("$10^{", log10(.), "}$"))
        ) +
    facet_wrap(~name, nrow=2, scales="free") +
    labs(x="Frequency (Hz)")
}