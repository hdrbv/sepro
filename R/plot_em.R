#' plot_em Function 
#'
#' This function created plot of mixture distributions as a result of separation process via EM algorithm (Expectation-Maximization algorithm) 
#'
#' @param x0 input data (vector)
#' @param EM application of EM function to input data
#' @return Plot of mixture 
#' @examples For example we want to separate 2 Gaussian distributions, 
#' @examples estimate parameters of each one and plot distributions.
#' @examples Let us assume that vector x1 - mixture of these distributions. 
#' @examples So we can use EM algorithm here:   
#' @examples EM1 <- sepro::EM(x0 = x0, k = 2).
#' @examples And create plot:
#' @examples plot_em(x0, EM1)
#' @author hdrbv
#' @export

plot_em <- function(x0, EM){
	plot_mix_comps <- function(x, mu, sigma, w) {
	  w * dnorm(x, mu, sigma)
	}	
	data.frame(x = x0) %>%
	ggplot() +
	geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
		   fill = "white") +
	stat_function(geom = "line", fun = plot_mix_comps,
		  args = list(EM$mu[1], sqrt(EM$var[1]), w = EM$w[1]),
		  colour = "red", lwd = 1.5) +
	stat_function(geom = "line", fun = plot_mix_comps,
		  args = list(EM1$mu[2], sqrt(EM$var[2]), w = EM$w[2]),
		  colour = "blue", lwd = 1.5) +
	ylab("Density") +
	xlab("Values") +
	ggtitle("Final plot of mixture")
}

