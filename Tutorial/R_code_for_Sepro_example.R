################################
### SePRo - Example of usage ###
################################

#Import additional packages:
if("pacman" %in% rownames(installed.packages()) == FALSE){
	install.packages("pacman")}
pacman:: p_load(dplyr, ggplot2, stargazer)

install_github("hdrbv/sepro", ref = "main", force = T)
library(sepro)


#Firstly, let's create mixture of two distributions which we will separate
set.seed(1) #fix results of randomization
cond <- sample(c(0, 1), size = 500, 
						   replace = TRUE, prob = c(0.4, 0.6))
# Sample from two different Gaussian distributions
mix <- ifelse(cond == 1, rnorm(n = 500, mean = 5, sd = 1.5), 
					    rnorm(n = 500, mean = 0, sd = 1))
plot(mix)

#Apply EM function from sepro package to our mixture
vect <- as.numeric(mix)
EM1 <- EM(vect, 2)

#And use $ plot\_em $ function from package to see results of separation process
plot_em(vect, EM1)

#Let's create final table
in_table <- data.frame(c(0, 5), c(1, 1.5))
colnames(in_table) <- c("Expected value", "Variance")
btable <- data.frame(c(EM1$mu[1], EM1$mu[2]), c(EM1$var[1], EM1$var[2])) %>% round(0)
colnames(btable) <- c("Expected value", "Variance")
btable <- rbind(in_table, btable)
rownames(btable) <- c("Distribution 1 (from initial dataset)", 
		      "Distribution 2 (from initial dataset)",
		      "Distribution 1 (after separation process)", 
		      "Distribution 2 (after separation process)")
colnames(btable) <- c("Expected value", "Variance")
btable <- as.data.frame(btable)
stargazer(btable, type = "latex", summary = F)


