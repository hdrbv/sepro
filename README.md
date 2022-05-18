# SePRo: Separation Process with R

Estimating of parameters in statistical distributions using R. Package
provides an implementation of the EM algorithm (Expectation Maximization
Algorithm) in the R language. Today, EM and its variants are regularly
used to solve a broad range of today’s estimation problems, from the
multiple EM for motif elicitation (MEME) algorithm for motif-finding in
DNA sequences, to fitting mixture models to disambiguate targets from
clutter in radar. Hope that you, too, will find EM useful. You can see
details of usage
[here](https://github.com/hdrbv/sepro/blob/main/Tutorial/sepro-details.pdf/).

![*Example of mixture - 2 Gaussian functions*](graphics/ex1.png)

## Installation in R

**sepro** is a GitHub package so you can use ‘install\_github()‘ from
[devtools](https://cran.r-project.org/web/packages/devtools/index.html)
package.  Install **devtools** first:

    if("devtools" %in% rownames(installed.packages()) == FALSE){
    install.packages("devtools")
    }
    library(devtools)

 After that you can install **sepro** package:

    install_github("hdrbv/sepro", ref = "main")
    library(sepro)

## Theory

Let’s assume that we have some observed data
<span class="math inline">\\(y\\)</span>, a parametric density
<span class="math inline">\\(p(y|\\theta)\\)</span>, a description of
some complete data <span class="math inline">\\(x\\)</span> that we wish
we had, and the parametric density
<span class="math inline">\\(p(x|\\theta)\\)</span>. We assume that the
complete data can be modeled as a continuous random vector
<span class="math inline">\\(X\\)</span> with density
<span class="math inline">\\(p(x|\\theta)\\)</span>, where
<span class="math inline">\\(\\theta \\in \\Omega\\)</span> for some set
<span class="math inline">\\(\\Omega\\)</span>. We do not observe
<span class="math inline">\\(X\\)</span> directly, instead we observe a
realization <span class="math inline">\\(y\\)</span> of the random
vector <span class="math inline">\\(Y\\)</span> that depends on
<span class="math inline">\\(X\\)</span>. For example,
<span class="math inline">\\(Y\\)</span> might be the first component of
the vector <span class="math inline">\\(X\\)</span>.

Given that we only have <span class="math inline">\\(y\\)</span>, the
main goal here is to find the maximum likelihood estimate (MLE) of
<span class="math inline">\\(\\theta\\)</span>:

<span class="math display">\\\[\\hat{\\theta}\_{MLE} = argmax\\
p(y|\\theta)\\\]</span>

Is’s often easier to calculate the
<span class="math inline">\\(\\theta\\)</span> that maximizes the
log-likelihood of <span class="math inline">\\(y\\)</span>:

<span class="math display">\\\[\\hat{\\theta}\_{MLE} = argmax\\ log\\
p(y|\\theta)\\\]</span>

Because <span class="math inline">\\(log()\\)</span> is a monotonically
increasing function, solutions will be the same for both equations. But
sometimes it’s difficult to solve them. Then we can try EM: we make a
guess about the complete data <span class="math inline">\\(X\\)</span>
and solve for the <span class="math inline">\\(\\theta\\)</span> that
maximizes the (expected) log-likelihood of X. And once we have an
estimate for <span class="math inline">\\(\\theta\\)</span>, we can make
a better guess about the complete data
<span class="math inline">\\(X\\)</span>, and iterate.

Let’s break E-step and M-step of algorithm down into five steps:

1.  Let <span class="math inline">\\(m = 0\\)</span> and make initial
    estimate <span class="math inline">\\(\\theta^{(m)}\\)</span> for
    <span class="math inline">\\(\\theta\\)</span>

2.  Given the observed data <span class="math inline">\\(y\\)</span> and
    pretending for the moment that your current guess
    <span class="math inline">\\(\\theta^{(m)}\\)</span> is correct,
    formulate the conditional probability distribution
    <span class="math inline">\\(p(x|y, \\theta^{(m)}\\)</span> for the
    complete data <span class="math inline">\\(x\\)</span>

3.  Using the conditional probability distribution
    <span class="math inline">\\(p(x|y, \\theta^{(m)})\\)</span>
    calculated in step 2, form the conditional expected log-likelihood,
    which is called the Q-function:
    
    <span class="math display">\\\[Q(\\theta | \\theta^{m}) = \\int
    logp(x|\\theta)p(x|y, \\theta^{m})dx =\\\]</span>
    <span class="math display">\\\[= E\_{X|y,
    \\theta^{m}}(logp(X|\\theta))\\\]</span>
    
    where the integral is over set
    <span class="math inline">\\(\\chi(y)\\)</span>, which is the
    closure of the set <span class="math inline">\\({x|p(x|y, \\theta)\>
    0}\\)</span>, and assume that
    <span class="math inline">\\(\\chi(y)\\)</span> does not depend on
    <span class="math inline">\\(\\theta\\)</span>.

4.  Find the <span class="math inline">\\(\\theta\\)</span> that
    maximizes <span class="math inline">\\(Q\\)</span> - function;
    result = our new estimate =
    <span class="math inline">\\(\\theta^{(m+1)}\\)</span>

5.  Let m = m + 1 and go back to Step №2. EM algorithm does not specify
    a stopping criterion; standard criteria are to iterate until the
    estimate stops changing:
    <span class="math inline">\\(|\\theta^{(m+1)} - \\theta^{(m)}| \<
    \\epsilon\\)</span> for some
    <span class="math inline">\\(\\epsilon\\)</span> \> 0, or to iterate
    until the log-likelihood
    <span class="math inline">\\(l(\\theta^{m+1}) - l(\\theta^{m}) \<
    \\epsilon\\)</span> for some
    <span class="math inline">\\(\\epsilon\\)</span> \> 0

EM algorithm is only guaranteed to never get worse. Usually, it will
find a peak in the likelihood
<span class="math inline">\\(p(y|\\theta)\\)</span>, but if the
likelihood function <span class="math inline">\\(p(y|\\theta)\\)</span>
has multiple peaks, EM will not necessarily find the global maximum of
the likelihood. In practise, it’s common to start EM from multiple
random initial guesses, and choose the one with the largest likelihood
as the final guess for <span class="math inline">\\(\\theta\\)</span>

## Practice

Let’s apply theory to practise and also check the work of
[sepro](https://github.com/hdrbv/sepro) package. Firstly, let’s create
mixture of two distributions which we will separate:

    set.seed(1) #fix results of randomization
    cond <- sample(c(0, 1), size = 500, 
    replace = TRUE, prob = c(0.4, 0.6))
    # Sample from two different Gaussian distributions
    mix <- ifelse(cond == 1, rnorm(n = 500, mean = 5, sd = 1.5), 
    rnorm(n = 500, mean = 0, sd = 1))
    plot(mix)

![*Plot of created mixture (result of plot(mix))*](graphics/mix.png)

Apply <span class="math inline">\\(EM\\)</span> function from
[sepro](https://github.com/hdrbv/sepro) package to our mixture:

    vect <- as.numeric(mix)
    EM1 <- EM(vect, 2)

And use <span class="math inline">\\(plot\\\_em\\)</span> function from
package to see results of separation process:

    plot_em(vect, EM1)

![*Plot of separated mixture (result of plot\_em)*](graphics/plot_em.png)

That’s it. We have a fairly accurate parameter estimation of our
distributions - it’s really close to
real:

| :----------------------------------------------: | :--------------------------------------: | :----------------------------------------: |
|                    \[-1.8ex\]                    |              Expected value              |                  Variance                  |
|                                                  |                                          |                                            |
| \[-1.8ex\] Distribution 1 (from initial dataset) | <span class="math inline">\\(0\\)</span> |  <span class="math inline">\\(1\\)</span>  |
|      Distribution 2 (from initial dataset)       | <span class="math inline">\\(5\\)</span> | <span class="math inline">\\(1.5\\)</span> |
|    Distribution 1 (after separation process)     | <span class="math inline">\\(0\\)</span> |  <span class="math inline">\\(1\\)</span>  |
|    Distribution 2 (after separation process)     | <span class="math inline">\\(5\\)</span> |  <span class="math inline">\\(2\\)</span>  |

<span label=""></span>

