Fits the lifespan datasets of biological systems such as yeast, fruit flies, and other similar biological units with well-known finite mixture models introduced by Farewell et al. (1982) <doi:10.2307/2529885> and Al-Hussaini et al. (2000) <doi:10.1080/00949650008812033>. Estimates parameter space fitting of a lifespan dataset with finite mixtures of parametric distributions. Computes the following tasks; 1) Estimates parameter space of the finite mixture model
by implementing the expectation maximization (EM) algorithm. 2) Finds a sequence of four goodness-of-fit measures consist of Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), Kolmogorov-Smirnov (KS), and log-likelihood (log-likelihood) statistics. 3)The initial values is determined by k-means clustering.

Mixture Model Fitting of Lifespan (Survival) Data
fitmix: Fits lifespan (survival) data with one of the models; Gompertz, Log-logistic, Log-normal, and Weibull.
To install the package fitmix:
require(“devtools”)
install.packages("devtools")
library("devtools")
install_github("guven-code/fitmix")
https://github.com/guven-code/fitmix
