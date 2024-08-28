## bisec


### Overview
-------


The `bisec` package allows you to estimate eigenvalues and eigenvectors and fit a principal component regression (PCR) using a Bayesian approach  under the assumption of a spiked high-dimensional covariance matrix. A key advantage of the proposed models is that they can provide uncertainty estimates. The algorithm is also straightforward to implement, as it requires only inverse-Wishart sampling and spectral decomposition.


### Installation
-------

You can install the package in two ways: 

#### Option 1: Install the most up-to-date version via `devtools`.

To install the most up-to-date version from GitHub, use the following commands:


``` R
if (!("devtools" %in% installed.packages()[,"Package"])) {
    install.packages("devtools")
}
devtools::install_github("swpark0413/bisec")
```

#### Option 2: Install from a bundled package.

First, download a bundled package from the [Releases](https://github.com/swpark0413/bisec/release) page. Then, install it using the command below:

``` R
# Replace "~/path/bisec_0.0.0.tar.gz" with the path to your downloaded file
install.packages("~/path/bisec_0.0.0.tar.gz", type = "source", repos = NULL)
```


### Example
-------

Here is an example of how to use the `bisec` package:

``` r
library(bisec)

# generate a spiked covariance matrix:
n <- 50
p <- 100
K <- 4
leading <- c(1000,500,200,100)
remaining <- rep(0.1, p - K)
Sigma0 <- diag(c(leading, remaining), p)

# generate data
set.seed(413)
X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma0)

# estimate eigenvalues and eigenvectors:
res <- bisec::spikedEIG(X = X, K = 4,prior = list(nu = (p+2), A = diag(2, p)), nsample = 100)
est <- bisec::estimate(res)
```



### Citation
-------

If you use the bisec package in your research, please cite the following paper:

- Kwangmin Lee, Sewon Park, and Jaeyong Lee.
  Bayesian inference on spiked eigenstructure of high-dimensional covariances.
  arXiv preprint arXiv:24xx.xxxxx (2024).

<!-- BibTeX citation:
``` bibtex
@Article{ZhangRD2022gps,
  author        = {Zhang, Ruda and Mak, Simon and Dunson, David},
  title         = {Gaussian Process Subspace Prediction for Model Reduction},
  journal       = {SIAM Journal on Scientific Computing},
  year          = {2022},
  volume        = {44},
  number        = {3},
  pages         = {A1428-A1449},
  doi           = {10.1137/21M1432739},
}
``` -->


### License
-------

The package is distributed under the GPL-3.0 license. See the `LICENSE` file for more details

