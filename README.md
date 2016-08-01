[![Build Status](https://travis-ci.org/Faskally/ef.svg?branch=master)](https://travis-ci.org/Faskally/ef)
[![CRAN Status](http://www.r-pkg.org/badges/version/ef)](https://cran.r-project.org/package=ef)


ef
======

ef implements functions to estimate capture probabilities and densities from multipass pass removal data.

ef is implemented as an [R](https://www.r-project.org) package and available on
[github](https://https://github.com/Faskally/ef).



Installation
------------

the most recent release of ef can be installed from github using the `install_github` command from the `devtools` library:

```R
devtools::install_github('faskally/ef@1.1-0')
```


Usage
-----

For a summary of the package
```R
library(ef)
?ef
```

First lets create some data and fit a basic model, this will introduce the format of the data required and also the basics of model fitting.  The most basic sample is one in which there are 3 passes and 2 lifestages; this constitutes a single sample.

```{r}
# create a single electrofishing site visit with 3 passes and 2 lifestages
ef_data <- data.frame(n      = c(100, 53, 24, 50, 26, 12),
                      pass   = c(  1,  2,  3,  1,  2,  3),
                      stage  = c(  1,  1,  1,  2,  2,  2),
                      sample = c(  1,  1,  1,  2,  2,  2))
ef_data
```

A sensible model to fit to this is one where the capture probability is constant accross passes but varies by lifestage.  This can be acheived by

```{r}
# Fit a simple model
m1 <- efp(n ~ 1 + factor(stage), data = ef_data, pass = pass, id = sample)
m1
```

Notice that it is nessicary to supply the information on which electrifishing pass via the argument 'pass', and the passes have to be linked using the argument 'id'. 




References
----------

The freshwater fisheries laboratory:

[http://www.gov.scot/Topics/marine/Salmon-Trout-Coarse/Freshwater](http://www.gov.scot/Topics/marine/Salmon-Trout-Coarse/Freshwater).


Development
-----------

ef is developed openly on [GitHub](https://github.com/faskally/ef).
Feel free to open an [issue](https://github.com/faskally/ef/issues) there if you encounter problems or have suggestions for future versions.

The current development version can be installed using:

```R
devtools::install_github('faskally/ef')
```
